#include <Scene.h>
#include <RenderContext.h>
#include <RenderUtils.h>
#include <CanvasTriangle.h>
#include <RenderStrategy.h>
#include <limits>

Scene::Scene(int width, int height, glm::vec3 cameraPos, float focalLength)
    : renderContext_(std::make_shared<RenderContext>(RenderContext{cameraPos, focalLength, width, height}))
{
	renderStrategyFactory_ = std::make_unique<PointCloundStrategyFactor>();
	renderStrategyList_.emplace_back(std::move(renderStrategyFactory_->CreateStrategy()));

	renderStrategyFactory_ = std::make_unique<WireFrameStrategyFactor>();
	renderStrategyList_.emplace_back(std::move(renderStrategyFactory_->CreateStrategy()));

	renderStrategyFactory_ = std::make_unique<RasterStrategyFactor>();
	renderStrategyList_.emplace_back(std::move(renderStrategyFactory_->CreateStrategy()));

	renderStrategyFactory_ = std::make_unique<RayTracingStrategyFactor>();
	renderStrategyList_.emplace_back(std::move(renderStrategyFactory_->CreateStrategy()));

	renderStrategyFactory_ = std::make_unique<SoftShadowStrategyFactor>();
	renderStrategyList_.emplace_back(std::move(renderStrategyFactory_->CreateStrategy()));

	renderStrategyFactory_ = std::make_unique<MirrorReflectStrategyFactor>();
	renderStrategyList_.emplace_back(std::move(renderStrategyFactory_->CreateStrategy()));

	renderStrategyFactory_ = std::make_unique<GouraudShadingStrategyFactor>();
	renderStrategyList_.emplace_back(std::move(renderStrategyFactory_->CreateStrategy()));

	renderStrategyFactory_ = std::make_unique<PhongShadingStrategyFactor>();
	renderStrategyList_.emplace_back(std::move(renderStrategyFactory_->CreateStrategy()));
}

void Scene::update()
{
	renderContext_->window()->clearPixels();
    processOrbit();
	if (renderContext_->window()->pollForInputEvents(event_)){
		handleEvent(event_, *renderContext_->window_);
	} 

    for (const auto &model : models_) {
		if (renderStrategyList_.size() > 0) {
        	renderStrategyList_.at(renderStrategyIndex_)->render(model, renderContext_);
		}
    }
	renderContext_->window()->renderFrame();
}

void Scene::processOrbit()
{
	if (renderContext_->orbit()) {
		auto angle = -glm::radians(0.2f);
		auto r = RenderUtils::Utils::Rotate(RenderUtils::rotateAxis::Y_AXIS, angle);
		renderContext_->cameraPos_ = renderContext_->cameraPos_ * r;
		renderContext_->setCameraOrientation(RenderUtils::Utils::lookAt(renderContext_->cameraPos()));
	}
}

void Scene::addModel(Model model)
{
    models_.emplace_back(model);
}

Scene &Scene::getScene(int width, int height, glm::vec3 cameraPos, float focalLength)
{
    static Scene scene(width, height, cameraPos, focalLength);
    return scene;
}

void Scene::addRenderStrategy(std::unique_ptr<IRenderStrategy> renderStartegy)
{
	renderStrategyList_.emplace_back(std::move(renderStartegy));
}

void Scene::start()
{
    while(true) {
        update();
    }
}

void ModelRender::drawFillTriangle(CanvasTriangle &triangle, Colour color, std::shared_ptr<RenderContext> context) {
	auto v0 = triangle.v0();
	auto v1 = triangle.v1();
	auto v2 = triangle.v2();
	if (v0.y > v1.y) { std::swap(v0, v1); }
	if (v0.y > v2.y) { std::swap(v0, v2); }
	if (v1.y > v2.y) { std::swap(v1, v2); }
	auto vm = v0 + (v2-v0) * ((v1.y - v0.y)/(v2.y - v0.y));
	//top
	uint32_t colour = (255 << 24) + (color.red << 16) + (color.green << 8) + int(color.blue);
	float xd0m = vm.x - v0.x;
	float xd01 = v1.x - v0.x;
	float zd0m = vm.depth - v0.depth;
	float zd01 = v1.depth - v0.depth;
	float yd = vm.y - v0.y;
	for (int i = 0; i < yd; i++) {
		float rate = i/(float)yd;
		float xFrom = int(v0.x + xd0m * rate + 0.5);
		float xTo = int(v0.x + xd01 * rate + 0.5);
		float zFrom = v0.depth + zd0m * rate;
		float zTo = v0.depth + zd01 * rate;
		line(CanvasPoint{xFrom, v0.y+i, zFrom}, CanvasPoint{xTo, v0.y+i, zTo}, colour, context);
	}
	//bottom
	yd = v2.y - vm.y;
	float xdm2 = v2.x - vm.x;
	float xd12 = v2.x - v1.x;
	float zdm2 = v2.depth - vm.depth;
	float zd12 = v2.depth - v1.depth;
	for (int i = 0; i < yd; i++) {
		float rate = i/(float)yd;
		float xFrom = int(vm.x + xdm2 * rate + 0.5);
		float xTo = int(v1.x + xd12 * rate + 0.5);
		float zFrom = vm.depth + zdm2 * rate;
		float zTo = v1.depth + zd12 * rate;
		line(CanvasPoint{xFrom, vm.y+i, zFrom}, CanvasPoint{xTo, vm.y+i, zTo}, colour, context);
	}
}

void ModelRender::drawTriangle(CanvasTriangle &triangle, Colour color, std::shared_ptr<RenderContext> context)
{
	uint32_t colour = (255 << 24) + (color.red << 16) + (color.green << 8) + int(color.blue);
	line(triangle[0], triangle[1], colour, context);
	line(triangle[1], triangle[2], colour, context);
	line(triangle[2], triangle[0], colour, context);
}

void ModelRender::line(CanvasPoint &from, CanvasPoint &to, uint32_t color, std::shared_ptr<RenderContext> context) {
    auto window = context->window();
    auto &zbuffer = context->zbuffer();

	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float zDiff = to.depth - from.depth;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/(numberOfSteps + 1e-6);
	float yStepSize = yDiff/(numberOfSteps + 1e-6);
	float zStepSize = zDiff/(numberOfSteps + 1e-6);
	for (float i = 0.0; i <= numberOfSteps; i++) {
		int x = int(from.x + (xStepSize * i) + 0.5);
		int y = int(from.y + (yStepSize * i) + 0.5);
		float z = from.depth + (zStepSize * i);

		if ( (y-1) * window->height + x >= 0 && (y-1) * window->height + x < window->width * window->height && zbuffer[(y-1) * window->height + x] < z) {
			zbuffer[(y-1) * window->height + x] = z;
			context->window()->setPixelColour(x, y, color);
		} 
	}
}

void Scene::handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		switch (event.key.keysym.sym)
		{
		case SDLK_LEFT:
			renderContext_->cameraPos_.x += 1.0;
			break;
	
		case SDLK_RIGHT:
			renderContext_->cameraPos_.x -= 1.0;
			break;

		case SDLK_LSHIFT:
			renderContext_->cameraPos_.y -= 1.0;
			break;
	
		case SDLK_SPACE:
			renderContext_->cameraPos_.y += 1.0;
			break;

		case SDLK_UP:
			renderContext_->cameraPos_.z -= 1.0;
			break;
	
		case SDLK_DOWN:
			renderContext_->cameraPos_.z += 1.0;
			break;

		case SDLK_u:
			// addTriangle();
			break;
	
		case SDLK_l:
			RenderUtils::Utils::lookAt(renderContext_->cameraPos_);
			break;

		case SDLK_o:
			setOrbit();
			break;
		
		case SDLK_w: 
		{
			auto angle = glm::radians(1.0f);
			auto r = RenderUtils::Utils::Rotate(RenderUtils::rotateAxis::X_AXIS, angle);
			renderContext_->cameraOrientation_ =  renderContext_->cameraOrientation_ * r;
			break;
		}
	
		case SDLK_s:
		{
			auto angle = glm::radians(-1.0f);
			auto r = RenderUtils::Utils::Rotate(RenderUtils::rotateAxis::X_AXIS, angle);
			renderContext_->cameraOrientation_ =  renderContext_->cameraOrientation_ * r;
			break;
		}

		case SDLK_a:
		{
			auto angle = glm::radians(1.0f);
			auto r = RenderUtils::Utils::Rotate(RenderUtils::rotateAxis::Y_AXIS, angle);
			renderContext_->cameraOrientation_ =  renderContext_->cameraOrientation_ * r;
			break;
		}
	
		case SDLK_d:
		{
			auto angle = glm::radians(-1.0f);
			auto r = RenderUtils::Utils::Rotate(RenderUtils::rotateAxis::Y_AXIS, angle);
			renderContext_->cameraOrientation_ =  renderContext_->cameraOrientation_ * r;
			break;
		}

		case SDLK_m:
		{
			renderStrategyIndex_ = (renderStrategyIndex_ + 1) % renderStrategyList_.size();
			break;
		}

		default:
			break;
		}
	}else if (event.type == SDL_MOUSEBUTTONDOWN){
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

void Scene::setOrbit()
{
	renderContext_->orbit_ = !renderContext_->orbit_;
}

void Scene::setLightPos(glm::vec3 lightPos)
{
	renderContext_->lightPos_ = lightPos;
}
