#include <RenderStrategy.h>
#include <Model.h>
#include <RenderUtils.h>
#include <Scene.h>

constexpr float PI = 3.14;

void RasterStrategy::render(const Model &model, std::shared_ptr<RenderContext> renderContext)
{
    const auto &triangles = model.triangles();
    renderContext->clearZbuffer();

    auto window = renderContext->window();
    auto width = window->width;
    auto height = window->height;
    auto cameraOrient = renderContext->cameraOrientation();
	for (const auto &triangle : triangles) {
		auto v0 = triangle.vertices[0];
		v0 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v0, width, height, renderContext->scale());
		auto v1 = triangle.vertices[1];
		v1 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v1, width, height, renderContext->scale());
		auto v2 = triangle.vertices[2];
		v2 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v2, width, height, renderContext->scale());
		CanvasTriangle t{v0, v1, v2};
		renderContext->modelRender()->drawFillTriangle(t, triangle.colour, renderContext);
    }
}

void RayTracingStrategy::render(const Model &model, std::shared_ptr<RenderContext> renderContext)
{
	auto black = RenderUtils::Utils::packColor(Colour{0, 0, 0});

    auto &triangles = model.triangles();
	for (int x = 0; x < renderContext->window()->width; x++) {
		for (int y = 0; y < renderContext->window()->height; y++) {
			float x_cam = (x - static_cast<float>(renderContext->window()->width) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		float y_cam = -(y - static_cast<float>(renderContext->window()->height) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		glm::vec3 dir_camera(x_cam, y_cam, -1.0f);
			dir_camera = glm::normalize(dir_camera);
    		glm::vec3 dir_world = renderContext->cameraOrientation() * dir_camera;
			auto intersection = RenderUtils::Utils::getClosestIntersection(renderContext->cameraPos(), dir_world, triangles);
			if (intersection.distanceFromCamera > 0.0) {
				auto color = intersection.intersectedTriangle.colour;
				if (RenderUtils::Utils::lightVisible(intersection.intersectionPoint, renderContext->lightPos(), triangles, intersection.triangleIndex)) {
					renderContext->window()->setPixelColour(x, y, RenderUtils::Utils::packColor(color));
				}else {
					renderContext->window()->setPixelColour(x, y, black);
				}
			}
		}
	}
}

void PointCloundStrategy::render(const Model &model, std::shared_ptr<RenderContext> renderContext)
{
    auto &triangles = model.triangles();
	uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + int(255);
    auto window = renderContext->window();
    auto width = window->width;
    auto height = window->height;
	for (const auto &triangle : triangles) {
		auto v0 = triangle.vertices[0];
		v0 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v0, width, height, renderContext->scale());
		auto v1 = triangle.vertices[1];
		v1 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v1, width, height, renderContext->scale());
		auto v2 = triangle.vertices[2];
		v2 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v2, width, height, renderContext->scale());

		window->setPixelColour(v0.x, v0.y, colour);
		window->setPixelColour(v1.x, v1.y, colour);
		window->setPixelColour(v2.x, v2.y, colour);
	}
}


void WireFrameStrategy::render(const Model &model, std::shared_ptr<RenderContext> renderContext)
{
    renderContext->clearZbuffer();
    auto triangles = model.triangles();
    auto window = renderContext->window();
    auto width = window->width;
    auto height = window->height;
	for (const auto &triangle : triangles) {
		auto v0 = triangle.vertices[0];
		v0 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v0, width, height, renderContext->scale());
		auto v1 = triangle.vertices[1];
		v1 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v1, width, height, renderContext->scale());
		auto v2 = triangle.vertices[2];
		v2 = RenderUtils::Utils::projectVertiexOntoCanvasPoint(renderContext->cameraPos(), renderContext->cameraOrientation(), renderContext->focalLength(), v2, width, height, renderContext->scale());
		CanvasTriangle t{v0, v1, v2};
		renderContext->modelRender()->drawTriangle(t, triangle.colour, renderContext);
	}
}

void SoftShadowStrategy::render(const Model &model, std::shared_ptr<RenderContext> renderContext)
{
    auto &triangles = model.triangles();
	auto black = RenderUtils::Utils::packColor(Colour{0, 0, 0});
 	auto samplePos = RenderUtils::Utils::sampleDiskOnXZPlane(renderContext->lightPos(), renderContext->sampleRadius(), renderContext->sampleCount());
	for (int x = 0; x < renderContext->window()->width; x++) {
		for (int y = 0; y < renderContext->window()->height; y++) {
			float x_cam = (x - static_cast<float>(renderContext->window()->width) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		float y_cam = -(y - static_cast<float>(renderContext->window()->height) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		glm::vec3 dir_camera(x_cam, y_cam, -1.0f);
			dir_camera = glm::normalize(dir_camera);
    		glm::vec3 dir_world = renderContext->cameraOrientation() * dir_camera;
			auto intersection = RenderUtils::Utils::getClosestIntersection(renderContext->cameraPos(), dir_world, triangles);
			if (intersection.distanceFromCamera > 0.0) {
				auto color = intersection.intersectedTriangle.colour;
				int hitTime = 0;
				for (const auto &pos : samplePos) {
					if (RenderUtils::Utils::lightVisible(intersection.intersectionPoint, pos, triangles, intersection.triangleIndex)) {
						hitTime++;
					}
				}
				float shadingFactor = float(hitTime) / float(renderContext->sampleCount());
				Colour c{int(color.red * shadingFactor), int(color.green * shadingFactor), int(color.blue * shadingFactor)};
				renderContext->window()->setPixelColour(x, y, RenderUtils::Utils::packColor(c));
			}
		}
	}
}

void MirrorReflectStrategy::render(const Model &model, std::shared_ptr<RenderContext> renderContext)
{
    auto window = renderContext->window();
    auto width = window->width;
    auto height = window->height;
    auto &mirrorFaces = model.mirrorFaceIndex();
    auto &triangles = model.triangles();
 	auto samplePos = RenderUtils::Utils::sampleDiskOnXZPlane(renderContext->lightPos(), renderContext->sampleRadius(), renderContext->sampleCount());
    for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			float x_cam = (x - static_cast<float>(width) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		float y_cam = -(y - static_cast<float>(height) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		glm::vec3 dir_camera(x_cam, y_cam, -1.0f);
			dir_camera = glm::normalize(dir_camera);
    		glm::vec3 dir_world = dir_camera * glm::inverse(renderContext->cameraOrientation());
			auto intersection = RenderUtils::Utils::getClosestIntersection(renderContext->cameraPos(), dir_world, triangles);
			if (intersection.distanceFromCamera > 0.0) {
				auto pos = intersection.intersectionPoint;
				auto color = intersection.intersectedTriangle.colour;
				auto dir_light = renderContext->lightPos() - pos; 
				auto dist = glm::length(dir_light);
				auto proximityFactor = std::min(888 / (3 * PI * dist * dist), 1.0f);
				auto ambientFactor = 0.2;

				int hitTime = 0;
				for (const auto &pos : samplePos) {
					if (RenderUtils::Utils::lightVisible(intersection.intersectionPoint, pos, triangles, intersection.triangleIndex)) {
						hitTime++;
					}
				}
				float shadingFactor = float(hitTime) / float(renderContext->sampleCount());

				if (RenderUtils::Utils::lightVisible(intersection.intersectionPoint, renderContext->lightPos(), triangles, intersection.triangleIndex)) {
					auto triangle = intersection.intersectedTriangle;
					auto normal = triangle.normal;
					auto reflectFactor = std::max(glm::dot(normal, glm::normalize(dir_light)), 0.0f);
					auto diffuseFactor = 0.6f * std::min(reflectFactor * proximityFactor, 1.0f);
					auto h = glm::normalize((dir_world + dir_light)/glm::length(dir_light + dir_world));
					auto exp = 16.0f;
					auto specularFactor =  0.2f * std::pow(std::max(glm::dot(normal, h), 0.0f), exp) * proximityFactor;
					auto lightFactor = (diffuseFactor + specularFactor + ambientFactor) * shadingFactor;
					Colour c{int(color.red * lightFactor), int(color.green * lightFactor), int(color.blue * lightFactor)};
					window->setPixelColour(x, y, RenderUtils::Utils::packColor(c));
				}else {
					Colour c{int(color.red * ambientFactor), int(color.green * ambientFactor), int(color.blue * ambientFactor)};
					window->setPixelColour(x, y, RenderUtils::Utils::packColor(c));
				}

                //镜面成像
                if (std::find(mirrorFaces.begin(), mirrorFaces.end(), intersection.triangleIndex) != mirrorFaces.end()) {
				    auto intersectColor = intersection.intersectedTriangle.colour;
				    auto pos = intersection.intersectionPoint;
				    auto I = dir_world; 
				    auto N = intersection.intersectedTriangle.normal;
				    if (glm::dot(I, N) > 0.0f) {
                    	N = -N;
                	}
				    const float EPS = 1e-4f;
                	auto reflectOrigin = pos + EPS * N;
				    auto R = glm::normalize(I - 2* glm::dot(I, N) * N);
				    auto reflectIntersection = RenderUtils::Utils::getClosestIntersection(reflectOrigin, R, triangles);
				    if (reflectIntersection.distanceFromCamera > 0) {
				    	auto pos = reflectIntersection.intersectionPoint;
				    	auto color = reflectIntersection.intersectedTriangle.colour;
				    	color = RenderUtils::Utils::mixColor(color, intersectColor);
				    	auto dir_light = renderContext->lightPos() - pos; 
				    	auto dist = glm::length(dir_light);
				    	auto proximityFactor = std::min(888 / (3 * PI * dist * dist), 1.0f);
				    	int hitTime = 0;
				    	auto ambientFactor = 0.2;
				    	for (const auto &pos : samplePos) {
				    		if (RenderUtils::Utils::lightVisible(reflectIntersection.intersectionPoint, pos, triangles, reflectIntersection.triangleIndex)) {
				    			hitTime++;
				    		}
				    	}
				    	float reflectShadingFactor = float(hitTime) / float(renderContext->sampleCount());

				    	Colour c{int(color.red * ambientFactor), int(color.green * ambientFactor), int(color.blue * ambientFactor)};
				    	if (RenderUtils::Utils::lightVisible(reflectIntersection.intersectionPoint, renderContext->lightPos(), triangles, reflectIntersection.triangleIndex)) {
				    		auto triangle = reflectIntersection.intersectedTriangle;
				    		auto normal = triangle.normal;
				    		auto reflectFactor = std::max(glm::dot(normal, glm::normalize(dir_light)), 0.0f);
				    		auto diffuseFactor = 0.6f * std::min(reflectFactor * proximityFactor, 1.0f);
				    		auto h = glm::normalize((dir_world + dir_light)/glm::length(dir_light + dir_world));
				    		auto exp = 16.0f;
				    		auto specularFactor =  0.2f * std::pow(std::max(glm::dot(normal, h), 0.0f), exp) * proximityFactor;
				    		double lightFactor = 1.0;
				    		if (RenderUtils::Utils::lightVisible(intersection.intersectionPoint, renderContext->lightPos(), triangles, intersection.triangleIndex)) {
				    			lightFactor = (diffuseFactor + specularFactor + ambientFactor) * shadingFactor * reflectShadingFactor;
				    			if (lightFactor < ambientFactor) {
				    				lightFactor = ambientFactor;
				    			}
				    			c = Colour{int(color.red * lightFactor), int(color.green * lightFactor), int(color.blue * lightFactor)};
				    		}
				    	}
				    	window->setPixelColour(x, y, RenderUtils::Utils::packColor(c));
				    }
                }
			}
		}
	}
}

void GouraudShadingStrategy::render(const Model &model, std::shared_ptr<RenderContext> renderContext)
{
	decltype(auto) calculateLight = [](float proximityFactor, glm::vec3 normal, double ambientFactor, glm::vec3 dir_light, glm::vec3 dir_world){
		auto reflectFactor = std::max(glm::dot(normal, glm::normalize(dir_light)), 0.0f);
		auto diffuseFactor = 0.6f * std::min(reflectFactor * proximityFactor, 1.0f);
		auto h = glm::normalize((dir_world + dir_light)/glm::length(dir_light + dir_world));
		auto exp = 16.0f;
		auto specularFactor =  0.2f * std::pow(std::max(glm::dot(normal, h), 0.0f), exp) * proximityFactor;
		auto lightFactor = diffuseFactor + specularFactor + ambientFactor;

		return lightFactor;
	}; 

    auto window = renderContext->window();
    auto width = window->width;
    auto height = window->height;
    auto &mirrorFaces = model.mirrorFaceIndex();
    auto &triangles = model.triangles();
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			float x_cam = (x - static_cast<float>(width) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		float y_cam = -(y - static_cast<float>(height) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		glm::vec3 dir_camera(x_cam, y_cam, -1.0f);
			dir_camera = glm::normalize(dir_camera);
    		glm::vec3 dir_world = dir_camera * glm::inverse(renderContext->cameraOrientation());
			auto intersection = RenderUtils::Utils::getClosestIntersection(renderContext->cameraPos(), dir_world, triangles);
			if (intersection.distanceFromCamera > 0.0) {
				auto pos = intersection.intersectionPoint;
				auto color = intersection.intersectedTriangle.colour;
				auto dir_light = renderContext->lightPos() - pos; 
				auto dist = glm::length(dir_light);
				auto proximityFactor = std::min(888 / (3 * PI * dist * dist), 1.0f);
				auto ambientFactor = 0.2;
				if (RenderUtils::Utils::lightVisible(intersection.intersectionPoint, renderContext->lightPos(), triangles, intersection.triangleIndex)) {
					auto triangle = intersection.intersectedTriangle;
					auto normal = triangle.normal;
					auto triangleIndex = intersection.triangleIndex;
					auto lightFactor = 1.0;
					Colour color_w{0, 0, 0};
					if (triangleIndex >= 0) {
						auto verList = model.faceIndex2verInd().at(triangleIndex);
						auto bcCoord = RenderUtils::Utils::TDBarycentric(triangle.vertices[0], triangle.vertices[1], triangle.vertices[2], intersection.intersectionPoint);
							auto lightFactor_1 = calculateLight(proximityFactor, model.vertexNormal()[verList[0]], ambientFactor, dir_light, dir_world);
							auto lightFactor_2 = calculateLight(proximityFactor, model.vertexNormal()[verList[1]], ambientFactor, dir_light, dir_world);
							auto lightFactor_3 = calculateLight(proximityFactor, model.vertexNormal()[verList[2]], ambientFactor, dir_light, dir_world);
							Colour c_1{int(color.red * lightFactor_1), int(color.green * lightFactor_1), int(color.blue * lightFactor_1)};
							Colour c_2{int(color.red * lightFactor_2), int(color.green * lightFactor_2), int(color.blue * lightFactor_2)};
							Colour c_3{int(color.red * lightFactor_3), int(color.green * lightFactor_3), int(color.blue * lightFactor_3)};
							color_w = c_1 * bcCoord[0] + c_2 * bcCoord[1] + c_3 * bcCoord[2];
					}
					window->setPixelColour(x, y, RenderUtils::Utils::packColor(color_w));
				}else {
					Colour c{int(color.red * ambientFactor), int(color.green * ambientFactor), int(color.blue * ambientFactor)};
					window->setPixelColour(x, y, RenderUtils::Utils::packColor(c));
				}
			}
		}
	}
}

void PhongShadingStrategy::render(const Model &model, std::shared_ptr<RenderContext> renderContext)
{
    auto window = renderContext->window();
    auto width = window->width;
    auto height = window->height;
    auto &mirrorFaces = model.mirrorFaceIndex();
    auto &triangles = model.triangles();
	for (int x = 0; x < width; x++) {
		for (int y = 0; y < height; y++) {
			float x_cam = (x - static_cast<float>(width) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		float y_cam = -(y - static_cast<float>(height) / 2.0f) / (renderContext->focalLength() * renderContext->scale());
    		glm::vec3 dir_camera(x_cam, y_cam, -1.0f);
			dir_camera = glm::normalize(dir_camera);
    		glm::vec3 dir_world = dir_camera * glm::inverse(renderContext->cameraOrientation());
			auto intersection = RenderUtils::Utils::getClosestIntersection(renderContext->cameraPos(), dir_world, triangles);
			if (intersection.distanceFromCamera > 0.0) {
				auto pos = intersection.intersectionPoint;
				auto color = intersection.intersectedTriangle.colour;
				auto dir_light = renderContext->lightPos() - pos; 
				auto dist = glm::length(dir_light);
				auto proximityFactor = std::min(888 / (3 * PI * dist * dist), 1.0f);
				auto ambientFactor = 0.2;
				if (RenderUtils::Utils::lightVisible(intersection.intersectionPoint, renderContext->lightPos(), triangles, intersection.triangleIndex)) {
					auto triangle = intersection.intersectedTriangle;
					auto normal = triangle.normal;
					//计算向量插值
					auto triangleIndex = intersection.triangleIndex;
					if (triangleIndex >= 0) {
						auto verList = model.faceIndex2verInd()[triangleIndex];
						auto bcCoord = RenderUtils::Utils::TDBarycentric(triangle.vertices[0], triangle.vertices[1], triangle.vertices[2], intersection.intersectionPoint);
						normal = glm::normalize(model.vertexNormal()[verList[0]] * bcCoord[0] + model.vertexNormal()[verList[1]] * bcCoord[1] + model.vertexNormal()[verList[2]] * bcCoord[2]);
					}
					auto reflectFactor = std::max(glm::dot(normal, glm::normalize(dir_light)), 0.0f);
					auto diffuseFactor = 0.6f * std::min(reflectFactor * proximityFactor, 1.0f);
					auto h = glm::normalize((dir_world + dir_light)/glm::length(dir_light + dir_world));
					auto exp = 16.0f;
					auto specularFactor =  0.2f * std::pow(std::max(glm::dot(normal, h), 0.0f), exp) * proximityFactor;
					auto lightFactor = diffuseFactor + specularFactor + ambientFactor;
					Colour c{int(color.red * lightFactor), int(color.green * lightFactor), int(color.blue * lightFactor)};
					window->setPixelColour(x, y, RenderUtils::Utils::packColor(c));
				}else {
					Colour c{int(color.red * ambientFactor), int(color.green * ambientFactor), int(color.blue * ambientFactor)};
					window->setPixelColour(x, y, RenderUtils::Utils::packColor(c));
				}
			}
		}
	}
}
