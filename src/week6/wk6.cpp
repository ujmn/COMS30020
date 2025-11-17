#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>
#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h>
#include <cstdlib>
#include <tuple>
#include <TextureMap.h>
#include <toolkit.hpp>
#include <algorithm>
#include <limits>

// #define WIDTH 320
// #define HEIGHT 240

extern const int WIDTH;
extern const int HEIGHT;
extern glm::mat3 cameraOrientation;

enum class DisplayModel : short {
	POINT_CLOUND = 0,
	WIRE_FRAME,
	RASTERISED,
	RAY_TRACED
};

DisplayModel displayModel{DisplayModel::POINT_CLOUND};
glm::vec3 right{1.0, 0.0, 0.0};
glm::vec3 up{0.0, 1.0, 0.0};
glm::vec3 forward{0.0, 0.0, 1.0};
// glm::vec3 lightPos{0.0, 2.0, 0.0};
// glm::vec3 lightPos( -0.64901096, 1.839334, 0.532032);
// glm::vec3 lightPos{0.1, 2.0, 0.1};
glm::vec3 lightPos(0.001f, 2.3389f, 0.007f);

glm::vec3 cameraPos{0.0, 0.0, 15.0};
float focalLength{10.0};
std::vector<std::tuple<CanvasTriangle, Colour>> triangles;
float zbuffer[WIDTH * HEIGHT]{0.0};

bool orbit = false;

void addTriangle() {
	auto v0 = std::rand() % 255;
	auto v1 = std::rand() % 255;
	auto v2 = std::rand() % 255;
	CanvasTriangle t{CanvasPoint{float(std::rand()%255), float(std::rand()%255)}, CanvasPoint{float(std::rand()%255), float(std::rand()%255)}, CanvasPoint{float(std::rand()%255), float(std::rand()%255)}};
	Colour c{std::rand()%255, std::rand()%255, std::rand()%255};

	triangles.emplace_back(std::make_tuple(t, c));
}

uint32_t sampleOnPoint(TexturePoint pt, TextureMap const* tex) {
	uint32_t index = round(pt.y) * tex->width + round(pt.x);
	return tex->pixels[index];
}

uint32_t sampleOnUV(float u, float v, TextureMap const* tex) {
	uint32_t col = round(u * tex->width);
	col %= tex->width;
	uint32_t row = round(v * tex->height);
	row %= tex->height;
	return sampleOnPoint({(float)col, (float)row}, tex);
}

void drawLineSampling(glm::vec2 from, glm::vec2 to, TextureMap const* tex, glm::vec2 uv1, glm::vec2 uv2, DrawingWindow &window) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	glm::vec2 uvDiff = uv2 - uv1;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	auto uvStep = uvDiff / (float)numberOfSteps;
	for(float i = 0.0f; i<numberOfSteps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		auto uv = uv1 + uvStep * i;
		// auto colour = sampleOnUV(uv.x, uv.y, tex);
		auto colour = sampleOnPoint(TexturePoint(uv.x, uv.y), tex);
		window.setPixelColour(round(x), round(y), colour);
	}
}

void drawLineSampling(glm::vec2 from, glm::vec2 to, TextureMap const* tex, TexturePoint tp1, TexturePoint tp2, DrawingWindow &window) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	TexturePoint tpDiff = {tp2.x - tp1.x, tp2.y - tp1.y};
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	TexturePoint tpStep = {tpDiff.x / (float)numberOfSteps, tpDiff.y / (float)numberOfSteps};
	for(float i = 0.0f; i<numberOfSteps; i++) {
		float x = from.x + (xStepSize*i);
		float y = from.y + (yStepSize*i);
		TexturePoint pt = {tp1.x +tpStep.x * i, tp1.y + tpStep.y * i };
		auto colour = sampleOnPoint(pt, tex);
		window.setPixelColour(round(x), round(y), colour);
	}
}


void drawTextureTriangle(CanvasTriangle triangle, TextureMap &text, DrawingWindow &window) {
	if (triangle.v0().y > triangle.v1().y) std::swap(triangle[0], triangle[1]);
    if (triangle.v0().y > triangle.v2().y) std::swap(triangle[0], triangle[2]);
    if (triangle.v1().y > triangle.v2().y) std::swap(triangle[1], triangle[2]);
	auto v0 = triangle.v0();
	auto v1 = triangle.v1();
	auto v2 = triangle.v2();

	auto vm = v0 + (v2-v0) * ((v1.y - v0.y)/(v2.y - v0.y));
	auto tm = v0.texturePoint + (v2.texturePoint-v0.texturePoint) * ((v1.y - v0.y)/(v2.y - v0.y));
	//top
	float xd0m = vm.x - v0.x;
	float xd01 = v1.x - v0.x;
	float yd = vm.y - v0.y;
	float yda = v2.y - v0.y;
	auto td0m = tm - v0.texturePoint;
	auto td02 = v2.texturePoint - v0.texturePoint;
	auto td01 = v1.texturePoint - v0.texturePoint;
	auto uvSteps1 = (v1.texturePoint - v0.texturePoint) / (float)(v1.y - v0.y);
	auto uvSteps2 = (v2.texturePoint - v0.texturePoint) / float(v2.y - v0.y);
	// auto uvFrom = v0.texturePoint;
	// auto uvTo = v0.texturePoint;
	for (int i = 0; i < yd; i++) {
		float rate1 = i/(float)yd;
		float rate2 = i/(float)yda;
		float xFrom = int(v0.x + xd0m * rate1 + 0.5);
		float xTo = int(v0.x + xd01 * rate1 + 0.5);
		// if (xTo < xFrom) { std::swap(xTo, xFrom); }
		// auto uvFrom = td0m * rate1 + v0.texturePoint;
		// auto uvTo = td01 * rate1 + v0.texturePoint;
		auto uvFrom = td02 * rate2 + v0.texturePoint;
		auto uvTo = td01 * rate1 + v0.texturePoint;
		CanvasPoint from{xFrom, v0.y + i};
		CanvasPoint to{xTo, v0.y + i};
		// TexturePoint uvFrom = interpolationUV(triangle, from);
		// TexturePoint uvTo = interpolationUV(triangle, to);
		// lineSampling(from, to, text, triangle, window);
		lineSamplingLine(from, to, uvFrom, uvTo, text, triangle, window);
		// uvTo = uvTo + uvSteps1;
		// uvFrom = uvFrom + uvSteps2;
	}
	//bottom
	yd = v2.y - vm.y;
	float xdm2 = v2.x - vm.x;
	float xd12 = v2.x - v1.x;
	for (int i = 0; i < yd; i++) {
		float rate = i/(float)yd;
		float xFrom = int(vm.x + xdm2 * rate + 0.5);
		float xTo = int(v1.x + xd12 * rate + 0.5);
		if (xTo < xFrom) { std::swap(xTo, xFrom); }
		CanvasPoint from{xFrom, vm.y + i};
		CanvasPoint to{xTo, vm.y + i};
		TexturePoint uvFrom = interpolationUV(triangle, from);
		TexturePoint uvTo = interpolationUV(triangle, to);
		// lineSampling(from, to, text, triangle, window);
		lineSamplingLine(from, to, uvFrom, uvTo, text, triangle, window);
	}

	drawTriangle(triangle, Colour{255, 255, 255}, window);
}



template<typename T>
void drawModel(T model, DrawingWindow &window, float scale = 1.0) {
	auto vertics = std::get<0>(model);
	auto faces = std::get<1>(model);

	for (const auto &face : faces) {
		auto p0 = vertics[face.x - 1] * scale;
		auto p1 = vertics[face.y - 1] * scale;
		auto p2 = vertics[face.z - 1] * scale;
		p0.x += 150;
		p1.x += 150;
		p2.x += 150;
		p0.y += 150;
		p1.y += 150;
		p2.y += 150;
		CanvasTriangle t{p0, p1, p2};
		Colour c{255, 255, 255};
		drawTriangle(t, c, window);
	}
}

void drawWireFrame(DrawingWindow &window, const std::vector<ModelTriangle> &model, float scale = 80.0) {
	for (int i = 0; i < WIDTH * HEIGHT; i++) {
		zbuffer[i] = 0.0;
	}

	for (const auto &triangle : model) {
		auto v0 = triangle.vertices[0];
		v0 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v0, scale);
		auto v1 = triangle.vertices[1];
		v1 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v1, scale);
		auto v2 = triangle.vertices[2];
		v2 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v2, scale);
		CanvasTriangle t{v0, v1, v2};
		drawTriangle(t, triangle.colour, window);
	}
}

void drawFillModel(const std::vector<ModelTriangle> &model, DrawingWindow &window, float scale = 1.0) {
	for (int i = 0; i < WIDTH * HEIGHT; i++) {
		zbuffer[i] = 0.0;
	}

	// cameraPos = cameraToVector * cameraOrientation;
	for (const auto &triangle : model) {
		auto v0 = triangle.vertices[0];
		v0 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v0, scale);
		auto v1 = triangle.vertices[1];
		v1 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v1, scale);
		auto v2 = triangle.vertices[2];
		v2 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v2, scale);
		CanvasTriangle t{v0, v1, v2};

		//test intersection
		// auto intersection = getClosestIntersection(cameraPos, -cameraPos, triangle);
		// if (intersection.distanceFromCamera > 0) {
		// 	std::cout << triangle.colour << std::endl;
		// }
		
		drawFillTriangle(t, triangle.colour, window);
	}
}

void lookAt() {
	glm::vec3 forward = glm::normalize(cameraPos - glm::vec3(0.0, 0.0, 0.0));
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3{0.0, 1.0, 0.0}, forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));

	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
}

void drawPointcloudModel(const std::vector<ModelTriangle> &model, DrawingWindow &window, float scale = 80.0) {
	uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + int(255);
	for (const auto &triangle : model) {
		auto v0 = triangle.vertices[0];
		v0 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v0, scale);
		auto v1 = triangle.vertices[1];
		v1 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v1, scale);
		auto v2 = triangle.vertices[2];
		v2 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, v2, scale);

		// CanvasTriangle t{v0, v1, v2};
		// drawFillTriangle(t, triangle.colour, window);
		window.setPixelColour(v0.x, v0.y, colour);
		window.setPixelColour(v1.x, v1.y, colour);
		window.setPixelColour(v2.x, v2.y, colour);
	}
}

// template<typename T>
// void drawPointcloudModel(T &model, DrawingWindow &window, float scale = 1.0) {
// 	auto vertics = std::get<0>(model);
// 	auto color = packColor(Colour{255, 255, 255});

// 	for (const auto &vertex : vertics) {
// 		auto v = projectVertiexOntoCanvasPoint(cameraPos, focalLength, vertex, scale);
// 		window.setPixelColour(v.x, v.y, color);
// 	}
// }

void setOrbit()
{
	orbit = !orbit;
	// lookAt();
}

template<typename T>
void drawWireframeModel(T &model, DrawingWindow &window, float scale = 1.0) {
	auto [vertics, faces] = model;
	// auto vertics = std::get<0>(model);
	auto color = packColor(Colour{255, 255, 255});

	for (const auto &face : faces) {
		auto p0 = vertics[face.x - 1];
		auto p1 = vertics[face.y - 1];
		auto p2 = vertics[face.z - 1];
		p0 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, p0, scale);
		p1 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, p1, scale);
		p2 = projectVertiexOntoCanvasPoint(cameraPos, focalLength, p2, scale);
		CanvasTriangle t{p0, p1, p2};
		Colour c{255, 255, 255};
		drawTriangle(t, c, window);
	}
}

void changModel() {

}


void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		switch (event.key.keysym.sym)
		{
		case SDLK_LEFT:
			cameraPos.x += 1.0;
			break;
	
		case SDLK_RIGHT:
			cameraPos.x -= 1.0;
			break;

		case SDLK_LSHIFT:
			cameraPos.y -= 1.0;
			break;
	
		case SDLK_SPACE:
			cameraPos.y += 1.0;
			break;

		case SDLK_UP:
			cameraPos.z -= 1.0;
			break;
	
		case SDLK_DOWN:
			cameraPos.z += 1.0;
			break;

		case SDLK_u:
			addTriangle();
			break;
	
		case SDLK_l:
			lookAt();
			break;

		case SDLK_o:
			setOrbit();
			break;
		
		case SDLK_w: 
		{
			auto angle = glm::radians(1.0f);
			auto r = Rotate(rotateAxis::X_AXIS, angle);
			cameraOrientation =  cameraOrientation * r;
			break;
		}
	
		case SDLK_s:
		{
			auto angle = glm::radians(-1.0f);
			auto r = Rotate(rotateAxis::X_AXIS, angle);
			cameraOrientation = cameraOrientation * r;
			break;
		}

		case SDLK_a:
		{
			auto angle = glm::radians(1.0f);
			auto r = Rotate(rotateAxis::Y_AXIS, angle);
			cameraOrientation = r * cameraOrientation;
			break;
		}
	
		case SDLK_d:
		{
			auto angle = glm::radians(-1.0f);
			auto r = Rotate(rotateAxis::Y_AXIS, angle);
			cameraOrientation = r * cameraOrientation;
			break;
		}

		case SDLK_m:
		{
			switch (displayModel)
			{
			case DisplayModel::POINT_CLOUND:
				displayModel = DisplayModel::WIRE_FRAME;
				break;

			case DisplayModel::WIRE_FRAME:
				displayModel = DisplayModel::RASTERISED;
				break;
			
			case DisplayModel::RASTERISED:
				displayModel = DisplayModel::RAY_TRACED;
				break;

			case DisplayModel::RAY_TRACED:
				displayModel = DisplayModel::POINT_CLOUND;
				break;

			default:
				break;
			}
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

void drawRasterScene(DrawingWindow &window, const std::vector<ModelTriangle> &model) {
	window.clearPixels();
	if (orbit) {
		auto angle = -glm::radians(0.2f);
		auto r = Rotate(rotateAxis::Y_AXIS, angle);
		cameraPos = cameraPos * r;
		lookAt();
	}
	drawFillModel(model, window, 80.0);
}

void drawRayTracScene(DrawingWindow &window, std::vector<ModelTriangle> &model) {
	window.clearPixels();
	if (orbit) {
		auto angle = -glm::radians(0.2f);
		auto r = Rotate(rotateAxis::Y_AXIS, angle);
		cameraPos = cameraPos * r;
		lookAt();
	}

	for (int x = 0; x < WIDTH; x++) {
		for (int y = 0; y < HEIGHT; y++) {
			float x_cam = (x - static_cast<float>(WIDTH) / 2.0f) / (focalLength * 80);
    		float y_cam = -(y - static_cast<float>(HEIGHT) / 2.0f) / (focalLength * 80);
    		glm::vec3 dir_camera(x_cam, y_cam, -1.0f);
			dir_camera = glm::normalize(dir_camera);
    		glm::vec3 dir_world = cameraOrientation * dir_camera;
			auto intersection = getClosestIntersection(cameraPos, dir_world, model);
			if (intersection.distanceFromCamera > 0.0) {
				auto color = intersection.intersectedTriangle.colour;
				window.setPixelColour(x, y, packColor(color));
			}
		}
	}
}

void drawRayTracShadowScene(DrawingWindow &window, const std::vector<ModelTriangle> &model) {
	window.clearPixels();
	if (orbit) {
		auto angle = -glm::radians(0.2f);
		auto r = Rotate(rotateAxis::Y_AXIS, angle);
		cameraPos = cameraPos * r;
		lookAt();
	}

	auto black = packColor(Colour{0, 0, 0});

	for (int x = 0; x < WIDTH; x++) {
		for (int y = 0; y < HEIGHT; y++) {
			float x_cam = (x - static_cast<float>(WIDTH) / 2.0f) / (focalLength * 80);
    		float y_cam = -(y - static_cast<float>(HEIGHT) / 2.0f) / (focalLength * 80);
    		glm::vec3 dir_camera(x_cam, y_cam, -1.0f);
			dir_camera = glm::normalize(dir_camera);
    		glm::vec3 dir_world = cameraOrientation * dir_camera;
			auto intersection = getClosestIntersection(cameraPos, dir_world, model);
			if (intersection.distanceFromCamera > 0.0) {
				auto color = intersection.intersectedTriangle.colour;
				if (lightVisible(intersection.intersectionPoint, lightPos, model, intersection.triangleIndex)) {
					window.setPixelColour(x, y, packColor(color));
				}else {
					window.setPixelColour(x, y, black);
				}
			}
		}
	}
}

void draw(DrawingWindow &window, const std::vector<ModelTriangle> &model) {
	window.clearPixels();
	switch (displayModel)
	{
	case DisplayModel::POINT_CLOUND:
		drawPointcloudModel(model, window);
		break;

	case DisplayModel::WIRE_FRAME:
		drawWireFrame(window, model);
		break;
	
	case DisplayModel::RASTERISED:
		drawRasterScene(window, model);
		break;

	case DisplayModel::RAY_TRACED:
		drawRayTracShadowScene(window, model);
		break;
	default:
		break;
	}
}

int main(int argc, char *argv[]) {
	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<CanvasTriangle> triangles;
	auto model = parseOBJ2ModelTriangle("models/cornell-box.obj", "models/cornell-box.mtl");
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		// drawRasterScene(window, model);
		// drawRayTracScene(window, model);
		// drawRayTracShadowScene(window, model);
		draw(window, model);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}

	return 0;
}
