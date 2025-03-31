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

#define WIDTH 320
#define HEIGHT 240

std::vector<std::tuple<CanvasTriangle, Colour>> triangles;

std::vector<float> interpolateSingleFloats(float from, float to, int numberOfValues) {
	std::vector<float> ret;
	float num = (to - from) / (float)(numberOfValues - 1);
	float tmp = from;
	for (int i = 0; i < numberOfValues; i++) {
		ret.emplace_back(tmp);
		tmp += num;
	}

	return ret;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, int numberOfValues) {
	std::vector<glm::vec3> ret;
	float xNum = (to.x - from.x) / (float)(numberOfValues - 1);
	float yNum = (to.y - from.y) / (float)(numberOfValues - 1);
	float zNum = (to.z - from.z) / (float)(numberOfValues - 1);
	glm::vec3 tmp = from;
	for (int i = 0; i < numberOfValues; i++) {
		ret.emplace_back(tmp);
		tmp.x += xNum;
		tmp.y += yNum;
		tmp.z += zNum;
	}

	return ret;
}

void line(CanvasPoint from, CanvasPoint to, uint32_t color, DrawingWindow &window) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff));
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	for (float i = 0.0; i <= numberOfSteps; i++) {
		float x = int(from.x + (xStepSize * i) + 0.5);
		float y = int(from.y + (yStepSize * i) + 0.5);

		window.setPixelColour(x, y, color);
	}
}

void addTriangle() {
	auto v0 = std::rand() % 255;
	auto v1 = std::rand() % 255;
	auto v2 = std::rand() % 255;
	CanvasTriangle t{CanvasPoint{float(std::rand()%255), float(std::rand()%255)}, CanvasPoint{float(std::rand()%255), float(std::rand()%255)}, CanvasPoint{float(std::rand()%255), float(std::rand()%255)}};
	Colour c{std::rand()%255, std::rand()%255, std::rand()%255};

	triangles.emplace_back(std::make_tuple(t, c));
}

void drawTriangle(CanvasTriangle triangle, Colour color, DrawingWindow &window) {
	uint32_t colour = (255 << 24) + (color.red << 16) + (color.green << 8) + int(color.blue);
	line(triangle[0], triangle[1], colour, window);
	line(triangle[1], triangle[2], colour, window);
	line(triangle[2], triangle[0], colour, window);
}


void fillTriangle(CanvasTriangle triangle, Colour color, DrawingWindow &window) {
	auto v0 = triangle.v0();
	auto v1 = triangle.v1();
	auto v2 = triangle.v2();
	if (v0.y > v1.y) { std::swap(v0, v1); }
	if (v0.y > v2.y) { std::swap(v0, v2); }
	if (v1.y > v2.y) { std::swap(v1, v2); }
	auto vm = v0 + (v2-v0) * ((v1.y - v0.y)/(v2.y - v0.y));
	// std::cout << vm.x << ", " << vm.y << std::endl;
	// uint32_t colour = (255 << 24) + (255 << 16) + (255 << 8) + int(0);
	// window.setPixelColour(vm.x + 3, vm.y, colour);
	//top
	uint32_t colour = (255 << 24) + (color.red << 16) + (color.green << 8) + int(color.blue);
	float xd0m = vm.x - v0.x;
	float xd01 = v1.x - v0.x;
	float yd = vm.y - v0.y;
	for (int i = 0; i < yd; i++) {
		float rate = i/(float)yd;
		float xFrom = int(v0.x + xd0m * rate + 0.5);
		float xTo = int(v0.x + xd01 * rate + 0.5);
		line(CanvasPoint{xFrom, v0.y+i}, CanvasPoint{xTo, v0.y+i}, colour, window);
	}
	//bottom
	yd = v2.y - vm.y;
	float xdm2 = v2.x - vm.x;
	float xd12 = v2.x - v1.x;
	for (int i = 0; i < yd; i++) {
		float rate = i/(float)yd;
		float xFrom = int(vm.x + xdm2 * rate + 0.5);
		float xTo = int(v1.x + xd12 * rate + 0.5);
		line(CanvasPoint{xFrom, vm.y+i}, CanvasPoint{xTo, vm.y+i}, colour, window);
	}

	drawTriangle(triangle, Colour{255, 255, 255}, window);
}
std::vector<float> barycentric(const CanvasPoint &A, const CanvasPoint &B, const CanvasPoint &C,
                                     const CanvasPoint &P) {
	std::vector<float> res;
    // 提取顶点坐标
    float x1 = A.x, y1 = A.y;
    float x2 = B.x, y2 = B.y;
    float x3 = C.x, y3 = C.y;
    float x = P.x, y = P.y;

    // 计算分母（总面积的两倍）
    float denominator = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3);

    // 如果分母为零，说明三点共线或重合，无法形成三角形
    if (std::abs(denominator) < 1e-6) {
        std::cerr << "Error: The points are collinear or coincident!" << std::endl;
		res.emplace_back(0.0);
		res.emplace_back(0.0);
		res.emplace_back(0.0);
        return res;
    }

    // 计算 λ1 和 λ2
    auto lambda1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denominator;
    auto lambda2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denominator;

    // 计算 λ3
    auto lambda3 = 1.0f - lambda1 - lambda2;
	res.emplace_back(lambda1);
	res.emplace_back(lambda2);
	res.emplace_back(lambda3);

	return res;
}

// TexturePoint interpolationUV(CanvasTriangle &triangle, CanvasPoint &p) {
// 	std::vector<float> baryCoord = barycentric(triangle.v0(), triangle.v1(), triangle.v2(), p);
// 	return triangle.v0().texturePoint * baryCoord[0] + triangle.v1().texturePoint * baryCoord[1] + triangle.v2().texturePoint * baryCoord[2];
// }

TexturePoint interpolationUV(CanvasTriangle &triangle, CanvasPoint &p) {
    CanvasPoint &v0 = triangle[0];
    CanvasPoint &v1 = triangle[1];
    CanvasPoint &v2 = triangle[2];

    float denominator = (v1.y - v2.y) * (v0.x - v2.x) + (v2.x - v1.x) * (v0.y - v2.y);
    if (denominator == 0.0f) {
        return TexturePoint(0.0f, 0.0f);
    }

    float alpha_numerator = (v1.y - v2.y) * (p.x - v2.x) + (v2.x - v1.x) * (p.y - v2.y);
    float alpha = alpha_numerator / denominator;

    float beta_numerator = (v2.y - v0.y) * (p.x - v2.x) + (v0.x - v2.x) * (p.y - v2.y);
    float beta = beta_numerator / denominator;

    float gamma = 1.0f - alpha - beta;

    float u = alpha * v0.texturePoint.x + beta * v1.texturePoint.x + gamma * v2.texturePoint.x;
    float v = alpha * v0.texturePoint.y + beta * v1.texturePoint.y + gamma * v2.texturePoint.y;

    return TexturePoint(u, v);
    // return TexturePoint(alpha, beta);
}

void drawTextureTriangle(CanvasTriangle triangle, TextureMap &text, DrawingWindow &window) {
	if (triangle.v0().y > triangle.v1().y) std::swap(triangle[0], triangle[1]);
    if (triangle.v0().y > triangle.v2().y) std::swap(triangle[0], triangle[2]);
    if (triangle.v1().y > triangle.v2().y) std::swap(triangle[1], triangle[2]);
	auto v0 = triangle.v0();
	auto v1 = triangle.v1();
	auto v2 = triangle.v2();
	// if (v0.y > v1.y) { std::swap(v0, v1); }
	// if (v0.y > v2.y) { std::swap(v0, v2); }
	// if (v1.y > v2.y) { std::swap(v1, v2); }

	auto vm = v0 + (v2-v0) * ((v1.y - v0.y)/(v2.y - v0.y));
	//top
	float xd0m = vm.x - v0.x;
	float xd01 = v1.x - v0.x;
	float yd = vm.y - v0.y;
	for (int i = 0; i < yd; i++) {
		float rate = i/(float)yd;
		float xFrom = int(v0.x + xd0m * rate + 0.5);
		float xTo = int(v0.x + xd01 * rate + 0.5);
		if (xTo < xFrom) { std::swap(xTo, xFrom); }
		float xDiff = xTo - xFrom;
		float numberOfSteps = abs(xDiff) + 1e-6;
		for (float j = 0.0; j <= numberOfSteps; j++) {
			float x = int(xFrom + j + 0.5);
			float rate = j / numberOfSteps;
			// auto uv = uvTo + uvd * rate;
			auto uv = interpolationUV(triangle, CanvasPoint(x, v0.y + i));
			auto color = text.at(uv.x, uv.y);
			// int red = uv.x * 255.0;
			// int green = uv.y * 255.0;
			// int blue = (1 - uv.x - uv.y) * 255.0;
			// auto color = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue); 
			window.setPixelColour(x, v0.y + i, color);
		}
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
		float xDiff = xTo - xFrom;
		float numberOfSteps = abs(xDiff) + 1e-6;
		for (float j = 0.0; j <= numberOfSteps; j++) {
			float x = int(xFrom + j + 0.5);
			auto uv = interpolationUV(triangle, CanvasPoint(x, vm.y + i));
			auto color = text.at(uv.x, uv.y);
			// int red = uv.x * 255.0;
			// int green = uv.y * 255.0;
			// int blue = (1 - uv.x - uv.y) * 255.0;
			// auto color = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue); 
			window.setPixelColour(x, vm.y + i, color);
		}
	}

	drawTriangle(triangle, Colour{255, 255, 255}, window);
}

//task3:single dimension greyscale interpolation
// void draw(DrawingWindow &window) {
// 	window.clearPixels();
// 	std::vector<float> row = interpolateSingleFloats(255.0, 0.0, WIDTH);
// 	std::cout << row[0] << std::endl;
// 	for (size_t y = 0; y < window.height; y++) {
// 		for (size_t x = 0; x < window.width; x++) {
// 			float red = row[x];
// 			float green = row[x];
// 			float blue = row[x];
// 			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
// 			window.setPixelColour(x, y, colour);
// 		}
// 	}
// }

//task5 : two dimensional colour interpolation
// void draw(DrawingWindow &window) {
// 	window.clearPixels();
// 	glm::vec3 topLeft(255, 0, 0);        // red 
// 	glm::vec3 topRight(0, 0, 255);       // blue 
// 	glm::vec3 bottomRight(0, 255, 0);    // green 
// 	glm::vec3 bottomLeft(255, 255, 0);   // yellow
// 	std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft, bottomLeft, HEIGHT);
// 	std::vector<glm::vec3> right = interpolateThreeElementValues(topRight, bottomRight, HEIGHT);
// 	for (int i = 0; i < HEIGHT; i++) {
// 		float xNum = (right[i].x - left[i].x) / (float)(WIDTH - 1);
// 		float yNum = (right[i].y - left[i].y) / (float)(WIDTH - 1);
// 		float zNum = (right[i].z - left[i].z) / (float)(WIDTH - 1);
// 		glm::vec3 tmp = left[i];
// 		for (int j = 0; j < WIDTH; j++) {
// 			float red = tmp.x;
// 			float green = tmp.y;
// 			float blue = tmp.z;
// 			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
// 			window.setPixelColour(j, i, colour);
// 			tmp.x += xNum;
// 			tmp.y += yNum;
// 			tmp.z += zNum;
// 		}
// 	}
// }

void draw(DrawingWindow &window) {
	window.clearPixels();
	CanvasPoint v0{160, 10};
	v0.texturePoint = TexturePoint{195, 5};
	CanvasPoint v2{300, 230};
	v2.texturePoint = TexturePoint{395, 380};
	CanvasPoint v1{10, 150};
	v1.texturePoint = TexturePoint{65, 330};
	CanvasTriangle triangle{v0, v1, v2};
	TextureMap text("./texture.ppm");
	drawTextureTriangle(triangle, text, window);

	// for (auto triangle : triangles) {
	// 	CanvasTriangle t = std::get<0>(triangle);
	// 	Colour c = std::get<1>(triangle);
	// 	// drawTriangle(t, c, window);
	// 	fillTriangle(t, c, window);
	// }
	// uint32_t colour = (255 << 24) + (int(255) << 16) + (int(255) << 8) + int(255);
	// line(CanvasPoint(0.0, 0.0), CanvasPoint(WIDTH/2, HEIGHT/2), colour, window);
	// line(CanvasPoint(WIDTH/2, HEIGHT/2), CanvasPoint(WIDTH, 0.0), colour, window);
	// line(CanvasPoint(WIDTH/2, 0.0), CanvasPoint(WIDTH/2, HEIGHT), colour, window);
	// line(CanvasPoint(WIDTH/2 - 100.0, HEIGHT/2), CanvasPoint(WIDTH/2 + 100.0, HEIGHT/2), colour, window);
}


void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
		else if (event.key.keysym.sym == SDLK_u) {	addTriangle(); std::cout << "u" << std::endl;	}
		else if (event.key.keysym.sym == SDLK_f) {	addTriangle(); std::cout << "f" << std::endl;	}
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
	std::vector<CanvasTriangle> triangles;
	while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
	//test interpolateSingleFloats
	// std::vector<float> result;
	// result = interpolateSingleFloats(2.2, 8.5, 7);
	// for(size_t i=0; i<result.size(); i++) std::cout << result[i] << " ";
	// std::cout << std::endl;


	//test interpolateThreeElementValues
	// glm::vec3 from(1.0, 4.0, 9.2);
	// glm::vec3 to(4.0, 1.0, 9.8);
	// std::vector<glm::vec3> result = interpolateThreeElementValues(from, to, 4);
	// for (const auto &point : result) {
	// 	std::cout << point.x << ", " << point.y << ", " << point.z << std::endl;
	// }

	return 0;
}
