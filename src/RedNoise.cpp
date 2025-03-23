#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <Utils.h>
#include <fstream>
#include <vector>
#include <glm/glm.hpp>

#define WIDTH 320
#define HEIGHT 240

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
void draw(DrawingWindow &window) {
	window.clearPixels();
	glm::vec3 topLeft(255, 0, 0);        // red 
	glm::vec3 topRight(0, 0, 255);       // blue 
	glm::vec3 bottomRight(0, 255, 0);    // green 
	glm::vec3 bottomLeft(255, 255, 0);   // yellow
	std::vector<glm::vec3> left = interpolateThreeElementValues(topLeft, bottomLeft, HEIGHT);
	std::vector<glm::vec3> right = interpolateThreeElementValues(topRight, bottomRight, HEIGHT);
	for (int i = 0; i < HEIGHT; i++) {
		float xNum = (right[i].x - left[i].x) / (float)(WIDTH - 1);
		float yNum = (right[i].y - left[i].y) / (float)(WIDTH - 1);
		float zNum = (right[i].z - left[i].z) / (float)(WIDTH - 1);
		glm::vec3 tmp = left[i];
		for (int j = 0; j < WIDTH; j++) {
			float red = tmp.x;
			float green = tmp.y;
			float blue = tmp.z;
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(j, i, colour);
			tmp.x += xNum;
			tmp.y += yNum;
			tmp.z += zNum;
		}
	}
}

void handleEvent(SDL_Event event, DrawingWindow &window) {
	if (event.type == SDL_KEYDOWN) {
		if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
		else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
		else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
		else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
	} else if (event.type == SDL_MOUSEBUTTONDOWN) {
		window.savePPM("output.ppm");
		window.saveBMP("output.bmp");
	}
}

int main(int argc, char *argv[]) {
	DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);
	SDL_Event event;
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
