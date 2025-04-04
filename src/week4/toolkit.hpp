#include <vector>
#include <glm/glm.hpp>
#include <Utils.h>
#include <iostream>
#include <ModelTriangle.h>
#include <map>

const float WIDTH = 640;
const float HEIGHT = 480;

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


void drawTriangle(CanvasTriangle triangle, Colour color, DrawingWindow &window) {
	uint32_t colour = (255 << 24) + (color.red << 16) + (color.green << 8) + int(color.blue);
	line(triangle[0], triangle[1], colour, window);
	line(triangle[1], triangle[2], colour, window);
	line(triangle[2], triangle[0], colour, window);
}

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
}


void lineSampling(CanvasPoint from, CanvasPoint to, TextureMap texture, CanvasTriangle &triangle, DrawingWindow &window) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff)) + 1e-6;
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
	for (float i = 0.0; i <= numberOfSteps; i++) {
		float x = int(from.x + (xStepSize * i) + 0.5);
		float y = int(from.y + (yStepSize * i) + 0.5);
		float rate = i / (float)numberOfSteps;
        auto uv = interpolationUV(triangle, CanvasPoint(x, y));
        // auto uv = tTo + tp * rate;
		auto color = texture.at(uv.x, uv.y);
		window.setPixelColour(x, y, color);
	}
}

void lineSamplingLine(CanvasPoint from, CanvasPoint to, TexturePoint tFrom, TexturePoint tTo, TextureMap texture, CanvasTriangle &triangle, DrawingWindow &window) {
	float xDiff = to.x - from.x;
	float yDiff = to.y - from.y;
	float numberOfSteps = std::max(abs(xDiff), abs(yDiff)) + 1e-6;
	float xStepSize = xDiff/numberOfSteps;
	float yStepSize = yDiff/numberOfSteps;
    auto tpStep = (tTo - tFrom) / numberOfSteps;
    auto tp = tTo - tFrom;
	auto uv = tFrom;
	for (float i = 0.0; i <= numberOfSteps; i++) {
		float x = int(from.x + (xStepSize * i) + 0.5);
		float y = int(from.y + (yStepSize * i) + 0.5);
		// float rate = i / (float)numberOfSteps;
        // auto uv = interpolationUV(triangle, CanvasPoint(x, y));
        // auto uv = tFrom + tp * rate;
		auto color = texture.at(uv.x, uv.y);
		window.setPixelColour(x, y, color);
        uv = uv + tpStep;
	}
}

uint32_t packColor(Colour color) {
	return (255 << 24) + (color.red << 16) + (color.green << 8) + int(color.blue);
}

Colour unpackColor(uint32_t color) {
	Colour c;
	c.red = (color & 0x00FF0000) >> 16;
	c.green = (color & 0x0000FF00) >> 8;
	c.blue = (color & 0x000000FF);

	return c;
}


void drawFillTriangle(CanvasTriangle triangle, Colour color, DrawingWindow &window) {
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

	// drawTriangle(triangle, Colour{255, 255, 255}, window);
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

std::tuple<std::vector<glm::vec3>, std::vector<glm::ivec3>> parseOBJ(std::string fileName) {
	std::ifstream file{fileName, std::ios::in};
	std::vector<glm::vec3> vertics;
	std::vector<glm::ivec3> faces;
	if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << fileName << '\n';
        return {vertics, faces};
    }

	for (std::string str; std::getline(file, str);) {
		if (!str.empty()) {
			if(str.front() == 'v') {
				std::vector<std::string> lineList = split(str, ' ');		
				glm::vec3 vertic{std::stof(lineList[1]), std::stof(lineList[2]), std::stof(lineList[3])};
				vertics.emplace_back(vertic);
			}else if (str.front() == 'f') {
				std::vector<std::string> lineList = split(str, ' ');		
				glm::ivec3 face{std::stoi(lineList[1]), std::stoi(lineList[2]), std::stoi(lineList[3])};
				faces.emplace_back(face);
			}
		}
	}

	return {vertics, faces};
}

std::vector<ModelTriangle> parseOBJ2ModelTriangle(std::string model, std::string mtl) {
	std::ifstream file{model, std::ios::in};
	std::ifstream mtlFile{mtl, std::ios::in};
	std::vector<ModelTriangle> res;
	std::vector<glm::vec3> vertics;
	if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << model << '\n';
        return res;
    }
	if (!mtlFile.is_open())
    {
        std::cerr << "Failed to open file: " << mtl << '\n';
        return res;
    }

	//parse mtl
	std::map<std::string, uint32_t> str2color;
	for (std::string str; std::getline(mtlFile, str);) {
		if (!str.empty()) {
			if(str.find("newmtl") != std::string::npos) {
				auto list = split(str, ' ');	
				auto colorName = list[1];
				std::getline(mtlFile, str);
				list = split(str, ' ');	
				uint32_t color = packColor(Colour{int(std::stof(list[1])*255.0), int(std::stof(list[2])*255.0), int(std::stof(list[3])*255.0)});
				str2color[colorName] = color;
			}
		}
	}

	//parse OBJ
	uint32_t color{0};
	std::string colorName{};
	for (std::string str; std::getline(file, str);) {
		if (!str.empty()) {
			if (str.find("usemtl") != std::string::npos) {
				std::vector<std::string> list = split(str, ' ');		
				colorName = list[1];
				color = str2color[colorName];
			}
			else if(str.front() == 'v') {
				std::vector<std::string> lineList = split(str, ' ');		
				glm::vec3 vertic{std::stof(lineList[1]), std::stof(lineList[2]), std::stof(lineList[3])};
				vertics.emplace_back(vertic);
			}else if (str.front() == 'f') {
				std::vector<std::string> lineList = split(str, ' ');		
				glm::ivec3 face{std::stoi(lineList[1]), std::stoi(lineList[2]), std::stoi(lineList[3])};
				auto c = unpackColor(color);
				c.name = colorName;
				ModelTriangle mt{vertics[face.x-1], vertics[face.y-1], vertics[face.z-1], c};
				res.emplace_back(mt);
			}
		}
	}

	return res;
}



glm::vec3 mapPoint2Screen(const glm::vec3 &p, float width, float height) {
	auto halfWidth = width / 2;
	auto halfHeight = height / 2;
	auto res{p};
	res.x = p.x * halfWidth + halfWidth;
	res.y = p.y * halfHeight + halfHeight;

	return res;
}

glm::vec3 projectVertiexOntoCanvasPoint(glm::vec3 cameraPos, float focalLength, glm::vec3 vertex, float scale = 1.0) {
	vertex = vertex - cameraPos;
	auto u =  focalLength * (vertex.x / abs(vertex.z)) * scale + WIDTH / 2;
	auto v =  focalLength * (vertex.y / vertex.z) * scale + HEIGHT / 2;

	return glm::vec3{u, v, vertex.z};
}