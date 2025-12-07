#include <vector>
#include <glm/glm.hpp>
#include <Utils.h>
#include <iostream>
#include <ModelTriangle.h>
#include <map>
#include <CanvasTriangle.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <cassert>
#include <limits>
#include <cmath>
#include <algorithm>


// const int WIDTH = 1280;
// const int HEIGHT = 960;
const int WIDTH = 640;
const int HEIGHT = 480;
glm::mat3 cameraOrientation;
extern float zbuffer[WIDTH * HEIGHT];

enum class rotateAxis{
	X_AXIS,
	Y_AXIS,
	Z_AXIS
};

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

// 同心映射：将 [0,1)^2 -> 单位圆盘 (返回 {dx, dz})
glm::vec3 squareToDiskConcentricXZ(double u, double v) {
    double ux = 2.0 * u - 1.0; // [-1, 1)
    double uz = 2.0 * v - 1.0;

    if (ux == 0.0 && uz == 0.0) {
        return glm::vec3(0.0, 0.0, 0.0); // y 分量无意义，设为 0
    }

    double r, theta;
    if (std::abs(ux) > std::abs(uz)) {
        r = ux;
        theta = (M_PI / 4.0) * (uz / ux);
    } else {
        r = uz;
        theta = (M_PI / 2.0) - (M_PI / 4.0) * (ux / uz);
    }

    double dx = r * std::cos(theta);
    double dz = r * std::sin(theta);
    return glm::vec3(dx, 0.0, dz); // y = 0，仅表示偏移
}

// 在 XZ 平面上的圆盘内确定性采样 N 个点
std::vector<glm::vec3> sampleDiskOnXZPlane(
    const glm::vec3& center, 
    double radius, 
    int numSamples
) {
    if (numSamples <= 0) return {};

    int n = static_cast<int>(std::ceil(std::sqrt(numSamples)));
    std::vector<glm::vec3> samples;
    samples.reserve(numSamples);

    for (int i = 0; i < n * n && samples.size() < numSamples; ++i) {
        int ix = i % n;
        int iz = i / n;

        double u = (ix + 0.5) / n; // [0,1)
        double v = (iz + 0.5) / n;

        // 映射到单位圆盘（XZ 偏移）
        glm::vec3 offset = squareToDiskConcentricXZ(u, v);

        // 缩放并平移到目标位置
        glm::vec3 point;
        point.x = center.x + radius * offset.x;
        point.y = center.y;
        point.z = center.z + radius * offset.z;

        samples.push_back(point);
    }

    return samples;
}


Colour operator*(Colour color, float num) {
	return Colour{int(color.red * num), int(color.green * num), int(color.blue *num)};
}

Colour operator+(Colour lhs, Colour rhs){
	return Colour{lhs.red + rhs.red, lhs.green + rhs.green, lhs.blue + rhs.blue};
}


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

void line(CanvasPoint &from, CanvasPoint &to, uint32_t color, DrawingWindow &window) {
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
		//test depth
		// auto cz = int(z) * 10;
		// auto c = packColor(Colour{cz, cz, cz});

		if ( (y-1) * HEIGHT + x >= 0 && (y-1) *HEIGHT + x < WIDTH * HEIGHT && zbuffer[(y-1) * HEIGHT + x] < z) {
			zbuffer[(y-1) * HEIGHT + x] = z;
			window.setPixelColour(x, y, color);
		} 
	}
}


void drawTriangle(CanvasTriangle triangle, Colour color, DrawingWindow &window) {
	for (int i = 0; i < WIDTH * HEIGHT; i++) {
		zbuffer[i] = 0.0;
	}

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

void drawFillTriangle(CanvasTriangle &triangle, Colour color, DrawingWindow &window) {
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
		line(CanvasPoint{xFrom, v0.y+i, zFrom}, CanvasPoint{xTo, v0.y+i, zTo}, colour, window);
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
		line(CanvasPoint{xFrom, vm.y+i, zFrom}, CanvasPoint{xTo, vm.y+i, zTo}, colour, window);
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

std::vector<float> barycentric(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C, const glm::vec3 &P) {
	const CanvasPoint a{A.x, A.y};
	const CanvasPoint b{B.x, B.y};
	const CanvasPoint c{C.x, C.y};
	const CanvasPoint p{P.x, P.y};

	return barycentric(a, b, c, p);
}

std::vector<float> TDBarycentric(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C, const glm::vec3 &P) {
	const float EPSILON = 1e-6f;
    std::vector<float> res(3, 0.0f); // 初始化结果为(0,0,0)

    // 1. 构建三角形局部坐标系的基向量（GLM直接支持向量减法）
    glm::vec3 v1 = A - C;  // 向量C->A
    glm::vec3 v2 = B - C;  // 向量C->B
    glm::vec3 p = P - C;   // 向量C->P

    // 2. 校验三角形是否退化（三点共线）：叉乘模长平方为0则共线
    glm::vec3 crossV1V2 = glm::cross(v1, v2);
    float crossLenSq = glm::dot(crossV1V2, crossV1V2); // 模长平方（避免开方）
    if (crossLenSq < EPSILON * EPSILON) {
        std::cerr << "Error: 三点共线，无法构成三角形！" << std::endl;
        return res; // 退化三角形返回(0,0,0)
    }

    // 3. 校验点P是否在三角形所在平面内
    float planeDot = glm::dot(p, crossV1V2);
    if (std::fabs(planeDot) > EPSILON) {
        std::cerr << "Warning: 点P不在三角形平面内，结果为近似值！" << std::endl;
        // 若需要严格校验，可直接返回(0,0,0)：return res;
    }

    // 4. 解线性方程组求重心坐标（GLM直接支持点乘）
    float a = glm::dot(v1, v1);  // v1·v1
    float b = glm::dot(v1, v2);  // v1·v2
    float c = glm::dot(v2, v2);  // v2·v2
    float d = glm::dot(p, v1);   // p·v1
    float e = glm::dot(p, v2);   // p·v2

    float denominator = a * c - b * b; // 方程组分母
    if (std::fabs(denominator) < EPSILON) {
        std::cerr << "Error: 计算分母为0，无法求解！" << std::endl;
        return res;
    }

    // 计算三个重心坐标分量
    float lambda1 = (d * c - b * e) / denominator; // 对应顶点A的权重
    float lambda2 = (a * e - b * d) / denominator; // 对应顶点B的权重
    float lambda3 = 1.0f - lambda1 - lambda2;      // 对应顶点C的权重

    // 赋值结果
    res[0] = lambda1;
    res[1] = lambda2;
    res[2] = lambda3;

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

decltype(auto) createVertex2FaceIndex(std::string modelFile) {
	std::ifstream file{modelFile, std::ios::in};
	if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << modelFile << '\n';
    }

	int verIndex = 0;
	int faceIndex = 0;
	std::vector<std::vector<int>> verInd2faceInd;
	std::vector<std::vector<int>> faceIndex2verInd;
	for (std::string str; std::getline(file, str);) {
		if (!str.empty()) {
			if(str.front() == 'v') {
				verInd2faceInd.emplace_back(std::vector<int>{});
			}else if (str.front() == 'f') {
				faceIndex2verInd.emplace_back(std::vector<int>{});
				std::vector<std::string> lineList = split(str, ' ');		
				verInd2faceInd[std::stoi(lineList[1]) - 1].emplace_back(faceIndex);
				verInd2faceInd[std::stoi(lineList[2]) - 1].emplace_back(faceIndex);
				verInd2faceInd[std::stoi(lineList[3]) - 1].emplace_back(faceIndex);
				faceIndex2verInd[faceIndex].emplace_back(std::stoi(lineList[1]) - 1);
				faceIndex2verInd[faceIndex].emplace_back(std::stoi(lineList[2]) - 1);
				faceIndex2verInd[faceIndex].emplace_back(std::stoi(lineList[3]) - 1);
				faceIndex++;
			}
		}
	}

	return std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>>{verInd2faceInd, faceIndex2verInd};
}

std::vector<glm::vec3> calculateVertexNormal(const std::vector<std::vector<int>> &verInd2faceInd, const std::vector<std::vector<int>> &faceInd2verInd, const std::vector<ModelTriangle> &model) {
	std::vector<glm::vec3> normalArr(verInd2faceInd.size(), glm::vec3{});

	for (int i = 0; i < verInd2faceInd.size(); i++) {
		auto faceList = verInd2faceInd[i];
		glm::vec3 aveNormal{0.0f, 0.0f, 0.0f};
		for (const auto face : faceList) {
			aveNormal += model[face].normal;
		}

		aveNormal = glm::normalize(aveNormal);
		normalArr[i] = aveNormal;
	}

	return normalArr;
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
				if (colorName.empty()) {
					colorName = "Red";
				}
				color = str2color[colorName];
				auto c = unpackColor(color);
				c.name = colorName;
				ModelTriangle mt{vertics[face.x-1], vertics[face.y-1], vertics[face.z-1], c};
				mt.normal = glm::normalize(glm::cross((mt.vertices[1] - mt.vertices[0]), (mt.vertices[2] - mt.vertices[0])));
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
	auto adjustVertex = vertex * cameraOrientation;
	auto u =  (focalLength * (adjustVertex.x) / abs(adjustVertex.z)) * scale + WIDTH / 2;
	auto v =  focalLength * (adjustVertex.y) / (adjustVertex.z) * scale + HEIGHT / 2;

	return glm::vec3{u, v, -1/adjustVertex.z};
}

decltype(auto) Rotate(rotateAxis axis, float radians) {
	switch (axis)
	{
	case rotateAxis::X_AXIS:{
		glm::vec3 col1(1.0f, 0.0f, 0.0f); 
		glm::vec3 col2(0.0f, std::cosf(radians), std::sinf(radians)); 
		glm::vec3 col3(0.0f, -std::sinf(radians), std::cosf(radians)); 
		return glm::mat3{col1, col2, col3};
		break;
	}
	
	case rotateAxis::Y_AXIS: {
		glm::vec3 col1(cosf(radians), 0.0f, -sinf(radians));
		glm::vec3 col2(0.0f, 1, 0.0); 
		glm::vec3 col3(sinf(radians), 0.0, std::cosf(radians)); 
		return glm::mat3{col1, col2, col3};
		break;
	}
	case rotateAxis::Z_AXIS: {
		glm::vec3 col1(cosf(radians), sinf(radians), 0.0f);
		glm::vec3 col2(-sinf(radians), cosf(radians), 0.0); 
		glm::vec3 col3(0.0, 0.0, 1.0); 
		return glm::mat3{col1, col2, col3};
		break;
	}
	default:
		return glm::mat3{};
		break;
	}

}

RayTriangleIntersection getClosestIntersection(glm::vec3 cameraPos, glm::vec3 rayDirect, ModelTriangle triangle) {
	glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
	glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
	glm::vec3 SPVector = cameraPos - triangle.vertices[0];
	glm::mat3 DEMatrix(-rayDirect, e0, e1);
	glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

	glm::vec3 intersectionPoint = triangle.vertices[0] + possibleSolution[1] * e0 + possibleSolution[2] * e1;
	glm::vec3 intersectionPointC = cameraPos + possibleSolution[0] * rayDirect;
	auto t = possibleSolution.x;
	auto u = possibleSolution.y;
	auto v = possibleSolution.z;

	if ((u >= 0.0 && u <= 1.0) && (v >= 0.0 && v <= 1.0) && ((u + v) <= 1.0) && (t > 0.0)) {
		assert(intersectionPoint.x - intersectionPointC.x < 0.000001 && intersectionPoint.y - intersectionPointC.y < 0.000001 && intersectionPoint.z - intersectionPointC.z < 0.000001);
		return RayTriangleIntersection{intersectionPoint, t, triangle, 0};
	}

	return RayTriangleIntersection{intersectionPoint, -1, triangle, 0};
}

RayTriangleIntersection getClosestIntersection(glm::vec3 startPos, glm::vec3 rayDirect, const std::vector<ModelTriangle> &triangles) {
	int index = 0;
	RayTriangleIntersection res;
	res.distanceFromCamera = -1.0;
	double maxDistance = std::numeric_limits<double>::infinity();
	for (const auto &triangle : triangles) {
		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = startPos - triangle.vertices[0];
		glm::mat3 DEMatrix(-rayDirect, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		glm::vec3 intersectionPoint = triangle.vertices[0] + possibleSolution[1] * e0 + possibleSolution[2] * e1;
		// glm::vec3 intersectionPointC = cameraPos + possibleSolution[0] * rayDirect;
		auto t = possibleSolution.x;
		auto u = possibleSolution.y;
		auto v = possibleSolution.z;
		if ((u >= 0.0 && u <= 1.0) && (v >= 0.0 && v <= 1.0) && ((u + v) <= 1.0) && (t > 0.0)) {
			if (maxDistance > t) {
				maxDistance = t;
				res.triangleIndex = index;
				res.distanceFromCamera = t;
				res.intersectedTriangle = triangle;
				res.intersectionPoint = intersectionPoint;
			}
		}

		index++;
	}

	return res;
}

bool lightVisible(glm::vec3 surfacePos, glm::vec3 lightPos, const std::vector<ModelTriangle> &triangles, size_t excludeIndex) {
	int index = 0;
	const float eps = 1e-6f;
	RayTriangleIntersection res;
	res.distanceFromCamera = std::numeric_limits<float>::infinity();
	auto distance = glm::distance(lightPos, surfacePos);
	for (const auto &triangle : triangles) {
		if (index == excludeIndex) {
			index++;
			continue;
		}

		glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
		glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
		glm::vec3 SPVector = surfacePos - triangle.vertices[0];
		auto rayDirect = glm::normalize(lightPos - surfacePos);
		glm::mat3 DEMatrix(-rayDirect, e0, e1);
		glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

		auto t = possibleSolution.x;
		auto u = possibleSolution.y;
		auto v = possibleSolution.z;
		if ((u >= 0.0 && u <= 1.0) && (v >= 0.0 && v <= 1.0) && ((u + v) <= 1.0) && (t > 0.0)) {
			if (res.distanceFromCamera > t) {
				res.triangleIndex = index;
				res.distanceFromCamera = t;
				res.intersectedTriangle = triangle;
			}
		}
		index++;
	}

	if (res.distanceFromCamera < distance && res.distanceFromCamera > eps) {
		return false;
	}

	return true;
}

Colour mixColor(Colour lhs, Colour rhs) {
	auto r = ((lhs.red/255.0f) * (rhs.red/255.0f)) * 255.0f;
	auto g = ((lhs.green/255.0f) * (rhs.green/255.0f)) * 255.0f;
	auto b = ((lhs.blue/255.0f) * (rhs.blue/255.0f)) * 255.0f;

	return Colour(int(r), int(g), int(b));
}
