#include <RenderUtils.h>

uint32_t RenderUtils::Utils::packColor(Colour color)
{
	return (255 << 24) + (color.red << 16) + (color.green << 8) + int(color.blue);
}

Colour RenderUtils::Utils::unpackColor(uint32_t color)
{
	Colour c;
	c.red = (color & 0x00FF0000) >> 16;
	c.green = (color & 0x0000FF00) >> 8;
	c.blue = (color & 0x000000FF);

	return c;
}

glm::mat3 RenderUtils::Utils::Rotate(RenderUtils::rotateAxis axis, float radians)
{
	switch (axis)
	{
	case RenderUtils::rotateAxis::X_AXIS:{
		glm::vec3 col1(1.0f, 0.0f, 0.0f); 
		glm::vec3 col2(0.0f, std::cosf(radians), std::sinf(radians)); 
		glm::vec3 col3(0.0f, -std::sinf(radians), std::cosf(radians)); 
		return glm::mat3{col1, col2, col3};
		break;
	}
	
	case RenderUtils::rotateAxis::Y_AXIS: {
		glm::vec3 col1(cosf(radians), 0.0f, -sinf(radians));
		glm::vec3 col2(0.0f, 1, 0.0); 
		glm::vec3 col3(sinf(radians), 0.0, std::cosf(radians)); 
		return glm::mat3{col1, col2, col3};
		break;
	}
	case RenderUtils::rotateAxis::Z_AXIS: {
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

glm::mat3 RenderUtils::Utils::lookAt(glm::vec3 cameraPos)
{
	glm::vec3 forward = glm::normalize(cameraPos - glm::vec3(0.0, 0.0, 0.0));
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3{0.0, 1.0, 0.0}, forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));

    glm::mat3 cameraOrientation;
	cameraOrientation[0] = right;
	cameraOrientation[1] = up;
	cameraOrientation[2] = forward;
    return cameraOrientation;
}


glm::vec3 RenderUtils::Utils::projectVertiexOntoCanvasPoint(glm::vec3 cameraPos, glm::mat3 cameraOrientation, float focalLength, glm::vec3 vertex, int width, int height, float scale = 1.0) {
	vertex = vertex - cameraPos;
	auto adjustVertex = vertex * cameraOrientation;
	auto u =  (focalLength * (adjustVertex.x) / abs(adjustVertex.z)) * scale + width / 2;
	auto v =  focalLength * (adjustVertex.y) / (adjustVertex.z) * scale + height / 2;

	return glm::vec3{u, v, -1/adjustVertex.z};
}

RayTriangleIntersection RenderUtils::Utils::getClosestIntersection(glm::vec3 startPos, glm::vec3 rayDirect, const std::vector<ModelTriangle> &triangles) {
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

bool RenderUtils::Utils::lightVisible(glm::vec3 surfacePos, glm::vec3 lightPos, const std::vector<ModelTriangle> &triangles, size_t excludeIndex) {
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

// 同心映射：将 [0,1)^2 -> 单位圆盘 (返回 {dx, dz})
glm::vec3 RenderUtils::Utils::squareToDiskConcentricXZ(double u, double v) {
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

std::vector<glm::vec3> RenderUtils::Utils::sampleDiskOnXZPlane(const glm::vec3& center, double radius, int numSamples) {
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

Colour RenderUtils::Utils::mixColor(Colour lhs, Colour rhs) {
	auto r = ((lhs.red/255.0f) * (rhs.red/255.0f)) * 255.0f;
	auto g = ((lhs.green/255.0f) * (rhs.green/255.0f)) * 255.0f;
	auto b = ((lhs.blue/255.0f) * (rhs.blue/255.0f)) * 255.0f;

	return Colour(int(r), int(g), int(b));
}

std::vector<float> RenderUtils::Utils::TDBarycentric(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C, const glm::vec3 &P) {
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
        std::cerr << "Error: 三点共线，无法构成三角形!" << std::endl;
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

Colour operator*(Colour color, float num) {
	return Colour{int(color.red * num), int(color.green * num), int(color.blue *num)};
}

Colour operator+(Colour lhs, Colour rhs){
	return Colour{lhs.red + rhs.red, lhs.green + rhs.green, lhs.blue + rhs.blue};
}