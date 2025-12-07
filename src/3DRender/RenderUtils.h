#pragma once

#include <glm/glm.hpp>
#include <Colour.h>
#include <DrawingWindow.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>

class RenderContext;

namespace RenderUtils
{
    enum class rotateAxis{
    	X_AXIS,
    	Y_AXIS,
    	Z_AXIS
    };

    class Utils{
        public:
        Utils() = delete;
        ~Utils() = delete;
        static uint32_t packColor(Colour color);
        static Colour unpackColor(uint32_t color);
        static glm::mat3 Rotate(RenderUtils::rotateAxis axis, float radians);
        static glm::mat3 lookAt(glm::vec3 cameraPos);
        static glm::vec3 projectVertiexOntoCanvasPoint(glm::vec3 cameraPos, glm::mat3 cameraOrientation, float focalLength, glm::vec3 vertex, int width, int height, float scale);
        static RayTriangleIntersection RenderUtils::Utils::getClosestIntersection(glm::vec3 startPos, glm::vec3 rayDirect, const std::vector<ModelTriangle> &triangles);
        static bool lightVisible(glm::vec3 surfacePos, glm::vec3 lightPos, const std::vector<ModelTriangle> &triangles, size_t excludeIndex);
        static glm::vec3 squareToDiskConcentricXZ(double u, double v);
        static std::vector<glm::vec3> sampleDiskOnXZPlane(const glm::vec3& center, double radius, int numSamples);
        static Colour mixColor(Colour lhs, Colour rhs);
        static std::vector<float> TDBarycentric(const glm::vec3 &A, const glm::vec3 &B, const glm::vec3 &C, const glm::vec3 &P);
    };

} // namespace RenderUtils


Colour operator*(Colour color, float num);
Colour operator+(Colour lhs, Colour rhs);

