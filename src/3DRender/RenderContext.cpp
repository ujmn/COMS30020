#include "RenderContext.h"
#include <Scene.h>

RenderContext::RenderContext(glm::vec3 cameraPos, float focalLength, int width, int height)
    : window_(std::make_shared<DrawingWindow>(DrawingWindow{width, height, false}))
    , cameraPos_(cameraPos)
    , focalLength_(focalLength)   
    , zbuffer_(std::vector<float>(width * height, 0.0))
    , modelRender_(std::make_shared<ModelRender>(ModelRender{}))
{
	glm::vec3 forward = glm::normalize(cameraPos - glm::vec3(0.0, 0.0, 0.0));
	glm::vec3 right = glm::normalize(glm::cross(glm::vec3{0.0, 1.0, 0.0}, forward));
	glm::vec3 up = glm::normalize(glm::cross(forward, right));

	cameraOrientation_[0] = right;
	cameraOrientation_[1] = up;
	cameraOrientation_[2] = forward;
}

std::shared_ptr<DrawingWindow> RenderContext::window() const
{
    return window_;
}

bool RenderContext::orbit() const
{
    return orbit_;
}

glm::vec3 RenderContext::cameraPos() const
{
    return cameraPos_;
}

void RenderContext::setCameraPos(glm::vec3 cameraPos)
{
    cameraPos_ = cameraPos;
}

void RenderContext::setCameraOrientation(glm::mat3 cameraOrientation)
{
    cameraOrientation_ = cameraOrientation;
}

glm::mat3 RenderContext::cameraOrientation() const
{
    return cameraOrientation_;
}

float RenderContext::focalLength() const
{
    return focalLength_;
}

std::vector<float>& RenderContext::zbuffer()
{
    return zbuffer_;
}

void RenderContext::clearZbuffer()
{
    for (auto &buffer : zbuffer_) {
        buffer = std::numeric_limits<float>::min();
    }
}

std::shared_ptr<ModelRender> RenderContext::modelRender() const
{
    return modelRender_;
}

float RenderContext::scale() const
{
    return scale_;
}

glm::vec3 RenderContext::lightPos() const
{
    return lightPos_;
}

size_t RenderContext::sampleCount() const
{
    return sampleCount_;
}

float RenderContext::sampleRadius() const
{
    return sampleRadius_;
}
