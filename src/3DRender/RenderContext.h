#pragma once
#include <DrawingWindow.h>
#include <iostream>
#include <memory>
#include <glm/glm.hpp>

class ModelRender;
class Scene;

class RenderContext{
    friend class Scene;
    public:
    explicit RenderContext() = default;
    explicit RenderContext(glm::vec3 cameraPos, float focalLength, int width, int height);
    ~RenderContext() = default;

    std::shared_ptr<DrawingWindow> window() const;
    bool orbit() const;
    glm::vec3 cameraPos() const;
    void setCameraPos(glm::vec3 cameraPos);
    void setCameraOrientation(glm::mat3 cameraOrientation);
    glm::mat3 cameraOrientation() const;
    float focalLength() const;
    std::vector<float>& zbuffer();
    void clearZbuffer();
    std::shared_ptr<ModelRender> modelRender() const;
    float scale() const;
    glm::vec3 lightPos() const;
    size_t sampleCount() const;
    float sampleRadius() const;


    private:
    std::shared_ptr<DrawingWindow>          window_;
    bool                                    orbit_{false};
    glm::vec3                               cameraPos_;
    glm::mat3                               cameraOrientation_;
    float                                   focalLength_;
    std::vector<float>                      zbuffer_;
    std::shared_ptr<ModelRender>            modelRender_;
    float                                   scale_{80.0f};
    glm::vec3                               lightPos_;
    size_t                                  sampleCount_{32};
    float                                   sampleRadius_{0.25};
};