#pragma once
#include <Model.h>
#include <DrawingWindow.h>
#include <memory>
#include <SDL_events.h>

class IRenderStrategy;
class RenderContext;
class IRenderStrategyFactory;
struct CanvasTriangle;
struct CanvasPoint;

class IScene{
    public:
    explicit IScene() = default;
    virtual ~IScene() = default;
    virtual void update() = 0;
};

class Scene : public IScene{
    private:
    explicit Scene(int width, int height, glm::vec3 cameraPos, float focalLength);
    public:
    virtual ~Scene() = default;
    virtual void update() override;
    void processOrbit();
    void addModel(Model model);
    static Scene& getScene(int width, int height, glm::vec3 cameraPos, float focalLength);
    void addRenderStrategy(std::unique_ptr<IRenderStrategy> renderStartegy);
    void start();
    void handleEvent(SDL_Event event, DrawingWindow &window);
	void setOrbit();
    void setLightPos(glm::vec3 lightPos);

    private:
    std::shared_ptr<RenderContext>                  renderContext_;
	SDL_Event                                       event_;
    std::unique_ptr<IRenderStrategyFactory>         renderStrategyFactory_;         
    std::vector<std::unique_ptr<IRenderStrategy>>   renderStrategyList_;
    std::vector<Model>                              models_;
    size_t                                          renderStrategyIndex_{0};
};

class ModelRender{
    public:
    explicit ModelRender() = default;
    void drawFillTriangle(CanvasTriangle &triangle, Colour color, std::shared_ptr<RenderContext> context);
    void drawTriangle(CanvasTriangle &triangle, Colour color, std::shared_ptr<RenderContext> context);
    void line(CanvasPoint &from, CanvasPoint &to, uint32_t color,  std::shared_ptr<RenderContext> context);
}; 
