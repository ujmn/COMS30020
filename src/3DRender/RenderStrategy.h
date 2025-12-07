#pragma once
#include <glm/glm.hpp>
#include <CanvasTriangle.h>
#include <DrawingWindow.h>
#include <RenderContext.h>

class Model;

class IRenderStrategy {
    public:
    explicit IRenderStrategy() = default;
    virtual ~IRenderStrategy() = default;
    virtual void render(const Model& model, std::shared_ptr<RenderContext> renderContext) = 0;
};


class IRenderStrategyFactory {
    public:
    IRenderStrategyFactory() = default;
    virtual ~IRenderStrategyFactory() = default;
    virtual std::unique_ptr<IRenderStrategy> CreateStrategy() const = 0;
};

template<typename T>
class RenderFactoryCRTP : public IRenderStrategyFactory{
    public:
    std::unique_ptr<IRenderStrategy> CreateStrategy() const override{
        static_assert(std::is_base_of_v<IRenderStrategy, T>, "StrategyT must inherit from IRenderStrategy");
        return std::make_unique<T>();
    }
};

class RasterStrategy : public IRenderStrategy {
    public:
    explicit RasterStrategy() = default;
    void render(const Model &model, std::shared_ptr<RenderContext> renderContext) override;
};
using RasterStrategyFactor = RenderFactoryCRTP<RasterStrategy>;

class RayTracingStrategy : public IRenderStrategy{
    public:
    explicit RayTracingStrategy() = default;
    virtual void render(const Model& model, std::shared_ptr<RenderContext> renderContext) override;
};
using RayTracingStrategyFactor = RenderFactoryCRTP<RayTracingStrategy>;

class PointCloundStrategy : public IRenderStrategy{
    public:
    explicit PointCloundStrategy() = default;
    virtual void render(const Model& model, std::shared_ptr<RenderContext> renderContext) override;
};
using PointCloundStrategyFactor = RenderFactoryCRTP<PointCloundStrategy>;

class WireFrameStrategy : public IRenderStrategy{
    public:
    explicit WireFrameStrategy() = default;
    virtual void render(const Model& model, std::shared_ptr<RenderContext> renderContext) override;
};
using WireFrameStrategyFactor = RenderFactoryCRTP<WireFrameStrategy>;

class SoftShadowStrategy : public IRenderStrategy{
    public:
    explicit SoftShadowStrategy() = default;
    virtual void render(const Model& model, std::shared_ptr<RenderContext> renderContext) override;
};
using SoftShadowStrategyFactor = RenderFactoryCRTP<SoftShadowStrategy>;

class MirrorReflectStrategy : public IRenderStrategy{
    public:
    explicit MirrorReflectStrategy() = default;
    virtual void render(const Model& model, std::shared_ptr<RenderContext> renderContext) override;
};
using MirrorReflectStrategyFactor = RenderFactoryCRTP<MirrorReflectStrategy>;

class GouraudShadingStrategy : public IRenderStrategy{
    public:
    explicit GouraudShadingStrategy() = default;
    virtual void render(const Model& model, std::shared_ptr<RenderContext> renderContext) override;
};
using GouraudShadingStrategyFactor = RenderFactoryCRTP<GouraudShadingStrategy>;

class PhongShadingStrategy : public IRenderStrategy{
    public:
    explicit PhongShadingStrategy() = default;
    virtual void render(const Model& model, std::shared_ptr<RenderContext> renderContext) override;
};
using PhongShadingStrategyFactor = RenderFactoryCRTP<PhongShadingStrategy>;