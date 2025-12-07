#pragma once

#include <glm/glm.hpp>
#include <iostream>
#include <string>
#include <vector>
#include <ModelTriangle.h>

class RenderContext;
class IRenderStrategy;

class RenderableObject{
    public:
        explicit RenderableObject() = default;
        virtual ~RenderableObject() = default;
        virtual void loadFile(std::string objFilePath, std::string mtlFilePath) = 0;
        virtual void draw(std::shared_ptr<RenderContext> context) = 0;
};

class Model : public RenderableObject{
    public:
        explicit Model(std::string objFilePath, std::string mtlFilePath);
        virtual ~Model();
        Model(const Model &rhs);
        virtual void loadFile(std::string objFilePath, std::string mtlFilePath) override;
        virtual void draw(std::shared_ptr<RenderContext> context) override;
        const std::vector<ModelTriangle>& triangles() const;
        void setMirrorFaces(std::vector<int> mirrorFaces);
        std::vector<int> mirrorFaceIndex() const;
        std::tuple<std::vector<std::vector<int>>, std::vector<std::vector<int>>> createVertex2FaceIndex();
        std::vector<glm::vec3> calculateVertexNormal();
        const std::vector<std::vector<int>>& verInd2faceInd() const;
        const std::vector<std::vector<int>>& faceIndex2verInd() const;
        const std::vector<glm::vec3>& vertexNormal() const;


    private:
    std::string                                     objFilePath_{};
    std::string                                     mtlFilePath_{};
    std::vector<ModelTriangle>                      modelTriangles_{};
    std::vector<int>                                mirrorFaces_{};
    std::shared_ptr<IRenderStrategy>                renderStrategy_{};
    std::vector<std::vector<int>>                   verInd2faceInd_{};
    std::vector<std::vector<int>>                   faceIndex2verInd_{};
    std::vector<glm::vec3>                          vertexNormal_{};
};