#include <Model.h>
#include <Scene.h>
#include <iostream>


int main(int argc, char *argv[])
{
    auto &scene = Scene::getScene(640, 480, glm::vec3(0.0f, 0.0f, 15.0f), 10.0);
    Model model{"models/cornell-box.obj", "models/cornell-box.mtl"}; //康奈尔盒子
	scene.setLightPos(glm::vec3{0.001f, 2.3389f, 1.507f}); //康奈尔盒子的灯坐标
    // Model model{"models/sphere.obj", "models/cornell-box.mtl"}; 
	// scene.setLightPos(glm::vec3{0.001f, 1.2f, 6.207f}); //球的坐标
    std::vector<int> mirrorFaces{26, 31};
    model.setMirrorFaces(mirrorFaces);
    scene.addModel(model);
    scene.start();
    return 0;
}