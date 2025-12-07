---
typora-root-url: resource
---

```markdown
<div align="center">[Chinese](#chinese) | [English](#english)</div>
```

<a id="chinese"></a>>

# 3D Renderer

3D Renderer是基于 **C++** 实现的简易 3D 渲染器，支持多种经典图形学渲染技术，可在 Windows 11 环境下使用 MSVC 2022 编译运行。

## ✨ 核心功能

- [点云渲染（Point Cloud）](#point-cloud)
- [线框模式（Wireframe）](#wire-frame)
- [光栅化（Rasterization）](#rasterization)
- [光线追踪（Ray Tracing）](#ray-tracing)
- [软阴影（Soft Shadows）](#soft-shadows)
- [光照模型（Lighting Model）](#lighting-model)
- [镜面反射（Mirror Reflection）](#mirror-reflection)
- [球面光源模型（Sphere Lighting）](#sphere-lighting)
- [Gouraud着色（Gouraud Shading）](#gouraud-shading)
- [Phong着色（Phong Shading）](#phong-shading)

## 🛠️ 开发环境

- 操作系统：Windows 11
- 编译器：MSVC 2022
- 语言：C++

> 项目无外部图形库依赖（如 OpenGL），所有渲染逻辑均为自主实现。

<a id="english"></a>

# 3D Renderer

The 3D renderer is a simple implementation in **C++**, supporting multiple classic computer graphics rendering techniques, and can be compiled and run on Windows 11 using MSVC 2022.

## ✨ Features

- [Point Cloud](#point-cloud)
- [Wireframe](#wire-frame)
- [Rasterization](#rasterization)
- [Ray Tracing](#ray-tracing)
- [Soft Shadows](#soft-shadows)
- [Lighting Model](#lighting-model)
- [Mirror Reflection](#mirror-reflection)
- [Sphere Lighting](#sphere-lighting)
- [Gouraud Shading](#gouraud-shading)
-  [Phong Shading](#phong-shading)

## 🛠️ Environment

- OS: Windows 11
- Compiler: MSVC 2022
- Language: C++

> The project has no dependencies on external graphics libraries (such as  OpenGL); all rendering logic is implemented from scratch.



<a id="point-cloud"></a>

## Point Cloud

![](pointCloud.png)

<a id="wire-frame"></a>

## Wire Frame

![](wireFrame.png)

<a id="rasterization"></a>

## Rasterization

![](raster.png)

<a id="ray-tracing"></a>

## Ray Tracing

![](rayTracing.png)

<a id="soft-shadows"></a>

## Soft Shadows

![](softShadow.png)

<a id="lighting-model"></a>

## Lighting Model

![](lightModel.png)

<a id="mirror-reflection"></a>

## Mirror Reflection

![](mirrorReflect.png)

<a id="sphere-lighting-model"></a>

## Sphere Lighting Model

![](sphere_lightModel.png)

<a id="gouraud-shading"></a>

## Gouraud Shading

![](sphere_gouraudShading.png)

<a id="phong-shading"></a>

## Phong Shading

![](sphere_phongShading.png)