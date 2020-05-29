
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_CAMERA_H
#define PBRT_CORE_CAMERA_H

// core/camera.h*
#include "pbrt.h"
#include "geometry.h"
#include "transform.h"
#include "film.h"

namespace pbrt {

// Camera Declarations
class Camera {
  public:
    // Camera Interface
	// 相机的初始化，最重要的就是 CameraToWorld 这个参数代表坐标系的转换，而且它是一个 AnimatedTransform，可以拆解成 TRS，进行动态插值
    Camera(const AnimatedTransform &CameraToWorld, Float shutterOpen,
           Float shutterClose, Film *film, const Medium *medium);
    virtual ~Camera();
	// 对给定的采样点 CameraSample
	// returns a floating-point weight associated with the ray
	// 返回一个浮点值，这个值会影响最终的权重，简单的摄像机返回的是1，真实的摄像机，会根据光穿透镜片携带多少能量，来计算这个数
	// 在简单照相机中, Ray 的权重是一样的
	// 05-28：给入 Sample， 计算 Ray，这里的 Ray.d 会做 normalized
    virtual Float GenerateRay(const CameraSample &sample, Ray *ray) const = 0;
	// 根据采样点 CameraSample 生成 Ray, 还包含了 在 Pixel 在 Image 平面上 关于 x,y方向 的微分
	// 微分的数据 -> pixel spacing -> texture antialiasing 抗锯齿算法
    virtual Float GenerateRayDifferential(const CameraSample &sample,
                                          RayDifferential *rd) const;
    virtual Spectrum We(const Ray &ray, Point2f *pRaster2 = nullptr) const;
    virtual void Pdf_We(const Ray &ray, Float *pdfPos, Float *pdfDir) const;
    virtual Spectrum Sample_Wi(const Interaction &ref, const Point2f &u,
                               Vector3f *wi, Float *pdf, Point2f *pRaster,
                               VisibilityTester *vis) const;

    // Camera Public Data
    AnimatedTransform CameraToWorld;
	// 快门的开关时间点，用于实现动态模糊的效果
    const Float shutterOpen, shutterClose;
	// handles image storage
	// 指向一张 film，用于显示最终的图像
    Film *film;
	// 介质，照相机所处的散射环境的介质（???? 可能是水下拍摄的做法
    const Medium *medium;
};

// 这里面的信息，足够用来生成一套射线了
struct CameraSample {
	// the position on the film for which the camera should generate the corresponding ray
	// film 上的位置
	// 在 film 上，能找到对应的点，和对应的 radiance
    Point2f pFilm;
	// for camera models that simulate non-pinhole apertures
	// ???? 用于非针孔模型的摄像机
	// ???? 啥意思 The point on the lens the ray passes through is in pLens (for cameras that include the notion of lenses),
    Point2f pLens;
	// used when rendering scenes with moving objects
	// 渲染移动物体时用
	// interpolate within the shutterOpen–shutterClose time range
    Float time;
};

inline std::ostream &operator<<(std::ostream &os, const CameraSample &cs) {
    os << "[ pFilm: " << cs.pFilm << " , pLens: " << cs.pLens <<
        StringPrintf(", time %f ]", cs.time);
    return os;
}

class ProjectiveCamera : public Camera {
  public:
    // ProjectiveCamera Public Methods
	// screenWindow: 屏幕大小范围
	// focald: 聚焦距离，用于模拟真实相机的失焦
    ProjectiveCamera(const AnimatedTransform &CameraToWorld,
                     const Transform &CameraToScreen, // 投影矩阵
                     const Bounds2f &screenWindow, Float shutterOpen,
                     Float shutterClose, Float lensr, Float focald, Film *film,
                     const Medium *medium)
        : Camera(CameraToWorld, shutterOpen, shutterClose, film, medium),
          CameraToScreen(CameraToScreen) {
        // Initialize depth of field parameters
        lensRadius = lensr;
        focalDistance = focald;

        // Compute projective camera transformations

        // Compute projective camera screen transformations
        ScreenToRaster =
            Scale(film->fullResolution.x, film->fullResolution.y, 1) *
            Scale(1 / (screenWindow.pMax.x - screenWindow.pMin.x),
                  1 / (screenWindow.pMin.y - screenWindow.pMax.y), 1) *
            Translate(Vector3f(-screenWindow.pMin.x, -screenWindow.pMax.y, 0));
		// 这里是右乘的乘法，所以从下往上读
        RasterToScreen = Inverse(ScreenToRaster);
        RasterToCamera = Inverse(CameraToScreen) * RasterToScreen; // TODO 这里是线性代数知识，左手系+行优先矩阵，就是这种乘法，右乘的乘法
    }

  protected:
    // ProjectiveCamera Protected Data
    Transform CameraToScreen, RasterToCamera;
    Transform ScreenToRaster, RasterToScreen;
    Float lensRadius, focalDistance;
};

}  // namespace pbrt

#endif  // PBRT_CORE_CAMERA_H
