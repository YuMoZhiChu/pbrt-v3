
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


// cameras/orthographic.cpp*
#include "cameras/orthographic.h"
#include "paramset.h"
#include "sampler.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// OrthographicCamera Definitions
Float OrthographicCamera::GenerateRay(const CameraSample &sample,
                                      Ray *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute raster and camera sample positions
	// 在 pFilm 上采样
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
	// 转换到 Camera 空间
    Point3f pCamera = RasterToCamera(pFilm);
	// 因为是正交坐标系，所以射线方向默认是  0,0,1
    *ray = Ray(pCamera, Vector3f(0, 0, 1));
    // Modify ray for depth of field
	// 如果有景深，会对 Ray 的 O 和 d 做调整
    if (lensRadius > 0) {
        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // Compute point on plane of focus
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }
	// 时间用参数做插值
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->medium = medium;
	// 切换到世界坐标系
    *ray = CameraToWorld(*ray);
    return 1;
}

Float OrthographicCamera::GenerateRayDifferential(const CameraSample &sample,
                                                  RayDifferential *ray) const {
    ProfilePhase prof(Prof::GenerateCameraRay);
    // Compute main orthographic viewing ray

    // Compute raster and camera sample positions
    Point3f pFilm = Point3f(sample.pFilm.x, sample.pFilm.y, 0);
    Point3f pCamera = RasterToCamera(pFilm);
    *ray = RayDifferential(pCamera, Vector3f(0, 0, 1));

    // Modify ray for depth of field
    if (lensRadius > 0) {
        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);

        // Compute point on plane of focus
        Float ft = focalDistance / ray->d.z;
        Point3f pFocus = (*ray)(ft);

        // Update ray for effect of lens
        ray->o = Point3f(pLens.x, pLens.y, 0);
        ray->d = Normalize(pFocus - ray->o);
    }

    // Compute ray differentials for _OrthographicCamera_
    if (lensRadius > 0) {
        // Compute _OrthographicCamera_ ray differentials accounting for lens

        // Sample point on lens
        Point2f pLens = lensRadius * ConcentricSampleDisk(sample.pLens);
        Float ft = focalDistance / ray->d.z;

        Point3f pFocus = pCamera + dxCamera + (ft * Vector3f(0, 0, 1));
        ray->rxOrigin = Point3f(pLens.x, pLens.y, 0);
        ray->rxDirection = Normalize(pFocus - ray->rxOrigin);

        pFocus = pCamera + dyCamera + (ft * Vector3f(0, 0, 1));
        ray->ryOrigin = Point3f(pLens.x, pLens.y, 0);
        ray->ryDirection = Normalize(pFocus - ray->ryOrigin);
    } else {
		// 这种情况下，微分不管是多少倍的采样，ray->rxDirection 都不会变，只会变 ray->rxOrigin xy方向的微分距离
		// ray->rxDirection = ray->ryDirection = ray->d; 所以 ray->任何D 都是 (0,0,1)
        ray->rxOrigin = ray->o + dxCamera;
        ray->ryOrigin = ray->o + dyCamera;
        ray->rxDirection = ray->ryDirection = ray->d;
    }
    ray->time = Lerp(sample.time, shutterOpen, shutterClose);
    ray->hasDifferentials = true;
    ray->medium = medium;
    *ray = CameraToWorld(*ray);
    return 1;
}

OrthographicCamera *CreateOrthographicCamera(const ParamSet &params,
                                             const AnimatedTransform &cam2world,
                                             Film *film, const Medium *medium) {
    // Extract common camera parameters from _ParamSet_
    Float shutteropen = params.FindOneFloat("shutteropen", 0.f);
    Float shutterclose = params.FindOneFloat("shutterclose", 1.f);
    if (shutterclose < shutteropen) {
        Warning("Shutter close time [%f] < shutter open [%f].  Swapping them.",
                shutterclose, shutteropen);
        std::swap(shutterclose, shutteropen);
    }
    Float lensradius = params.FindOneFloat("lensradius", 0.f);
    Float focaldistance = params.FindOneFloat("focaldistance", 1e6f);
    Float frame = params.FindOneFloat(
        "frameaspectratio",
        Float(film->fullResolution.x) / Float(film->fullResolution.y));
    Bounds2f screen;
    if (frame > 1.f) {
        screen.pMin.x = -frame;
        screen.pMax.x = frame;
        screen.pMin.y = -1.f;
        screen.pMax.y = 1.f;
    } else {
        screen.pMin.x = -1.f;
        screen.pMax.x = 1.f;
        screen.pMin.y = -1.f / frame;
        screen.pMax.y = 1.f / frame;
    }
    int swi;
    const Float *sw = params.FindFloat("screenwindow", &swi);
    if (sw) {
        if (swi == 4) {
            screen.pMin.x = sw[0];
            screen.pMax.x = sw[1];
            screen.pMin.y = sw[2];
            screen.pMax.y = sw[3];
        } else
            Error("\"screenwindow\" should have four values");
    }
    return new OrthographicCamera(cam2world, screen, shutteropen, shutterclose,
                                  lensradius, focaldistance, film, medium);
}

}  // namespace pbrt
