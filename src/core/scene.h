
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

#ifndef PBRT_CORE_SCENE_H
#define PBRT_CORE_SCENE_H

// core/scene.h*
#include "pbrt.h"
#include "geometry.h"
#include "primitive.h"
#include "light.h"

namespace pbrt {

// Scene Declarations
class Scene {
  public:
    // Scene Public Methods
    Scene(std::shared_ptr<Primitive> aggregate,
          const std::vector<std::shared_ptr<Light>> &lights)
        : lights(lights), aggregate(aggregate) {
        // Scene Constructor Implementation
        worldBound = aggregate->WorldBound();
        for (const auto &light : lights) {
			// 对一些 光照Light 进行 预处理Preprocess , 在 场景Scene 定义之后, 开始渲染之前
            light->Preprocess(*this);
            if (light->flags & (int)LightFlags::Infinite)
                infiniteLights.push_back(light);
        }
    }
    const Bounds3f &WorldBound() const { return worldBound; }
	// 获得 Ray 和 Primitive 的 Intersect 函数, 如果 True, 参数 isect 会被填入数据
    bool Intersect(const Ray &ray, SurfaceInteraction *isect) const;
	// 获得 Ray 和 Primitive 的 Intersect 函数, 不去计算最近的 Intersection, 效率更高
    bool IntersectP(const Ray &ray) const;
    bool IntersectTr(Ray ray, Sampler &sampler, SurfaceInteraction *isect,
                     Spectrum *transmittance) const;

    // Scene Public Data
	// 这里的 光照Light 是全局的, 光会对场景中的所有物体造成影响(没有剔除, 分块等优化操作)
    std::vector<std::shared_ptr<Light>> lights;
    // Store infinite light sources separately for cases where we only want
    // to loop over them.
    std::vector<std::shared_ptr<Light>> infiniteLights;

  private:
    // Scene Private Data
	// 所有在场景中的几何体都用 Primitive 表示, 结合了 Shape 和 Material
	// 类 Aggregate 是一个几何体, 包含其他 Primitive 的引用
    std::shared_ptr<Primitive> aggregate;
    Bounds3f worldBound;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SCENE_H
