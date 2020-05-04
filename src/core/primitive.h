
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

#ifndef PBRT_CORE_PRIMITIVE_H
#define PBRT_CORE_PRIMITIVE_H

// core/primitive.h*
#include "pbrt.h"
#include "shape.h"
#include "material.h"
#include "medium.h"
#include "transform.h"

namespace pbrt {

// Primitive Declarations
// Primitive 是 几何处理 和 渲染子系统 的链接桥梁
// 提供了许多和 Shape 相似的接口
class Primitive {
  public:
    // Primitive Interface
    virtual ~Primitive();
	// 返回一个包裹 Primitive 的包围框(世界坐标), 可以放置在 acceleration data structures(加速结构中, 比如kd树
    virtual Bounds3f WorldBound() const = 0;
	// 不同于 Shape 类的 Intersect 的 float* 值, Primitive 的 Intersect 会更新 Ray 的 tMax
	// 更新 Ray 的 tMax, 能够缩短射线, 跳过一些过远的物体
    virtual bool Intersect(const Ray &r, SurfaceInteraction *) const = 0;
    virtual bool IntersectP(const Ray &r) const = 0;
	// 如果这个基元是发光体, 这里会返回一个 AreaLight 来描述它的  emission distribution 自发光分布
    virtual const AreaLight *GetAreaLight() const = 0;
	// 返回该基元的材质, 如果是 nullptr, 会跳过 ray intersect
    virtual const Material *GetMaterial() const = 0;
	// 次表面散射相关
	// MemoryArena 为 BSDF/BSSRDF 分配的内存
	// TransportMode 枚举类型, 此光照射线的起点是 相机/光源
	// allowMultipleLobes BRDF 的表现方式
    virtual void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                            MemoryArena &arena,
                                            TransportMode mode,
                                            bool allowMultipleLobes) const = 0;
};

// GeometricPrimitive Declarations
// 用于表示场景中的 单个Shape (a single shape)
class GeometricPrimitive : public Primitive {
  public:
    // GeometricPrimitive Public Methods
    virtual Bounds3f WorldBound() const;
    virtual bool Intersect(const Ray &r, SurfaceInteraction *isect) const;
    virtual bool IntersectP(const Ray &r) const;
    GeometricPrimitive(const std::shared_ptr<Shape> &shape,
                       const std::shared_ptr<Material> &material,
                       const std::shared_ptr<AreaLight> &areaLight,
                       const MediumInterface &mediumInterface);
    const AreaLight *GetAreaLight() const;
    const Material *GetMaterial() const;
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const;

  private:
    // GeometricPrimitive Private Data
	// 指向 Shape 的指针
    std::shared_ptr<Shape> shape;
	// 指向 材质 的指针
    std::shared_ptr<Material> material;
	// 如果是自发光
    std::shared_ptr<AreaLight> areaLight;
	// 内部和外部的 介质
    MediumInterface mediumInterface;
};

// TransformedPrimitive Declarations
class TransformedPrimitive : public Primitive {
  public:
    // TransformedPrimitive Public Methods
    TransformedPrimitive(std::shared_ptr<Primitive> &primitive,
                         const AnimatedTransform &PrimitiveToWorld);
    bool Intersect(const Ray &r, SurfaceInteraction *in) const;
    bool IntersectP(const Ray &r) const;
    const AreaLight *GetAreaLight() const { return nullptr; }
    const Material *GetMaterial() const { return nullptr; }
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const {
        LOG(FATAL) <<
            "TransformedPrimitive::ComputeScatteringFunctions() shouldn't be "
            "called";
    }
    Bounds3f WorldBound() const {
		// 这里也需要做矩阵处理
        return PrimitiveToWorld.MotionBounds(primitive->WorldBound());
    }

  private:
    // TransformedPrimitive Private Data
	// 单个基元
    std::shared_ptr<Primitive> primitive;
	// 动画矩阵, 实现2个特性
	// 1. object instancing 物体实例化 (不去复制物体的 Vertex 等信息
	// 2. primitives with animated transformations 基元的动画变换 (不需要所有的Shape都支持 animated 挂载 AnimatedTransform 即可
	// 如果shape有 PrimitiveToWorld, 做了 PrimitiveToWorld 的矩阵变换 才是 shape 在世界坐标中的位置
    const AnimatedTransform PrimitiveToWorld;
};

// Aggregate Declarations
class Aggregate : public Primitive {
  public:
    // Aggregate Public Methods
	// 这些函数不应该被调用
    const AreaLight *GetAreaLight() const;
    const Material *GetMaterial() const;
	// isect 是返回具体的相交信息
    void ComputeScatteringFunctions(SurfaceInteraction *isect,
                                    MemoryArena &arena, TransportMode mode,
                                    bool allowMultipleLobes) const;
};

}  // namespace pbrt

#endif  // PBRT_CORE_PRIMITIVE_H
