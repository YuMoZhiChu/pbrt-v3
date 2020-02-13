
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

#ifndef PBRT_CORE_SHAPE_H
#define PBRT_CORE_SHAPE_H

// core/shape.h*
#include "pbrt.h"
#include "geometry.h"
#include "interaction.h"
#include "memory.h"
#include "transform.h"

namespace pbrt {

// Shape Declarations
// Sideness 说明: pbrt 中没有去做剪切物体的背面, 虽然剪切背面能优化渲染, 但 Ray-Intersect 逻辑会在设置正背面之前常常调用, 因此没必要
class Shape {
  public:
    // Shape Interface
    Shape(const Transform *ObjectToWorld, const Transform *WorldToObject,
          bool reverseOrientation);
    virtual ~Shape();
	// 返回该 Shape 的范围, 在本地坐标系
    virtual Bounds3f ObjectBound() const = 0;
	// 返回该 Shape 的范围, 在世界坐标系
    virtual Bounds3f WorldBound() const;
	// 纯虚函数, 但是子类应该提供与射线的交点 Interaction 以及 tHit
	// 该函数的注意点:
	// 1. Ray中包含 Ray::tMax 来表示末端, 在 t > tMax 的情况下, 这一次的 Interaction 需要被忽略
	// 2. tHit 用于储存 Intereaction 对应的 Ray 的 t值, 在多个交点的情况下, 存储最近的一个
	// 3. isect 用于储存 SurfaceInteraction, 表面交点
	// 4. 传入的 Ray 是在 世界坐标系中表示, Shape中的计算会把他们转换到 本地坐标系, 返回的 Interaction 信息则是世界坐标系的
	// 5. testAlphaTexture 表示是否采用: 使用一些纹理来切掉 Shape 的部分表面
    virtual bool Intersect(const Ray &ray, Float *tHit,
                           SurfaceInteraction *isect,
                           bool testAlphaTexture = true) const = 0;
	// 相比于 Intersect, 这个不需要计算 tHit, isect 是一个简化版本
    virtual bool IntersectP(const Ray &ray,
                            bool testAlphaTexture = true) const {
        return Intersect(ray, nullptr, nullptr, testAlphaTexture);
    }
    virtual Float Area() const = 0;
    // Sample a point on the surface of the shape and return the PDF with
    // respect to area on the surface.
    virtual Interaction Sample(const Point2f &u, Float *pdf) const = 0;
    virtual Float Pdf(const Interaction &) const { return 1 / Area(); }

    // Sample a point on the shape given a reference point |ref| and
    // return the PDF with respect to solid angle from |ref|.
    virtual Interaction Sample(const Interaction &ref, const Point2f &u,
                               Float *pdf) const;
    virtual Float Pdf(const Interaction &ref, const Vector3f &wi) const;

    // Returns the solid angle subtended by the shape w.r.t. the reference
    // point p, given in world space. Some shapes compute this value in
    // closed-form, while the default implementation uses Monte Carlo
    // integration; the nSamples parameter determines how many samples are
    // used in this case.
    virtual Float SolidAngle(const Point3f &p, int nSamples = 512) const;

    // Shape Public Data
	// Shape 定位在对象空间中, 所有Shape都定义在本地坐标系中, 提供2个矩阵支持相应的变化
    const Transform *ObjectToWorld, *WorldToObject;
	// 该参数表示Shape 的 Normal 是否翻转, 对于自发光的物体, 这个参数会影响光照的发射方向
    const bool reverseOrientation;
	// 是否因为 transform 而改变了手系, 这个参数会在构造函数中做计算并保存
    const bool transformSwapsHandedness;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SHAPE_H
