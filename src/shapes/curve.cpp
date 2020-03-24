
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


// shapes/curve.cpp*
#include "shapes/curve.h"
#include "paramset.h"
#include "stats.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Curves", curveBytes);
STAT_PERCENT("Intersections/Ray-curve intersection tests", nHits, nTests);
STAT_INT_DISTRIBUTION("Intersections/Curve refinement level", refinementLevel);
STAT_COUNTER("Scene/Curves", nCurves);
STAT_COUNTER("Scene/Split curves", nSplitCurves);

// Curve Utility Functions
static Point3f BlossomBezier(const Point3f p[4], Float u0, Float u1, Float u2) {
	// 贝塞尔一阶插值
    Point3f a[3] = {Lerp(u0, p[0], p[1]), Lerp(u0, p[1], p[2]),
                    Lerp(u0, p[2], p[3])};
	// 二阶
    Point3f b[2] = {Lerp(u1, a[0], a[1]), Lerp(u1, a[1], a[2])};
	// 三阶
    return Lerp(u2, b[0], b[1]);
}

// 贝塞尔曲线的2段拆分
// 将4个控制点分为 7个控制点, 分别控制2段 子曲线
inline void SubdivideBezier(const Point3f cp[4], Point3f cpSplit[7]) {
    cpSplit[0] = cp[0];
    cpSplit[1] = (cp[0] + cp[1]) / 2;
    cpSplit[2] = (cp[0] + 2 * cp[1] + cp[2]) / 4;
    cpSplit[3] = (cp[0] + 3 * cp[1] + 3 * cp[2] + cp[3]) / 8;
    cpSplit[4] = (cp[1] + 2 * cp[2] + cp[3]) / 4;
    cpSplit[5] = (cp[2] + cp[3]) / 2;
    cpSplit[6] = cp[3];
}

// 验证贝塞尔曲线, 带入 CP u 得到贝塞尔的当前点, 顺便算个微分
static Point3f EvalBezier(const Point3f cp[4], Float u,
                          Vector3f *deriv = nullptr) {
    Point3f cp1[3] = {Lerp(u, cp[0], cp[1]), Lerp(u, cp[1], cp[2]),
                      Lerp(u, cp[2], cp[3])};
    Point3f cp2[2] = {Lerp(u, cp1[0], cp1[1]), Lerp(u, cp1[1], cp1[2])};
    if (deriv) {
		// 下面是求积分的简写（算出来是一致的)
        if ((cp2[1] - cp2[0]).LengthSquared() > 0)
            *deriv = 3 * (cp2[1] - cp2[0]);
        else {
            // For a cubic Bezier, if the first three control points (say) are
            // coincident, then the derivative of the curve is legitimately (0,0,0)
            // at u=0.  This is problematic for us, though, since we'd like to be
            // able to compute a surface normal there.  In that case, just punt and
            // take the difference between the first and last control points, which
            // ain't great, but will hopefully do.
			// 个人理解, 这里的积分就是斜率, 当贝塞尔表示一条直线是, 斜率就是这个
            *deriv = cp[3] - cp[0];
        }
    }
    return Lerp(u, cp2[0], cp2[1]);
}

// Curve Method Definitions
CurveCommon::CurveCommon(const Point3f c[4], Float width0, Float width1,
                         CurveType type, const Normal3f *norm)
    : type(type) {
    width[0] = width0;
    width[1] = width1;
    for (int i = 0; i < 4; ++i)
        cpObj[i] = c[i];
    if (norm) {
        n[0] = Normalize(norm[0]);
        n[1] = Normalize(norm[1]);
        normalAngle = std::acos(Clamp(Dot(n[0], n[1]), 0, 1));
		// 这个 invSin 是预计算的
        invSinNormalAngle = 1 / std::sin(normalAngle);
    }
    ++nCurves;
}

std::vector<std::shared_ptr<Shape>> CreateCurve(
    const Transform *o2w, const Transform *w2o, bool reverseOrientation,
    const Point3f *c, Float w0, Float w1, CurveType type,
    const Normal3f *norm, int splitDepth) {
    std::vector<std::shared_ptr<Shape>> segments;
    std::shared_ptr<CurveCommon> common =
        std::make_shared<CurveCommon>(c, w0, w1, type, norm);
    const int nSegments = 1 << splitDepth;
    segments.reserve(nSegments);
    for (int i = 0; i < nSegments; ++i) {
        Float uMin = i / (Float)nSegments;
        Float uMax = (i + 1) / (Float)nSegments;
        segments.push_back(std::make_shared<Curve>(o2w, w2o, reverseOrientation,
                                                   common, uMin, uMax));
        ++nSplitCurves;
    }
    curveBytes += sizeof(CurveCommon) + nSegments * sizeof(Curve);
    return segments;
}

// 计算 贝塞尔曲线 的 Bound
// 这是基于 一段区间在 uMin uMax 的 bound
Bounds3f Curve::ObjectBound() const {
    // Compute object-space control points for curve segment, _cpObj_
	// 贝塞尔曲线的凸包性质, 可以计算出一段贝塞尔曲线的凸包，而贝塞尔曲线就在这个凸包中(凸包定义参考 计算几何内容)
    Point3f cpObj[4];
	// 算，对于 uMin,uMax 这一段距离中, 它的 4个 控制点(控制点组成凸包
    cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
    cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
    cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
    cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);
    Bounds3f b =
        Union(Bounds3f(cpObj[0], cpObj[1]), Bounds3f(cpObj[2], cpObj[3]));
	// 因为我们有宽度的定义, 所以需要做一下边缘拓展
    Float width[2] = {Lerp(uMin, common->width[0], common->width[1]),
                      Lerp(uMax, common->width[0], common->width[1])};
    return Expand(b, std::max(width[0], width[1]) * 0.5f);
}

bool Curve::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect,
                      bool testAlphaTexture) const {
    ProfilePhase p(isect ? Prof::CurveIntersect : Prof::CurveIntersectP);
    ++nTests;
    // Transform _Ray_ to object space
	// 第一步 把射线转化到 Object 空间
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute object-space control points for curve segment, _cpObj_
	// 计算这一段曲线 uMin-uMan 的4个控制点
    Point3f cpObj[4];
    cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
    cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
    cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
    cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);

    // Project curve control points to plane perpendicular to ray

    // Be careful to set the "up" direction passed to LookAt() to equal the
    // vector from the first to the last control points.  In turn, this
    // helps orient the curve to be roughly parallel to the x axis in the
    // ray coordinate system.
    //
    // In turn (especially for curves that are approaching stright lines),
    // we get curve bounds with minimal extent in y, which in turn lets us
    // early out more quickly in recursiveIntersect().
    Vector3f dx = Cross(ray.d, cpObj[3] - cpObj[0]);
	// dx 默认是垂直于 ray.d 和 cpObj[3] - cpObj[0] 的, 这是 Up 向量
    if (dx.LengthSquared() == 0) {
        // If the ray and the vector between the first and last control
        // points are parallel, dx will be zero.  Generate an arbitrary xy
        // orientation for the ray coordinate system so that intersection
        // tests can proceeed in this unusual case.
		// 如果 cpObj[3] - cpObj[0] 和 ray.d 平行的话, 那么我们随便选择一个 dy, 再构建一次
		// 总之目的是得到 dx
        Vector3f dy;
        CoordinateSystem(ray.d, &dx, &dy);
    }

	// 构建 worldToCamera
    Transform objectToRay = LookAt(ray.o, ray.o + ray.d, dx);
	// 映射控制点
    Point3f cp[4] = {objectToRay(cpObj[0]), objectToRay(cpObj[1]),
                     objectToRay(cpObj[2]), objectToRay(cpObj[3])};

    // Before going any further, see if the ray's bounding box intersects
    // the curve's bounding box. We start with the y dimension, since the y
    // extent is generally the smallest (and is often tiny) due to our
    // careful orientation of the ray coordinate ysstem above.
	// 这里是优化性的计算, 因为我们求交点, 可以先求 Bound + 拓展 和 射线的交点
	// 这里的 0 有关的判断, 是判断是否零点在 XOY 的投影内
    Float maxWidth = std::max(Lerp(uMin, common->width[0], common->width[1]),
                              Lerp(uMax, common->width[0], common->width[1]));
    if (std::max(std::max(cp[0].y, cp[1].y), std::max(cp[2].y, cp[3].y)) +
            0.5f * maxWidth < 0 ||
        std::min(std::min(cp[0].y, cp[1].y), std::min(cp[2].y, cp[3].y)) -
            0.5f * maxWidth > 0)
        return false;

    // Check for non-overlap in x.
    if (std::max(std::max(cp[0].x, cp[1].x), std::max(cp[2].x, cp[3].x)) +
            0.5f * maxWidth < 0 ||
        std::min(std::min(cp[0].x, cp[1].x), std::min(cp[2].x, cp[3].x)) -
            0.5f * maxWidth > 0)
        return false;

    // Check for non-overlap in z.
	// Z 的判断, 就在 Z 轴上
    Float rayLength = ray.d.Length();
    Float zMax = rayLength * ray.tMax;
    if (std::max(std::max(cp[0].z, cp[1].z), std::max(cp[2].z, cp[3].z)) +
            0.5f * maxWidth < 0 ||
        std::min(std::min(cp[0].z, cp[1].z), std::min(cp[2].z, cp[3].z)) -
            0.5f * maxWidth > zMax)
        return false;

	// 求曲线对射线的交点, 整体的思路是细分, 细分到每一小段曲线都接近线性，这样就可以用线性的高效方式来实现求交
	// ???? 这一段的实现原理还没懂, 目的就是动态的求出一个 最合适的 将这一段曲线分几段的算法
    // Compute refinement depth for curve, _maxDepth_
    Float L0 = 0;
    for (int i = 0; i < 2; ++i)
        L0 = std::max(
            L0, std::max(
                    std::max(std::abs(cp[i].x - 2 * cp[i + 1].x + cp[i + 2].x),
                             std::abs(cp[i].y - 2 * cp[i + 1].y + cp[i + 2].y)),
                    std::abs(cp[i].z - 2 * cp[i + 1].z + cp[i + 2].z)));

    Float eps =
        std::max(common->width[0], common->width[1]) * .05f;  // width / 20
    auto Log2 = [](float v) -> int {
        if (v < 1) return 0;
        uint32_t bits = FloatToBits(v);
        // https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
        // (With an additional add so get round-to-nearest rather than
        // round down.)
        return (bits >> 23) - 127 + (bits & (1 << 22) ? 1 : 0);
    };
    // Compute log base 4 by dividing log2 in half.
    int r0 = Log2(1.41421356237f * 6.f * L0 / (8.f * eps)) / 2;
    int maxDepth = Clamp(r0, 0, 10);
    ReportValue(refinementLevel, maxDepth);

    return recursiveIntersect(ray, tHit, isect, cp, Inverse(objectToRay), uMin,
                              uMax, maxDepth);
}

//给定线段, 给定曲线, 求交
bool Curve::recursiveIntersect(const Ray &ray, Float *tHit,
                               SurfaceInteraction *isect, const Point3f cp[4],
                               const Transform &rayToObject, Float u0, Float u1,
                               int depth) const {
    Float rayLength = ray.d.Length();

    if (depth > 0) {
        // Split curve segment into sub-segments and test for intersection
        Point3f cpSplit[7];
        SubdivideBezier(cp, cpSplit);

        // For each of the two segments, see if the ray's bounding box
        // overlaps the segment before recursively checking for
        // intersection with it.
        bool hit = false;
		// 2 分发, 这里先测 bounding
        Float u[3] = {u0, (u0 + u1) / 2.f, u1};
        // Pointer to the 4 control poitns for the current segment.
        const Point3f *cps = cpSplit;
        for (int seg = 0; seg < 2; ++seg, cps += 3) {
            Float maxWidth =
                std::max(Lerp(u[seg], common->width[0], common->width[1]),
                         Lerp(u[seg + 1], common->width[0], common->width[1]));

            // As above, check y first, since it most commonly lets us exit
            // out early.
			// 先测Y是因为Y容易第一个退出
            if (std::max(std::max(cps[0].y, cps[1].y),
                         std::max(cps[2].y, cps[3].y)) +
                        0.5 * maxWidth < 0 ||
                std::min(std::min(cps[0].y, cps[1].y),
                         std::min(cps[2].y, cps[3].y)) -
                        0.5 * maxWidth > 0)
                continue;

            if (std::max(std::max(cps[0].x, cps[1].x),
                         std::max(cps[2].x, cps[3].x)) +
                        0.5 * maxWidth < 0 ||
                std::min(std::min(cps[0].x, cps[1].x),
                         std::min(cps[2].x, cps[3].x)) -
                        0.5 * maxWidth > 0)
                continue;

            Float zMax = rayLength * ray.tMax;
            if (std::max(std::max(cps[0].z, cps[1].z),
                         std::max(cps[2].z, cps[3].z)) +
                        0.5 * maxWidth < 0 ||
                std::min(std::min(cps[0].z, cps[1].z),
                         std::min(cps[2].z, cps[3].z)) -
                        0.5 * maxWidth > zMax)
                continue;

			// 如果和 bounding 相交, 就继续递归测试
            hit |= recursiveIntersect(ray, tHit, isect, cps, rayToObject,
                                      u[seg], u[seg + 1], depth - 1);
            // If we found an intersection and this is a shadow ray,
            // we can exit out immediately.
			// 如果得到了答案, 就立刻退出
            if (hit && !tHit) return true;
        }
        return hit;
    } else {
        // Intersect ray with curve segment

        // Test ray against segment endpoint boundaries

		// 边缘判断逻辑
		// 这里的逻辑是这样的，通过做垂线判断和叉乘求积，得到原点 0,0 在控制点形成的连线的 右手侧还是左手侧
		// 注意这里使用的是 左手系
        // Test sample point against tangent perpendicular at curve start
        Float edge =
            (cp[1].y - cp[0].y) * -cp[0].y + cp[0].x * (cp[0].x - cp[1].x);
        if (edge < 0) return false;

        // Test sample point against tangent perpendicular at curve end
		// 同理
        edge = (cp[2].y - cp[3].y) * -cp[3].y + cp[3].x * (cp[3].x - cp[2].x);
        if (edge < 0) return false;

		// 这里计算出的，是 u0-u1 的比例 所以计算的点, 都是 Cruve 线上的点，也就是 Cruve 模型的中心点
        // Compute line $w$ that gives minimum distance to sample point
		// 在经过细分之后，直接求 原点(0,0) 在直线 CP[3]-CP[0] 上的投影的比例即可求出 交点
        Vector2f segmentDirection = Point2f(cp[3]) - Point2f(cp[0]);
        Float denom = segmentDirection.LengthSquared();
        if (denom == 0) return false;
        Float w = Dot(-Vector2f(cp[0]), segmentDirection) / denom;

        // Compute $u$ coordinate of curve intersection point and _hitWidth_
		// 通过 w, 计算出 u 和 hitWidth
        Float u = Clamp(Lerp(w, u0, u1), u0, u1);
        Float hitWidth = Lerp(u, common->width[0], common->width[1]);
        Normal3f nHit;
        if (common->type == CurveType::Ribbon) {
			// 如果是 ribbon 类型, 我们要使用 球形插值 来计算
            // Scale _hitWidth_ based on ribbon orientation
            Float sin0 = std::sin((1 - u) * common->normalAngle) *
                         common->invSinNormalAngle;
            Float sin1 =
                std::sin(u * common->normalAngle) * common->invSinNormalAngle;
			// 命中点 处的法线
            nHit = sin0 * common->n[0] + sin1 * common->n[1];
			// 那么 就要乘上 法线在这个方向 和 d 方向的 Cos
            hitWidth *= AbsDot(nHit, ray.d) / rayLength;
        }

        // Test intersection point against curve width
		// 验证性测试, 使用 w, 带入贝塞尔计算中, 进行一次验证, 顺便算关于 u 的积分 (也是那一点的切线方向
        Vector3f dpcdw;
        Point3f pc = EvalBezier(cp, Clamp(w, 0, 1), &dpcdw);
        Float ptCurveDist2 = pc.x * pc.x + pc.y * pc.y;
		// 宽度不够, 排除
        if (ptCurveDist2 > hitWidth * hitWidth * .25) return false;
        Float zMax = rayLength * ray.tMax;
		// 落在了 z的外面
        if (pc.z < 0 || pc.z > zMax) return false;

        // Compute $v$ coordinate of curve intersection point
		// Cruve 中心线和 Ray.d 他们的交点，之间的距离
        Float ptCurveDist = std::sqrt(ptCurveDist2);
		// 用切线判断是在 左还是右
        Float edgeFunc = dpcdw.x * -pc.y + pc.x * dpcdw.y;
		// v 表示的是宽度的 [0-1] 的取值点, 用 0.5是Cruve中心线 用 真实交点宽度/Cruve宽度 作为偏移
        Float v = (edgeFunc > 0) ? 0.5f + ptCurveDist / hitWidth
                                 : 0.5f - ptCurveDist / hitWidth;

        // Compute hit _t_ and partial derivatives for curve intersection
        if (tHit != nullptr) {
            // FIXME: this tHit isn't quite right for ribbons...
            *tHit = pc.z / rayLength;
            // Compute error bounds for curve intersection
            Vector3f pError(2 * hitWidth, 2 * hitWidth, 2 * hitWidth);

            // Compute $\dpdu$ and $\dpdv$ for curve intersection
            Vector3f dpdu, dpdv;
			// 用u计算正确的点, 直接上 common 中的大线段, u 也是 大线段中的值, 顺便求得 切线
			// 这里是物体的坐标系，也就是 Cruve 的坐标系！
            EvalBezier(common->cpObj, u, &dpdu);
            CHECK_NE(Vector3f(0, 0, 0), dpdu) << "u = " << u << ", cp = " <<
                common->cpObj[0] << ", " << common->cpObj[1] << ", " <<
                common->cpObj[2] << ", " << common->cpObj[3];

			// dv 的算法不一样
			if (common->type == CurveType::Ribbon)
			{
				// Ribbon 因为是面片形式, 所以 dv 直接用之前算出的 法线, 算个垂直即可
				// 切线的长度, 直接等于 Cruve 的宽度
				// 当然 dpdu 的长度, 因为有值带入，其长度就是那一段 近似直线 的Cruve 长度
				dpdv = Normalize(Cross(nHit, dpdu)) * hitWidth;
			}
			else {
                // Compute curve $\dpdv$ for flat and cylinder curves
				// 先转换到 Ray 坐标系
                Vector3f dpduPlane = (Inverse(rayToObject))(dpdu);
				// 先算出平切面的 切线 - Flat 类型就是这个值, 因为 Flat 的性质就是面向 射线
                Vector3f dpdvPlane =
                    Normalize(Vector3f(-dpduPlane.y, dpduPlane.x, 0)) *
                    hitWidth;
                if (common->type == CurveType::Cylinder) {
					// 这里是比较难想象的一点, 假设是一根圆柱，在射线坐标系上相交
					// 是用 v 的插值来模拟旋转角度
					// 但其实这里，使用弧度插值，会更好
                    // Rotate _dpdvPlane_ to give cylindrical appearance
                    Float theta = Lerp(v, -90., 90.);
					// 这个rot 是取逆矩阵
                    Transform rot = Rotate(-theta, dpduPlane);
                    dpdvPlane = rot(dpdvPlane);
                }
                dpdv = rayToObject(dpdvPlane);
            }
            *isect = (*ObjectToWorld)(SurfaceInteraction(
                ray(*tHit), pError, Point2f(u, v), -ray.d, dpdu, dpdv,
                Normal3f(0, 0, 0), Normal3f(0, 0, 0), ray.time, this));
        }
        ++nHits;
        return true;
    }
}

Float Curve::Area() const {
    // Compute object-space control points for curve segment, _cpObj_
    Point3f cpObj[4];
    cpObj[0] = BlossomBezier(common->cpObj, uMin, uMin, uMin);
    cpObj[1] = BlossomBezier(common->cpObj, uMin, uMin, uMax);
    cpObj[2] = BlossomBezier(common->cpObj, uMin, uMax, uMax);
    cpObj[3] = BlossomBezier(common->cpObj, uMax, uMax, uMax);
    Float width0 = Lerp(uMin, common->width[0], common->width[1]);
    Float width1 = Lerp(uMax, common->width[0], common->width[1]);
    Float avgWidth = (width0 + width1) * 0.5f;
    Float approxLength = 0.f;
    for (int i = 0; i < 3; ++i)
        approxLength += Distance(cpObj[i], cpObj[i + 1]);
    return approxLength * avgWidth;
}

Interaction Curve::Sample(const Point2f &u, Float *pdf) const {
    LOG(FATAL) << "Curve::Sample not implemented.";
    return Interaction();
}

std::vector<std::shared_ptr<Shape>> CreateCurveShape(const Transform *o2w,
                                                     const Transform *w2o,
                                                     bool reverseOrientation,
                                                     const ParamSet &params) {
    Float width = params.FindOneFloat("width", 1.f);
    Float width0 = params.FindOneFloat("width0", width);
    Float width1 = params.FindOneFloat("width1", width);

    int degree = params.FindOneInt("degree", 3);
    if (degree != 2 && degree != 3) {
        Error("Invalid degree %d: only degree 2 and 3 curves are supported.",
              degree);
        return {};
    }

    std::string basis = params.FindOneString("basis", "bezier");
    if (basis != "bezier" && basis != "bspline") {
        Error("Invalid basis \"%s\": only \"bezier\" and \"bspline\" are "
              "supported.", basis.c_str());
        return {};
    }

    int ncp;
    const Point3f *cp = params.FindPoint3f("P", &ncp);
    int nSegments;
    if (basis == "bezier") {
        // After the first segment, which uses degree+1 control points,
        // subsequent segments reuse the last control point of the previous
        // one and then use degree more control points.
        if (((ncp - 1 - degree) % degree) != 0) {
            Error("Invalid number of control points %d: for the degree %d "
                  "Bezier basis %d + n * %d are required, for n >= 0.", ncp,
                  degree, degree + 1, degree);
            return {};
        }
        nSegments = (ncp - 1) / degree;
    } else {
        if (ncp < degree + 1) {
            Error("Invalid number of control points %d: for the degree %d "
                  "b-spline basis, must have >= %d.", ncp, degree, degree + 1);
            return {};
        }
        nSegments = ncp - degree;
    }


    CurveType type;
    std::string curveType = params.FindOneString("type", "flat");
    if (curveType == "flat")
        type = CurveType::Flat;
    else if (curveType == "ribbon")
        type = CurveType::Ribbon;
    else if (curveType == "cylinder")
        type = CurveType::Cylinder;
    else {
        Error("Unknown curve type \"%s\".  Using \"cylinder\".", curveType.c_str());
        type = CurveType::Cylinder;
    }

    int nnorm;
    const Normal3f *n = params.FindNormal3f("N", &nnorm);
    if (n != nullptr) {
        if (type != CurveType::Ribbon) {
            Warning("Curve normals are only used with \"ribbon\" type curves.");
            n = nullptr;
        } else if (nnorm != nSegments + 1) {
            Error(
                "Invalid number of normals %d: must provide %d normals for ribbon "
                "curves with %d segments.", nnorm, nSegments + 1, nSegments);
            return {};
        }
    } else if (type == CurveType::Ribbon) {
        Error(
            "Must provide normals \"N\" at curve endpoints with ribbon "
            "curves.");
        return {};
    }

    int sd = params.FindOneInt("splitdepth",
                               int(params.FindOneFloat("splitdepth", 3)));

    std::vector<std::shared_ptr<Shape>> curves;
    // Pointer to the first control point for the current segment. This is
    // updated after each loop iteration depending on the current basis.
    const Point3f *cpBase = cp;
    for (int seg = 0; seg < nSegments; ++seg) {
        Point3f segCpBezier[4];

        // First, compute the cubic Bezier control points for the current
        // segment and store them in segCpBezier. (It is admittedly
        // wasteful storage-wise to turn b-splines into Bezier segments and
        // wasteful computationally to turn quadratic curves into cubics,
        // but yolo.)
        if (basis == "bezier") {
            if (degree == 2) {
                // Elevate to degree 3.
                segCpBezier[0] = cpBase[0];
                segCpBezier[1] = Lerp(2.f/3.f, cpBase[0], cpBase[1]);
                segCpBezier[2] = Lerp(1.f/3.f, cpBase[1], cpBase[2]);
                segCpBezier[3] = cpBase[2];
            } else {
                // Allset.
                for (int i = 0; i < 4; ++i)
                    segCpBezier[i] = cpBase[i];
            }
            cpBase += degree;
        } else {
            // Uniform b-spline.
            if (degree == 2) {
                // First compute equivalent Bezier control points via some
                // blossiming.  We have three control points and a uniform
                // knot vector; we'll label the points p01, p12, and p23.
                // We want the Bezier control points of the equivalent
                // curve, which are p11, p12, and p22.
                Point3f p01 = cpBase[0];
                Point3f p12 = cpBase[1];
                Point3f p23 = cpBase[2];

                // We already have p12.
                Point3f p11 = Lerp(0.5, p01, p12);
                Point3f p22 = Lerp(0.5, p12, p23);

                // Now elevate to degree 3.
                segCpBezier[0] = p11;
                segCpBezier[1] = Lerp(2.f/3.f, p11, p12);
                segCpBezier[2] = Lerp(1.f/3.f, p12, p22);
                segCpBezier[3] = p22;
            } else {
                // Otherwise we will blossom from p012, p123, p234, and p345
                // to the Bezier control points p222, p223, p233, and p333.
                // https://people.eecs.berkeley.edu/~sequin/CS284/IMGS/cubicbsplinepoints.gif
                Point3f p012 = cpBase[0];
                Point3f p123 = cpBase[1];
                Point3f p234 = cpBase[2];
                Point3f p345 = cpBase[3];

                Point3f p122 = Lerp(2.f/3.f, p012, p123);
                Point3f p223 = Lerp(1.f/3.f, p123, p234);
                Point3f p233 = Lerp(2.f/3.f, p123, p234);
                Point3f p334 = Lerp(1.f/3.f, p234, p345);

                Point3f p222 = Lerp(0.5f, p122, p223);
                Point3f p333 = Lerp(0.5f, p233, p334);

                segCpBezier[0] = p222;
                segCpBezier[1] = p223;
                segCpBezier[2] = p233;
                segCpBezier[3] = p333;
            }
            ++cpBase;
        }

        auto c = CreateCurve(o2w, w2o, reverseOrientation, segCpBezier,
                             Lerp(Float(seg) / Float(nSegments), width0, width1),
                             Lerp(Float(seg + 1) / Float(nSegments), width0, width1),
                             type, n ? &n[seg] : nullptr, sd);
        curves.insert(curves.end(), c.begin(), c.end());
    }
    return curves;
}

}  // namespace pbrt
