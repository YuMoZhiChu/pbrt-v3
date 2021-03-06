﻿
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

#ifndef PBRT_CORE_REFLECTION_H
#define PBRT_CORE_REFLECTION_H

// core/reflection.h*
#include "pbrt.h"
#include "geometry.h"
#include "microfacet.h"
#include "shape.h"
#include "spectrum.h"

namespace pbrt {

// Reflection Declarations
Float FrDielectric(Float cosThetaI, Float etaI, Float etaT);
Spectrum FrConductor(Float cosThetaI, const Spectrum &etaI,
                     const Spectrum &etaT, const Spectrum &k);

// BSDF Inline Functions
// inline 函数，计算传入 w 和 BRDF 坐标系的快速计算，这里返回的是对应 球坐标系的 (theta,phi) 对应的值
// 快速计算 cos
inline Float CosTheta(const Vector3f &w) { return w.z; }
// 快速计算 cos 的 平方
inline Float Cos2Theta(const Vector3f &w) { return w.z * w.z; }
// 快速计算 cos 的 绝对值
inline Float AbsCosTheta(const Vector3f &w) { return std::abs(w.z); }
// 快速计算 sin 的 平方
inline Float Sin2Theta(const Vector3f &w) {
    return std::max((Float)0, (Float)1 - Cos2Theta(w));
}
// 快速计算 sin
inline Float SinTheta(const Vector3f &w) { return std::sqrt(Sin2Theta(w)); }
// 快速计算 tan
inline Float TanTheta(const Vector3f &w) { return SinTheta(w) / CosTheta(w); }
// 快速计算 tan 的 平方
inline Float Tan2Theta(const Vector3f &w) {
    return Sin2Theta(w) / Cos2Theta(w);
}

// 快速计算 cos phi
inline Float CosPhi(const Vector3f &w) {
    Float sinTheta = SinTheta(w);
	// 考量到 /0 的情况，那么我们默认为 0 时，sin 是 0
    return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
}

inline Float SinPhi(const Vector3f &w) {
    Float sinTheta = SinTheta(w);
    // 考量到 /0 的情况，那么我们默认为 0 时，cos 是 1
    return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
}

inline Float Cos2Phi(const Vector3f &w) { return CosPhi(w) * CosPhi(w); }

inline Float Sin2Phi(const Vector3f &w) { return SinPhi(w) * SinPhi(w); }

// 传入两条射线，快速计算他们的 cos delta phi
// 这里只计算 平面角的差值 的 cos，所以和 z 没关系
// 这个公式可推
inline Float CosDPhi(const Vector3f &wa, const Vector3f &wb) {
    return Clamp(
        (wa.x * wb.x + wa.y * wb.y) / std::sqrt((wa.x * wa.x + wa.y * wa.y) *
                                                (wb.x * wb.x + wb.y * wb.y)),
        -1, 1);
}

// 根据 wo，计算 wr 的方法
inline Vector3f Reflect(const Vector3f &wo, const Vector3f &n) {
    return -wo + 2 * Dot(wo, n) * n;
}

// 计算，折射，的出射角的计算流程
// wi 入射角
// n 跟入射角同一半球体内的法线
// eta 折射率
// wt 计算出的 折射方向
inline bool Refract(const Vector3f &wi, const Normal3f &n, Float eta,
                    Vector3f *wt) {
    // Compute $\cos \theta_\roman{t}$ using Snell's law
    Float cosThetaI = Dot(n, wi);
    Float sin2ThetaI = std::max(Float(0), Float(1 - cosThetaI * cosThetaI));
    Float sin2ThetaT = eta * eta * sin2ThetaI;

    // Handle total internal reflection for transmission
    if (sin2ThetaT >= 1) return false;
    Float cosThetaT = std::sqrt(1 - sin2ThetaT);
    *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
    return true;
}

inline bool SameHemisphere(const Vector3f &w, const Vector3f &wp) {
    return w.z * wp.z > 0;
}

inline bool SameHemisphere(const Vector3f &w, const Normal3f &wp) {
    return w.z * wp.z > 0;
}

// BSDF Declarations
enum BxDFType {
	// 反射
    BSDF_REFLECTION = 1 << 0,
    // 透射
	BSDF_TRANSMISSION = 1 << 1,
    // 漫反射
	BSDF_DIFFUSE = 1 << 2,
    // 高光，其中 逆向反射，也视为 高光反色
	BSDF_GLOSSY = 1 << 3,
    // 镜面反射
	BSDF_SPECULAR = 1 << 4,
    BSDF_ALL = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR | BSDF_REFLECTION |
               BSDF_TRANSMISSION,
};

struct FourierBSDFTable {
    // FourierBSDFTable Public Data
    Float eta;
    int mMax;
    int nChannels;
    int nMu;
    Float *mu;
    int *m;
    int *aOffset;
    Float *a;
    Float *a0;
    Float *cdf;
    Float *recip;

    // FourierBSDFTable Public Methods
    static bool Read(const std::string &filename, FourierBSDFTable *table);
    const Float *GetAk(int offsetI, int offsetO, int *mptr) const {
        *mptr = m[offsetO * nMu + offsetI];
        return a + aOffset[offsetO * nMu + offsetI];
    }
    bool GetWeightsAndOffset(Float cosTheta, int *offset,
                             Float weights[4]) const;
};

// BSDF 类，他的特点
// 1. 整合了 BRDF 和 BTDF
// 2. 对上层接口，隐藏了细节法线，比如 阴影法线，来自模型或贴图的法线
class BSDF {
  public:
    // BSDF Public Methods
	// BSDF 的构造，需要参数
	// si ： 表面的微分信息
	// eta : 表面的相对折射率，对于不透明的表面，不使用 eta，默认值是 1
	// 
    BSDF(const SurfaceInteraction &si, Float eta = 1)
        : eta(eta),
          ns(si.shading.n),
          ng(si.n),
          ss(Normalize(si.shading.dpdu)),
          ts(Cross(ns, ss)) {}
    void Add(BxDF *b) {
        CHECK_LT(nBxDFs, MaxBxDFs);
        bxdfs[nBxDFs++] = b;
    }
    int NumComponents(BxDFType flags = BSDF_ALL) const;
    Vector3f WorldToLocal(const Vector3f &v) const {
        return Vector3f(Dot(v, ss), Dot(v, ts), Dot(v, ns));
    }
    Vector3f LocalToWorld(const Vector3f &v) const {
        return Vector3f(ss.x * v.x + ts.x * v.y + ns.x * v.z,
                        ss.y * v.x + ts.y * v.y + ns.y * v.z,
                        ss.z * v.x + ts.z * v.y + ns.z * v.z);
    }
    Spectrum f(const Vector3f &woW, const Vector3f &wiW,
               BxDFType flags = BSDF_ALL) const;
    Spectrum rho(int nSamples, const Point2f *samples1, const Point2f *samples2,
                 BxDFType flags = BSDF_ALL) const;
    Spectrum rho(const Vector3f &wo, int nSamples, const Point2f *samples,
                 BxDFType flags = BSDF_ALL) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType type = BSDF_ALL,
                      BxDFType *sampledType = nullptr) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi,
              BxDFType flags = BSDF_ALL) const;
    std::string ToString() const;

    // BSDF Public Data
    const Float eta;

  private:
    // BSDF Private Methods
    ~BSDF() {}

    // BSDF Private Data
    const Normal3f ns, ng;
    const Vector3f ss, ts;
    int nBxDFs = 0;
    static PBRT_CONSTEXPR int MaxBxDFs = 8;
    BxDF *bxdfs[MaxBxDFs];
    friend class MixMaterial;
};

inline std::ostream &operator<<(std::ostream &os, const BSDF &bsdf) {
    os << bsdf.ToString();
    return os;
}

// BxDF Declarations
// BRDF 和 BTDF 共有的接口
class BxDF {
  public:
    // BxDF Interface
    virtual ~BxDF() {}
	// 初始化时，必须初始化它的类型
    BxDF(BxDFType type) : type(type) {}
	// 类型匹配，也是满足 xxx 类型
    bool MatchesFlags(BxDFType t) const { return (type & t) == type; }
	// 核心函数， ????
	// 返回对应的 入射出射 方向的光谱表达式
	// 不同波长的光是解耦的，一个波长的能量将不会在另一个波长处反射
	// 反射的效果，可以直接用 光谱 来表示
	// 个人理解，给定一个 入射方向，出射方向，不同的材质会有不同的光谱分布，所以我们返回一个光谱来表示这种分布
	// 这里的光谱分布，应该是指， 出射光 的 分布
	// 对于支持 荧光材料 的材质，该方法会要求返回一个 举证，表示 光谱之间能量转移之间的 编码
    virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const = 0;

	// 并不是所有的 BxDF 都可以用这种方法，来得到 光谱 分布
	// 比如完美镜面反射，比如 镜子，玻璃，光就只能从 单个入射反向 反射到 单个出射方向
	// 这一类的 BxDF 最好使用 delta distributions 为 0 来描述
	// 这个函数，处理：
	// 1. 用 delta distributions 描述的散射
	// 2. 从 BxDF 沿多个方向，做散射光 的 随机采样 (14.1 的 Monte Carlo 的 BSDF 采样
	//
	// 同理，该函数传入 wo 和 wi，并返回 BxDF 的光照 分布
	// 对于 delta distributions，BxDF 必须选择这种方式传入 入射光，因为 调用者 没有其他机会生成正确的 wi ????
	// sample 和 pdf 参数会在 14.1 中说明
    virtual Spectrum Sample_f(const Vector3f &wo, Vector3f *wi,
                              const Point2f &sample, Float *pdf,
                              BxDFType *sampledType = nullptr) const;
	// 求半球全反射
	// 对于使用 Monte Carlo 的 BxDF 算法，我们需要 nSamples 和 samples 这两个参数
    virtual Spectrum rho(const Vector3f &wo, int nSamples,
                         const Point2f *samples) const;
	// 半球-半球 反射
	// 可以不传入 wo,wi
	// 参数是 Monte Carlo 的 BxDF 算法 的需求
    virtual Spectrum rho(int nSamples, const Point2f *samples1,
                         const Point2f *samples2) const;
    virtual Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    virtual std::string ToString() const = 0;

    // BxDF Public Data
	// 保存 BxDF 的类型
    const BxDFType type;
};

inline std::ostream &operator<<(std::ostream &os, const BxDF &bxdf) {
    os << bxdf.ToString();
    return os;
}

// 在用于表示 混合材质时，MixMaterial
// 我们可以用 一个 BxDF* 和 Spectrum 来实现此功能
// 其他的接口，就是直接乘上该参数即可
class ScaledBxDF : public BxDF {
  public:
    // ScaledBxDF Public Methods
    ScaledBxDF(BxDF *bxdf, const Spectrum &scale)
        : BxDF(BxDFType(bxdf->type)), bxdf(bxdf), scale(scale) {}
    Spectrum rho(const Vector3f &w, int nSamples,
                 const Point2f *samples) const {
        return scale * bxdf->rho(w, nSamples, samples);
    }
    Spectrum rho(int nSamples, const Point2f *samples1,
                 const Point2f *samples2) const {
        return scale * bxdf->rho(nSamples, samples1, samples2);
    }
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    BxDF *bxdf;
	// 缩放参数
    Spectrum scale;
};

// 这是一个抽象的 菲尼尔类，提供了用于菲涅尔计算的接口
// 该抽象父类，可以简化后续 BRDF 的实现
class Fresnel {
  public:
    // Fresnel Interface
    virtual ~Fresnel();
	// 传入入射方向 和 表面法线 的余弦值，返回 表面的 反射光
    virtual Spectrum Evaluate(Float cosI) const = 0;
    virtual std::string ToString() const = 0;
};

inline std::ostream &operator<<(std::ostream &os, const Fresnel &f) {
    os << f.ToString();
    return os;
}

// 导体的 菲涅尔 计算
class FresnelConductor : public Fresnel {
  public:
    // FresnelConductor Public Methods
	// 调用的是 FrConductor，对 cosThetal 做了 abs 处理
    Spectrum Evaluate(Float cosThetaI) const;
    FresnelConductor(const Spectrum &etaI, const Spectrum &etaT,
                     const Spectrum &k)
        : etaI(etaI), etaT(etaT), k(k) {}
    std::string ToString() const;

  private:
	// 入射材质的折射率，折射材质的折射率，以及，导体吸收参数 k
    Spectrum etaI, etaT, k;
};

// 电介质的 菲涅尔 计算
class FresnelDielectric : public Fresnel {
  public:
    // FresnelDielectric Public Methods
	// 调用的是 FrDielectric
    Spectrum Evaluate(Float cosThetaI) const;
    FresnelDielectric(Float etaI, Float etaT) : etaI(etaI), etaT(etaT) {}
    std::string ToString() const;

  private:
    Float etaI, etaT;
};

class FresnelNoOp : public Fresnel {
  public:
	// 直接全反射
    Spectrum Evaluate(Float) const { return Spectrum(1.); }
    std::string ToString() const { return "[ FresnelNoOp ]"; }
};

// 计算镜面反射的类
class SpecularReflection : public BxDF {
  public:
    // SpecularReflection Public Methods
    SpecularReflection(const Spectrum &R, Fresnel *fresnel)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),
          R(R),
          fresnel(fresnel) {}
	// 这里返回的是 0，因为镜面反射，在任何方向，不会返回任何散射
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return Spectrum(0.f);
    }
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; }
    std::string ToString() const;

  private:
    // SpecularReflection Private Data
	// 光谱，用于缩放反射的衍射
    const Spectrum R;
	// 菲涅尔对象指针，用于描述电介质或导体的菲涅尔特性
    const Fresnel *fresnel;
};

// 计算镜面透射的类
class SpecularTransmission : public BxDF {
  public:
    // SpecularTransmission Public Methods
    SpecularTransmission(const Spectrum &T, Float etaA, Float etaB,
                         TransportMode mode)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
          T(T),
          etaA(etaA),
          etaB(etaB),
          fresnel(etaA, etaB),
          mode(mode) {}
	// 因为 BTDF 是 缩放的增量分布 scaled delta distribution，所以这里返回的是 0（也可以理解成，没有散射
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return Spectrum(0.f);
    }
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &sample,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; }
    std::string ToString() const;

  private:
    // SpecularTransmission Private Data
	// 缩放（比如自身颜色
    const Spectrum T;
	// etaA 是 above 表面的折射率
	// etaB 是 below 表面的折射率
    const Float etaA, etaB;
    // 电介质的菲涅尔计算，因为导体通常不透光，所以始终用 电介质也就是非导体 的计算方式
	const FresnelDielectric fresnel;
	// 计算 Ray 的交点时，是从 光源 开始，还是从 摄像机 开始
    const TransportMode mode;
};

// 为了使得 14,15,16 章中的某些 Monte Carlo 更有效率，所以提供了一个调制版的镜面版本 ????
class FresnelSpecular : public BxDF {
  public:
    // FresnelSpecular Public Methods
    FresnelSpecular(const Spectrum &R, const Spectrum &T, Float etaA,
                    Float etaB, TransportMode mode)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR)),
          R(R),
          T(T),
          etaA(etaA),
          etaB(etaB),
          mode(mode) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const {
        return Spectrum(0.f);
    }
	// 14.1.3 的内容
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const { return 0; }
    std::string ToString() const;

  private:
    // FresnelSpecular Private Data
    const Spectrum R, T;
    const Float etaA, etaB;
    const TransportMode mode;
};

class LambertianReflection : public BxDF {
  public:
    // LambertianReflection Public Methods
    LambertianReflection(const Spectrum &R)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum rho(const Vector3f &, int, const Point2f *) const { return R; }
    Spectrum rho(int, const Point2f *, const Point2f *) const { return R; }
    std::string ToString() const;

  private:
    // LambertianReflection Private Data
    const Spectrum R;
};

class LambertianTransmission : public BxDF {
  public:
    // LambertianTransmission Public Methods
    LambertianTransmission(const Spectrum &T)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_DIFFUSE)), T(T) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum rho(const Vector3f &, int, const Point2f *) const { return T; }
    Spectrum rho(int, const Point2f *, const Point2f *) const { return T; }
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // LambertianTransmission Private Data
    Spectrum T;
};

class OrenNayar : public BxDF {
  public:
    // OrenNayar Public Methods
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    OrenNayar(const Spectrum &R, Float sigma)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(R) {
        sigma = Radians(sigma);
        Float sigma2 = sigma * sigma;
        A = 1.f - (sigma2 / (2.f * (sigma2 + 0.33f)));
        B = 0.45f * sigma2 / (sigma2 + 0.09f);
    }
    std::string ToString() const;

  private:
    // OrenNayar Private Data
    const Spectrum R;
    Float A, B;
};

class MicrofacetReflection : public BxDF {
  public:
    // MicrofacetReflection Public Methods
    MicrofacetReflection(const Spectrum &R,
                         MicrofacetDistribution *distribution, Fresnel *fresnel)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
          R(R),
          distribution(distribution),
          fresnel(fresnel) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // MicrofacetReflection Private Data
    const Spectrum R;
    const MicrofacetDistribution *distribution;
    const Fresnel *fresnel;
};

class MicrofacetTransmission : public BxDF {
  public:
    // MicrofacetTransmission Public Methods
    MicrofacetTransmission(const Spectrum &T,
                           MicrofacetDistribution *distribution, Float etaA,
                           Float etaB, TransportMode mode)
        : BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_GLOSSY)),
          T(T),
          distribution(distribution),
          etaA(etaA),
          etaB(etaB),
          fresnel(etaA, etaB),
          mode(mode) {}
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // MicrofacetTransmission Private Data
    const Spectrum T;
    const MicrofacetDistribution *distribution;
    const Float etaA, etaB;
    const FresnelDielectric fresnel;
    const TransportMode mode;
};

class FresnelBlend : public BxDF {
  public:
    // FresnelBlend Public Methods
    FresnelBlend(const Spectrum &Rd, const Spectrum &Rs,
                 MicrofacetDistribution *distrib);
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    Spectrum SchlickFresnel(Float cosTheta) const {
        auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
        return Rs + pow5(1 - cosTheta) * (Spectrum(1.) - Rs);
    }
    Spectrum Sample_f(const Vector3f &wi, Vector3f *sampled_f, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // FresnelBlend Private Data
    const Spectrum Rd, Rs;
    MicrofacetDistribution *distribution;
};

class FourierBSDF : public BxDF {
  public:
    // FourierBSDF Public Methods
    Spectrum f(const Vector3f &wo, const Vector3f &wi) const;
    FourierBSDF(const FourierBSDFTable &bsdfTable, TransportMode mode)
        : BxDF(BxDFType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)),
          bsdfTable(bsdfTable),
          mode(mode) {}
    Spectrum Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                      Float *pdf, BxDFType *sampledType) const;
    Float Pdf(const Vector3f &wo, const Vector3f &wi) const;
    std::string ToString() const;

  private:
    // FourierBSDF Private Data
    const FourierBSDFTable &bsdfTable;
    const TransportMode mode;
};

// BSDF Inline Method Definitions
inline int BSDF::NumComponents(BxDFType flags) const {
    int num = 0;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags)) ++num;
    return num;
}

}  // namespace pbrt

#endif  // PBRT_CORE_REFLECTION_H
