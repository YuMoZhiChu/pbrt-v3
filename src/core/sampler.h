
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

#ifndef PBRT_CORE_SAMPLER_H
#define PBRT_CORE_SAMPLER_H

// core/sampler.h*
#include "pbrt.h"
#include "geometry.h"
#include "rng.h"
#include <inttypes.h>

namespace pbrt {

// Sampler Declarations
// 采样基类，定义了接口，和一些功能性的函数实现
class Sampler {
  public:
    // Sampler Interface
    virtual ~Sampler();
	// Sampler 的初始化，一定要传入表示在这个像素区域，的采样样本个数，采样样本个数可以是 1 个(极少情况下)，可能有无数个，所以用 int64 来定义样本个数
    Sampler(int64_t samplesPerPixel);
	// 在开始渲染之前，先传入起点像素位置，用于一些算法可能需要用到该接口，来观察整个样本的分布情况
    virtual void StartPixel(const Point2i &p);
	// 返回当前样本向量，下一维度的浮点值
    virtual Float Get1D() = 0;
	// 返回当前样本向量，下两个维度的浮点值，两次 Get1D = Get2D
    virtual Point2f Get2D() = 0;
	// 便捷接口，获取照相机的采样内容
    CameraSample GetCameraSample(const Point2i &pRaster);
	// 申请将当前维度变成一个 一元的数组 [a, b, c, d, ...] 整个数组是一个维度的值，n 是数组长度
    void Request1DArray(int n);
	// 申请将当前维度变成一个 二元的数组 [Point2f(a, b), Point2f(c, d), ...] 整个数组是一个维度的值，n 是数组长度
    void Request2DArray(int n);
	// 调整 维度是数组时候，这个数组的长度，因为数组的长度会对部分算法的计算有帮助，比如保持2的幂次的长度
    virtual int RoundCount(int n) const { return n; }
	// 获取当前维度(该维度是一元数组) 的头指针，n 是长度
    const Float *Get1DArray(int n);
	// 获取当前维度(该维度是二元数组) 的头指针，n 是长度
    const Point2f *Get2DArray(int n);
	// 当完成一次采样时候调用，重新把维度设置为0，进行一次新的采样（对当前像素区域
	// 当完成了该像素区域的所有采样后，返回 True，否则返回 False
    virtual bool StartNextSample();
	// 采样器的克隆函数，对于带有一大堆当前采样状态的采样器，多线程显然不适合在采样器中实现
	// 通过传入随机数种子，生成不同的采样器，以便在多线程中运行 - 随机还有一个作用就是防止出现重复率高的噪声
    virtual std::unique_ptr<Sampler> Clone(int seed) = 0;
	// 对应奇特的采样算法，有一些算法可能会跳跃性的在某个像素区域周围采样，这里保证其采样的样本数量，小于等于初始定义的最大数量
    virtual bool SetSampleNumber(int64_t sampleNum);
    std::string StateString() const {
      return StringPrintf("(%d,%d), sample %" PRId64, currentPixel.x,
                          currentPixel.y, currentPixelSampleIndex);
    }
    int64_t CurrentSampleNumber() const { return currentPixelSampleIndex; }

    // Sampler Public Data
	// SPP 每一个像素区域的采样数
    const int64_t samplesPerPixel;

  protected:
    // Sampler Protected Data
	// 当前采样的像素区域的位置
    Point2i currentPixel;
	// 当前是给区域的第几个样本
    int64_t currentPixelSampleIndex;
	// 维度是数组时，大小记录
    std::vector<int> samples1DArraySizes, samples2DArraySizes;
	// 一维数据记录
	// 第一个系数是维度 第二个数组是 集中处理 即前 n 个, 是一个样本的一个维度的数据
    std::vector<std::vector<Float>> sampleArray1D;
	// 二维数据记录
    std::vector<std::vector<Point2f>> sampleArray2D;

  private:
    // Sampler Private Data
    size_t array1DOffset, array2DOffset;
};

class PixelSampler : public Sampler {
  public:
    // PixelSampler Public Methods
	// Pixel Sampler 会多传入一个最大维度，当超过这个最大维度时，返回一个归一化[0,1)的随机数
    PixelSampler(int64_t samplesPerPixel, int nSampledDimensions);
    bool StartNextSample();
    bool SetSampleNumber(int64_t);
    Float Get1D();
    Point2f Get2D();

  protected:
    // PixelSampler Protected Data
	// 用二次vector 来记录所有的 1D 和 2D 数据
	  // 第一个系数是维度 第二个数组是 就是哪个样本，对应下标就是哪个
    std::vector<std::vector<Float>> samples1D;
    std::vector<std::vector<Point2f>> samples2D;
	// 记录1D和2D 的当前偏移
    int current1DDimension = 0, current2DDimension = 0;
    RNG rng;
};

// Global 采样器的思想是
// 采样的样本 -> 采样的像素区域, 这是一个 多 -> 1 的映射
class GlobalSampler : public Sampler {
  public:
    // GlobalSampler Public Methods
    bool StartNextSample();
	// 为某个像素点，开始生成样本
    void StartPixel(const Point2i &);
    bool SetSampleNumber(int64_t sampleNum);
    Float Get1D();
    Point2f Get2D();
    GlobalSampler(int64_t samplesPerPixel) : Sampler(samplesPerPixel) {}
	// 根据当前的采样的像素区域，获取第 sampleNum 个样本，在样本总表的位置
    virtual int64_t GetIndexForSample(int64_t sampleNum) const = 0;
	// 获取在样本总表中，第 index 个样本，在样本原始数据中的，第 dimension 个维度的数据
    virtual Float SampleDimension(int64_t index, int dimension) const = 0;

  private:
    // GlobalSampler Private Data
	// 记录当前维度
    int dimension;
	// 记录当前样本，在样本表中的 index
    int64_t intervalSampleIndex;
	// 为照相机设计的前5个维度
    static const int arrayStartDim = 5;
    int arrayEndDim;
};

}  // namespace pbrt

#endif  // PBRT_CORE_SAMPLER_H
