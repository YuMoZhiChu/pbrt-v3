
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

#ifndef PBRT_CORE_LOWDISCREPANCY_H
#define PBRT_CORE_LOWDISCREPANCY_H

// core/lowdiscrepancy.h*
#include "pbrt.h"
#include "rng.h"
#include "sampling.h"
#include "sobolmatrices.h"

namespace pbrt {

// Low Discrepancy Declarations
Float RadicalInverse(int baseIndex, uint64_t a);
std::vector<uint16_t> ComputeRadicalInversePermutations(RNG &rng);
static PBRT_CONSTEXPR int PrimeTableSize = 1000;
extern const int Primes[PrimeTableSize];
Float ScrambledRadicalInverse(int baseIndex, uint64_t a, const uint16_t *perm);
extern const int PrimeSums[PrimeTableSize];
inline void Sobol2D(int nSamplesPerPixelSample, int nPixelSamples,
                    Point2f *samples, RNG &rng);
extern uint32_t CMaxMinDist[17][32];
inline uint64_t SobolIntervalToIndex(const uint32_t log2Resolution,
                                     uint64_t sampleNum, const Point2i &p);
inline float SobolSampleFloat(int64_t index, int dimension,
                              uint32_t scramble = 0);
inline double SobolSampleDouble(int64_t index, int dimension,
                                uint64_t scramble = 0);

// Low Discrepancy Inline Functions
inline uint32_t ReverseBits32(uint32_t n) {
    n = (n << 16) | (n >> 16); // 因为 n 是32位的，所以 n << 16 是会自动拓展成 64 位，能够保留数据，64 位不能这么做的原因是，不能再往前拓展了
    n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
    n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
    n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2); // 这里是把 abcd 改成了 cdab
    n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1); // 这里要把 cdab 改成 dcba，所以用的是 5(0101) 和 a(1010)
    return n;
}

inline uint64_t ReverseBits64(uint64_t n) {
    uint64_t n0 = ReverseBits32((uint32_t)n); // 截断, 取后面 32位
    uint64_t n1 = ReverseBits32((uint32_t)(n >> 32)); // 取前面 32位
    return (n0 << 32) | n1;
}

// 整数 基于 base，的翻转接口
// 这里都是整数运算，并且 1234 和 123400 都会转换成 4321，多于的 0 会被舍弃
template <int base>
inline uint64_t InverseRadicalInverse(uint64_t inverse, int nDigits) {
    uint64_t index = 0;
    for (int i = 0; i < nDigits; ++i) {
        uint64_t digit = inverse % base;
        inverse /= base;
        index = index * base + digit;
    }
    return index;
}

// 02序列的算法函数
// C矩阵的列向量，在参数 a 的每一位是否位 1 的情况下，判断是否做 相加再模2(这个操作，等同于 或)
inline uint32_t MultiplyGenerator(const uint32_t *C, uint32_t a) {
    uint32_t v = 0;
    for (int i = 0; a != 0; ++i, a >>= 1)
        if (a & 1) v ^= C[i];
    return v;
}

// 02序列用到的算法函数
// 用 异或 一个特质的数 scramble 来表示是否做 对称交换
inline Float SampleGeneratorMatrix(const uint32_t *C, uint32_t a,
                                   uint32_t scramble = 0) {
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
    return std::min((MultiplyGenerator(C, a) ^ scramble) * Float(2.3283064365386963e-10),
                    OneMinusEpsilon);
#else
    return std::min((MultiplyGenerator(C, a) ^ scramble) * Float(0x1p-32),
                    OneMinusEpsilon);
#endif
}

// 一个数，获得它对应的 格雷码
// 这个坑爹的函数，全局引用就这里有个定义，就是拿给你看的
inline uint32_t GrayCode(uint32_t v) { return (v >> 1) ^ v; }

// 非常关键的一点，这个算法，不是 02序列，就是简单的 2进制 Radical Inverse,也叫做 VanDerCorput
// 在这里， 同样有 C 矩阵的概念，不过是被大大优化了的，优化的方式主要是 异或计算 和 格雷码 这两点
// C: 2进制 逆转矩阵
// n: 传入的计算数量，这里的 n 是被补齐到 2^i 次方的
// scramble: 随机数
// p: 存放结果的地方
inline void GrayCodeSample(const uint32_t *C, uint32_t n, uint32_t scramble,
                           Float *p) {
    uint32_t v = scramble; // 这是一个随机数，随机初值
    for (uint32_t i = 0; i < n; ++i) {
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
        p[i] = std::min(v * Float(2.3283064365386963e-10) /* 1/2^32 */,
                        OneMinusEpsilon);
#else
		// 第一个数，就直接存随机数
        p[i] = std::min(v * Float(0x1p-32) /* 1/2^32 */,
                        OneMinusEpsilon);
#endif
		// 算尾巴有几个 0 ，从 1 开始，到 n 本身(也就是 2^xxx）不从 0 开始的原因是，0 的 CountTrailingZeros 是搞不了的
		// 1. 因为用的是格雷码的理念，所以当前这个数v(这里把它当做格雷码来用)，只在 CountTrailingZeros(i + 1) 位上，和前一个数不同
		// 2. 因为我们用的是 2进制的 Radical Inverse 也就是 VanDerCorput，所以，要用 C 矩阵做一个映射，映射到对应相对称的位置
		// 3. 因为不一样，所以做一次 异或
        v ^= C[CountTrailingZeros(i + 1)];
    }
}

inline void GrayCodeSample(const uint32_t *C0, const uint32_t *C1, uint32_t n,
                           const Point2i &scramble, Point2f *p) {
    uint32_t v[2] = {(uint32_t)scramble.x, (uint32_t)scramble.y};
    for (uint32_t i = 0; i < n; ++i) {
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
        p[i].x = std::min(v[0] * Float(2.3283064365386963e-10), OneMinusEpsilon);
        p[i].y = std::min(v[1] * Float(2.3283064365386963e-10), OneMinusEpsilon);
#else
        p[i].x = std::min(v[0] * Float(0x1p-32), OneMinusEpsilon);
        p[i].y = std::min(v[1] * Float(0x1p-32), OneMinusEpsilon);
#endif
        v[0] ^= C0[CountTrailingZeros(i + 1)];
        v[1] ^= C1[CountTrailingZeros(i + 1)];
    }
}

inline void VanDerCorput(int nSamplesPerPixelSample, int nPixelSamples,
                         Float *samples, RNG &rng) {
    uint32_t scramble = rng.UniformUInt32();
    // Define _CVanDerCorput_ Generator Matrix
    const uint32_t CVanDerCorput[32] = {
#ifdef PBRT_HAVE_BINARY_CONSTANTS
      // clang-format off
	  // 这里是一个逆转 0 对应 最高位的 1,1 对应 次高位的 1，以此类推
      0b10000000000000000000000000000000,
      0b1000000000000000000000000000000,
      0b100000000000000000000000000000,
      0b10000000000000000000000000000,
      // Remainder of Van Der Corput generator matrix entries
      0b1000000000000000000000000000,
      0b100000000000000000000000000,
      0b10000000000000000000000000,
      0b1000000000000000000000000,
      0b100000000000000000000000,
      0b10000000000000000000000,
      0b1000000000000000000000,
      0b100000000000000000000,
      0b10000000000000000000,
      0b1000000000000000000,
      0b100000000000000000,
      0b10000000000000000,
      0b1000000000000000,
      0b100000000000000,
      0b10000000000000,
      0b1000000000000,
      0b100000000000,
      0b10000000000,
      0b1000000000,
      0b100000000,
      0b10000000,
      0b1000000,
      0b100000,
      0b10000,
      0b1000,
      0b100,
      0b10,
      0b1,
      // clang-format on
#else
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000,
        0x2000000,  0x1000000,  0x800000,   0x400000,   0x200000,  0x100000,
        0x80000,    0x40000,    0x20000,    0x10000,    0x8000,    0x4000,
        0x2000,     0x1000,     0x800,      0x400,      0x200,     0x100,
        0x80,       0x40,       0x20,       0x10,       0x8,       0x4,
        0x2,        0x1
#endif  // PBRT_HAVE_BINARY_CONSTANTS
    };
    int totalSamples = nSamplesPerPixelSample * nPixelSamples;
    GrayCodeSample(CVanDerCorput, totalSamples, scramble, samples);
    // Randomly shuffle 1D sample points
	// 第一次洗牌，是针对 SPP，每一个样本的采样数据结构，不是数组就是 1，是数组就是 数组内乱序
    for (int i = 0; i < nPixelSamples; ++i)
        Shuffle(samples + i * nSamplesPerPixelSample, nSamplesPerPixelSample, 1,
                rng);
	// 第二次洗牌，是针对 当前 Pixel 的所有样本，是数组就是 1，在 nPixelSamples 下乱序
	// 是数组就保存 nSamplesPerPixelSample 也就是数组的长度的记录
    Shuffle(samples, nPixelSamples, nSamplesPerPixelSample, rng);
}

inline void Sobol2D(int nSamplesPerPixelSample, int nPixelSamples,
                    Point2f *samples, RNG &rng) {
    Point2i scramble;
    scramble[0] = rng.UniformUInt32();
    scramble[1] = rng.UniformUInt32();

    // Define 2D Sobol$'$ generator matrices _CSobol[2]_
	// 2D 的使用这里给出的矩阵 (为什么是这个矩阵，不清楚 TODO????
    const uint32_t CSobol[2][32] = {
        {0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000,
         0x2000000, 0x1000000, 0x800000, 0x400000, 0x200000, 0x100000, 0x80000,
         0x40000, 0x20000, 0x10000, 0x8000, 0x4000, 0x2000, 0x1000, 0x800,
         0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1},
        {0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000,
         0xaa000000, 0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000,
         0x88880000, 0xcccc0000, 0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000,
         0xa000a000, 0xf000f000, 0x88008800, 0xcc00cc00, 0xaa00aa00, 0xff00ff00,
         0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0, 0x88888888, 0xcccccccc,
         0xaaaaaaaa, 0xffffffff}};
    GrayCodeSample(CSobol[0], CSobol[1], nSamplesPerPixelSample * nPixelSamples,
                   scramble, samples);
    for (int i = 0; i < nPixelSamples; ++i)
        Shuffle(samples + i * nSamplesPerPixelSample, nSamplesPerPixelSample, 1,
                rng);
    Shuffle(samples, nPixelSamples, nSamplesPerPixelSample, rng);
}

inline uint64_t SobolIntervalToIndex(const uint32_t m, uint64_t frame,
                                     const Point2i &p) {
    if (m == 0) return 0;

    const uint32_t m2 = m << 1;
    uint64_t index = uint64_t(frame) << m2;

    uint64_t delta = 0;
    for (int c = 0; frame; frame >>= 1, ++c)
        if (frame & 1)  // Add flipped column m + c + 1.
            delta ^= VdCSobolMatrices[m - 1][c];

    // flipped b
    uint64_t b = (((uint64_t)((uint32_t)p.x) << m) | ((uint32_t)p.y)) ^ delta;

    for (int c = 0; b; b >>= 1, ++c)
        if (b & 1)  // Add column 2 * m - c.
            index ^= VdCSobolMatricesInv[m - 1][c];

    return index;
}

inline Float SobolSample(int64_t index, int dimension, uint64_t scramble = 0) {
#ifdef PBRT_FLOAT_AS_DOUBLE
    return SobolSampleDouble(index, dimension, scramble);
#else
    return SobolSampleFloat(index, dimension, (uint32_t)scramble);
#endif
}

inline float SobolSampleFloat(int64_t a, int dimension, uint32_t scramble) {
    CHECK_LT(dimension, NumSobolDimensions) <<
        "Integrator has consumed too many Sobol' dimensions; you "
        "may want to use a Sampler without a dimension limit like "
        "\"02sequence.\"";
    uint32_t v = scramble;
    for (int i = dimension * SobolMatrixSize; a != 0; a >>= 1, i++)
        if (a & 1) v ^= SobolMatrices32[i];
#ifndef PBRT_HAVE_HEX_FP_CONSTANTS
    return std::min(v * 2.3283064365386963e-10f /* 1/2^32 */,
                    FloatOneMinusEpsilon);
#else
    return std::min(v * 0x1p-32f /* 1/2^32 */,
                    FloatOneMinusEpsilon);
#endif
}

inline double SobolSampleDouble(int64_t a, int dimension, uint64_t scramble) {
  CHECK_LT(dimension, NumSobolDimensions) <<
      "Integrator has consumed too many Sobol' dimensions; you "
      "may want to use a Sampler without a dimension limit like "
      "\"02sequence\".";
    uint64_t result = scramble & ~ - (1LL << SobolMatrixSize);
    for (int i = dimension * SobolMatrixSize; a != 0; a >>= 1, i++)
        if (a & 1) result ^= SobolMatrices64[i];
    return std::min(result * (1.0 / (1ULL << SobolMatrixSize)),
                    DoubleOneMinusEpsilon);
}

}  // namespace pbrt

#endif  // PBRT_CORE_LOWDISCREPANCY_H
