
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


// samplers/halton.cpp*
#include "samplers/halton.h"
#include "paramset.h"
#include "rng.h"

namespace pbrt {

// HaltonSampler Local Constants
static PBRT_CONSTEXPR int kMaxResolution = 128;

// HaltonSampler Utility Functions
static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y);
static uint64_t multiplicativeInverse(int64_t a, int64_t n) {
    int64_t x, y;
    extendedGCD(a, n, &x, &y);
    return Mod(x, n);
}

static void extendedGCD(uint64_t a, uint64_t b, int64_t *x, int64_t *y) {
    if (b == 0) {
        *x = 1;
        *y = 0;
        return;
    }
    int64_t d = a / b, xp, yp;
    extendedGCD(b, a % b, &xp, &yp);
    *x = yp;
    *y = xp - (d * yp);
}

// HaltonSampler Method Definitions
HaltonSampler::HaltonSampler(int samplesPerPixel, const Bounds2i &sampleBounds,
                             bool sampleAtPixelCenter)
    : GlobalSampler(samplesPerPixel), sampleAtPixelCenter(sampleAtPixelCenter) {
    // Generate random digit permutations for Halton sampler
	// 第一次初始化的时候，创建一次 Permutation 数组
    if (radicalInversePermutations.empty()) {
        RNG rng;
        radicalInversePermutations = ComputeRadicalInversePermutations(rng);
    }

    // Find radical inverse base scales and exponents that cover sampling area
    Vector2i res = sampleBounds.pMax - sampleBounds.pMin;
	// 计算出，跟位置有关系的 base
    for (int i = 0; i < 2; ++i) {
        int base = (i == 0) ? 2 : 3;
        int scale = 1, exp = 0;
        while (scale < std::min(res[i], kMaxResolution)) {
            scale *= base;
            ++exp;
        }
        baseScales[i] = scale;
        baseExponents[i] = exp;
    }

    // Compute stride in samples for visiting each pixel area
    sampleStride = baseScales[0] * baseScales[1];

	// 计算欧几里得机制
	// 这里满足 a * x + b * y = 1，其中的 b 为负数的情况下，是对 a 取模
	// 满足一个特性 
    // Compute multiplicative inverses for _baseScales_
    multInverse[0] = multiplicativeInverse(baseScales[1], baseScales[0]);
    multInverse[1] = multiplicativeInverse(baseScales[0], baseScales[1]);
}

std::vector<uint16_t> HaltonSampler::radicalInversePermutations;
int64_t HaltonSampler::GetIndexForSample(int64_t sampleNum) const {
    if (currentPixel != pixelForOffset) {
        // Compute Halton sample offset for _currentPixel_
        offsetForCurrentPixel = 0;
        if (sampleStride > 1) {
			// 取模，在 2^j * 3^k 的范围内，做全局采样
            Point2i pm(Mod(currentPixel[0], kMaxResolution),
                       Mod(currentPixel[1], kMaxResolution));
            for (int i = 0; i < 2; ++i) {
				// 这里做 2,3 的基逆，就是为了把对应像素位置 -> 转换为 halton 随机点，然后再构造大的随机数 Index
                uint64_t dimOffset =
                    (i == 0)
                        ? InverseRadicalInverse<2>(pm[i], baseExponents[i])
                        : InverseRadicalInverse<3>(pm[i], baseExponents[i]);
				// 欧几里得的乘法，是满足一个特性的
				// offset0 * a * x + offset1 * b * y = result
				// result %= (a * b)
				// 这个 result 依然满足 result % a = offset0，result % b == offset1
				// 所以是 offset0，offset1 会在每隔 a*b 的区间内，都能一一对应一个独立的数
				// 所以这个 result ，可以作为 1. 对应的采样值 2. 因为其可以由每隔 a*b 对应生成，所以满足我们的复合采样需求
                offsetForCurrentPixel +=
                    dimOffset * (sampleStride / baseScales[i]) * multInverse[i];
            }
            offsetForCurrentPixel %= sampleStride;
        }
        pixelForOffset = currentPixel;
    }
	// 在这里，像素是根据 2^j * 3^k 进行跳跃采样的，一遍 2^j * 3^k 是对每个像素做一次采样
	// offsetForCurrentPixel 记录的是，当前像素的偏移
    return offsetForCurrentPixel + sampleNum * sampleStride; 
}

Float HaltonSampler::SampleDimension(int64_t index, int dim) const {
	// sampleAtPixelCenter 是指直接在像素中，取中间做采样
    if (sampleAtPixelCenter && (dim == 0 || dim == 1)) return 0.5f;
    if (dim == 0)
		// 这里要排除 用于进行定位的内容，因为 2^j 或者 3^k 都跟位置有关系
        return RadicalInverse(dim, index >> baseExponents[0]);
    else if (dim == 1)
        return RadicalInverse(dim, index / baseScales[1]);
    else
        return ScrambledRadicalInverse(dim, index,
                                       PermutationForDimension(dim));
}

std::unique_ptr<Sampler> HaltonSampler::Clone(int seed) {
    return std::unique_ptr<Sampler>(new HaltonSampler(*this));
}

HaltonSampler *CreateHaltonSampler(const ParamSet &params,
                                   const Bounds2i &sampleBounds) {
    int nsamp = params.FindOneInt("pixelsamples", 16);
    if (PbrtOptions.quickRender) nsamp = 1;
    bool sampleAtCenter = params.FindOneBool("samplepixelcenter", false);
    return new HaltonSampler(nsamp, sampleBounds, sampleAtCenter);
}

}  // namespace pbrt
