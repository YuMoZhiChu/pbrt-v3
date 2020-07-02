
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


// samplers/stratified.cpp*
#include "samplers/stratified.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// StratifiedSampler Method Definitions
void StratifiedSampler::StartPixel(const Point2i &p) {
    ProfilePhase _(Prof::StartPixel);
    // Generate single stratified samples for the pixel
    for (size_t i = 0; i < samples1D.size(); ++i) {
		// 对samples1D[i]，第i个维度做随机赋值
        StratifiedSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng,
                           jitterSamples);
		// 因为上面得到的数据 samples1D 会跟着 i 增大, 趋势是变大的，所以要洗一次牌
		// 因此，这里得到的 samp 数据，其实是有顺序的，为了打乱他们的顺序，要洗一次牌
        Shuffle(&samples1D[i][0], xPixelSamples * yPixelSamples, 1, rng);
    }
    for (size_t i = 0; i < samples2D.size(); ++i) {
        StratifiedSample2D(&samples2D[i][0], xPixelSamples, yPixelSamples, rng,
                           jitterSamples);
        Shuffle(&samples2D[i][0], xPixelSamples * yPixelSamples, 1, rng);
    }

    // Generate arrays of stratified samples for the pixel
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i)
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            int count = samples1DArraySizes[i];
			// 注意，这里和 1D，2D 不同，这里的分层是 count，我们是对 一个数据里面的数组，做 count 层，的分层采样
            StratifiedSample1D(&sampleArray1D[i][j * count], count, rng,
                               jitterSamples);
            Shuffle(&sampleArray1D[i][j * count], count, 1, rng);
        }
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i)
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            int count = samples2DArraySizes[i];
			// 所以这里就会有问题，因为 count 是任意的，但是我们要对 一个 2D 数组做分层采样，所以引入了 拉丁超Cube算法

			// 传 x 是因为 PointX 是 x在前面, &point.x 等于 &point
            LatinHypercube(&sampleArray2D[i][j * count].x, count, 2, rng);
        }
    PixelSampler::StartPixel(p);
}

std::unique_ptr<Sampler> StratifiedSampler::Clone(int seed) {
    StratifiedSampler *ss = new StratifiedSampler(*this);
    ss->rng.SetSequence(seed);
    return std::unique_ptr<Sampler>(ss);
}

StratifiedSampler *CreateStratifiedSampler(const ParamSet &params) {
    bool jitter = params.FindOneBool("jitter", true);
    int xsamp = params.FindOneInt("xsamples", 4);
    int ysamp = params.FindOneInt("ysamples", 4);
    int sd = params.FindOneInt("dimensions", 4);
    if (PbrtOptions.quickRender) xsamp = ysamp = 1;
    return new StratifiedSampler(xsamp, ysamp, jitter, sd);
}

}  // namespace pbrt
