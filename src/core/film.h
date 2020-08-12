
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

#ifndef PBRT_CORE_FILM_H
#define PBRT_CORE_FILM_H

// core/film.h*
#include "pbrt.h"
#include "geometry.h"
#include "spectrum.h"
#include "filter.h"
#include "stats.h"
#include "parallel.h"

namespace pbrt {

// FilmTilePixel Declarations
// FilmTile 的每一个像素，需要维护两个数，
struct FilmTilePixel {
	// the weighted contributions from the pixel samples (according to the reconstruction filter weights)
	// pixel 样本的权重和
    Spectrum contribSum = 0.f;
	// the filter weights
	// filter 滤波器的权重和
    Float filterWeightSum = 0.f;
};

// Film Declarations
class Film {
  public:
    // Film Public Methods
	// 构造函数传入的参数:
	// resolution：图像的整体分辨率
	// cropWindow：裁剪窗口，用于指定要渲染的图像子集
	// filter：滤波器指针
	// diagonal：Film 对应的物理区域的 对角线长度，传入的参数是毫米，但是在构造时会转换成米
	// filename：输出图片文件名
	// scale，maxSampleLuminance：图像的 Pixel Value 怎么在 files 中保存
    Film(const Point2i &resolution, const Bounds2f &cropWindow,
         std::unique_ptr<Filter> filter, Float diagonal,
         const std::string &filename, Float scale,
         Float maxSampleLuminance = Infinity);
    Bounds2i GetSampleBounds() const;
    Bounds2f GetPhysicalExtent() const;
    std::unique_ptr<FilmTile> GetFilmTile(const Bounds2i &sampleBounds);
    void MergeFilmTile(std::unique_ptr<FilmTile> tile);
    void SetImage(const Spectrum *img) const;
    void AddSplat(const Point2f &p, Spectrum v);
    void WriteImage(Float splatScale = 1);
    void Clear();

    // Film Public Data
    const Point2i fullResolution;
    const Float diagonal;
    std::unique_ptr<Filter> filter;
    const std::string filename;
	// 裁剪的像素范围
    Bounds2i croppedPixelBounds;

  private:
    // Film Private Data

	// 我们为在分辨率上的图像的每一个像素，分配了一个 Pixel 的结构
    struct Pixel {
        Pixel() { xyz[0] = xyz[1] = xyz[2] = filterWeightSum = 0; }
		// Spectrum 转换成的 XYZ 表示
        Float xyz[3];
		// 滤波器的权重值的和
        Float filterWeightSum;
		// splats 样本值
        AtomicFloat splatXYZ[3];
		// pad 不使用，为了结构上的优化，4 byte，凑够 32 bytes
        Float pad;
    };
    std::unique_ptr<Pixel[]> pixels;
    static PBRT_CONSTEXPR int filterTableWidth = 16;
	// 预计算的滤波器权重表
    Float filterTable[filterTableWidth * filterTableWidth];
    std::mutex mutex;
    const Float scale;
    const Float maxSampleLuminance;

    // Film Private Methods
    Pixel &GetPixel(const Point2i &p) {
        CHECK(InsideExclusive(p, croppedPixelBounds));
        int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
        int offset = (p.x - croppedPixelBounds.pMin.x) +
                     (p.y - croppedPixelBounds.pMin.y) * width;
        return pixels[offset];
    }
};

class FilmTile {
  public:
    // FilmTile Public Methods
	// 构造一个2D边界框，提供了最终的像素边界，以及做一些预计算
	// filterTable 一个 float 指针，里面是 filter 预计算的值
    FilmTile(const Bounds2i &pixelBounds, const Vector2f &filterRadius,
             const Float *filterTable, int filterTableSize,
             Float maxSampleLuminance)
        : pixelBounds(pixelBounds),
          filterRadius(filterRadius),
          invFilterRadius(1 / filterRadius.x, 1 / filterRadius.y),
          filterTable(filterTable),
          filterTableSize(filterTableSize),
          maxSampleLuminance(maxSampleLuminance) {
		// 初始化 pixelBounds 大小的 Pixel，每个 Piexl 里面有一个 FilmTilePixel 数据
        pixels = std::vector<FilmTilePixel>(std::max(0, pixelBounds.Area()));
    }
	// 当一条射线的一个样本对应辐射度被计算时，积分器能会调用 AddSample
	// 它的参数是
	// 获取一个样本 pFilm位置 其辐射值 L 以及 Camera::GenerateRayDifferential 计算出的权重
	// 这个函数的功能是，将该样本，对 Pixel 的共享结果，也就是滤波方程的计算结果，存在起来
    void AddSample(const Point2f &pFilm, Spectrum L,
                   Float sampleWeight = 1.) {
        ProfilePhase _(Prof::AddFilmSample);
        if (L.y() > maxSampleLuminance)
            L *= maxSampleLuminance / L.y();
        // Compute sample's raster bounds
		// 计算 该样本影响的像素范围
		// 连续转离散 -0.5 处理
        Point2f pFilmDiscrete = pFilm - Vector2f(0.5f, 0.5f);
		// 滤波器范围
        Point2i p0 = (Point2i)Ceil(pFilmDiscrete - filterRadius);
        Point2i p1 =
            (Point2i)Floor(pFilmDiscrete + filterRadius) + Point2i(1, 1);
		// 是否出界，出界就不用管了
        p0 = Max(p0, pixelBounds.pMin);
        p1 = Min(p1, pixelBounds.pMax);

        // Loop over filter support and add sample to pixel arrays

        // Precompute $x$ and $y$ filter table offsets
		// 滤波器权重表已经算好了，都在 filterTable 指针中
		// 这里做的是，我们做一个位置的映射表，这里节省的是
		// 我们做 2次 p1.x + p1.y 的循环的计算即可
		// 而不是在下面做 p1.x * p1.y 的循环次数 的 位置映射 的 计算
        int *ifx = ALLOCA(int, p1.x - p0.x);
        for (int x = p0.x; x < p1.x; ++x) {
            Float fx = std::abs((x - pFilmDiscrete.x) * invFilterRadius.x *
                                filterTableSize);
            ifx[x - p0.x] = std::min((int)std::floor(fx), filterTableSize - 1);
        }
        int *ify = ALLOCA(int, p1.y - p0.y);
        for (int y = p0.y; y < p1.y; ++y) {
            Float fy = std::abs((y - pFilmDiscrete.y) * invFilterRadius.y *
                                filterTableSize);
            ify[y - p0.y] = std::min((int)std::floor(fy), filterTableSize - 1);
        }
		// 计算在该样本影响下，该样本，对每个像素位置的贡献度
        for (int y = p0.y; y < p1.y; ++y) {
            for (int x = p0.x; x < p1.x; ++x) {
                // Evaluate filter value at $(x,y)$ pixel
                int offset = ify[y - p0.y] * filterTableSize + ifx[x - p0.x];
                Float filterWeight = filterTable[offset];

                // Update pixel values with filtered sample contribution
				// 计算流程，就是累加贡献度 的内容
                FilmTilePixel &pixel = GetPixel(Point2i(x, y));
                pixel.contribSum += L * sampleWeight * filterWeight;
                pixel.filterWeightSum += filterWeight;
            }
        }
    }
	// 获取某个像素坐标下，我们累计的 Piexl 值
    FilmTilePixel &GetPixel(const Point2i &p) {
        CHECK(InsideExclusive(p, pixelBounds));
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset =
            (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }
	// const 版本
    const FilmTilePixel &GetPixel(const Point2i &p) const {
        CHECK(InsideExclusive(p, pixelBounds));
        int width = pixelBounds.pMax.x - pixelBounds.pMin.x;
        int offset =
            (p.x - pixelBounds.pMin.x) + (p.y - pixelBounds.pMin.y) * width;
        return pixels[offset];
    }
	// 返回的是，它的 影响的像素 范围
    Bounds2i GetPixelBounds() const { return pixelBounds; }

  private:
    // FilmTile Private Data
    const Bounds2i pixelBounds;
    const Vector2f filterRadius, invFilterRadius;
    const Float *filterTable;
    const int filterTableSize;
    std::vector<FilmTilePixel> pixels;
    const Float maxSampleLuminance;
    friend class Film;
};

Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter);

}  // namespace pbrt

#endif  // PBRT_CORE_FILM_H
