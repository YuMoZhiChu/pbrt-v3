
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


// core/film.cpp*
#include "film.h"
#include "paramset.h"
#include "imageio.h"
#include "stats.h"

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/Film pixels", filmPixelMemory);

// Film Method Definitions
Film::Film(const Point2i &resolution, const Bounds2f &cropWindow,
           std::unique_ptr<Filter> filt, Float diagonal,
           const std::string &filename, Float scale, Float maxSampleLuminance)
    : fullResolution(resolution),
      diagonal(diagonal * .001),
      filter(std::move(filt)),
      filename(filename),
      scale(scale),
      maxSampleLuminance(maxSampleLuminance) {
    // Compute film image bounds
	// 初始化裁剪窗口大小，NDC空间->具体的像素位置
	// 使用 ceil 函数，保证如果有多个渲染流程，每一个像素出现在子图像中，只会出现一次（约定都使用 ceil
    croppedPixelBounds =
        Bounds2i(Point2i(std::ceil(fullResolution.x * cropWindow.pMin.x),
                         std::ceil(fullResolution.y * cropWindow.pMin.y)),
                 Point2i(std::ceil(fullResolution.x * cropWindow.pMax.x),
                         std::ceil(fullResolution.y * cropWindow.pMax.y)));
    LOG(INFO) << "Created film with full resolution " << resolution <<
        ". Crop window of " << cropWindow << " -> croppedPixelBounds " <<
        croppedPixelBounds;

    // Allocate film image storage
	// 为每一个像素，分配动态内存
	// 这是我们要显示的内容的 Pixel
    pixels = std::unique_ptr<Pixel[]>(new Pixel[croppedPixelBounds.Area()]);
    filmPixelMemory += croppedPixelBounds.Area() * sizeof(Pixel);

    // Precompute filter weight table
	// 预计算滤波器的表，不用每一次都调用 filter->Evaluate 来计算
	// 可能会导致 filter->Evaluate 传入的位置不够精确，但这个错误可以忽略
    int offset = 0;
    for (int y = 0; y < filterTableWidth; ++y) {
        for (int x = 0; x < filterTableWidth; ++x, ++offset) {
            Point2f p;
			// 从中心点开始，扩展，并进行预计算
            p.x = (x + 0.5f) * filter->radius.x / filterTableWidth;
            p.y = (y + 0.5f) * filter->radius.y / filterTableWidth;
            filterTable[offset] = filter->Evaluate(p);
        }
    }
}

// 这个函数负责返回给出 采样器的范围，因为我们的实际范围，要拓展到滤波器的过滤范围
// 我们在采样的时候，经常跨越多个像素，所以采样范围，是大于实际的输出范围的
// 所以图像不仅仅收到内部的值的影响，也会收到外部的影响
Bounds2i Film::GetSampleBounds() const {
	// 考虑到离散像素坐标 转换为 连续像素坐标时，需要做 半像素偏移，再拓展滤波器半径，向外舍入
    Bounds2f floatBounds(Floor(Point2f(croppedPixelBounds.pMin) +
                               Vector2f(0.5f, 0.5f) - filter->radius),
                         Ceil(Point2f(croppedPixelBounds.pMax) -
                              Vector2f(0.5f, 0.5f) + filter->radius));
    return (Bounds2i)floatBounds;
}

// 这个函数，返回在 film 在场景中的实际范围，可以理解为这一块 幕布 在实际场景中的大小
// RealisticCamera 会需要这个参数
// 给定的 a 是长宽比
// 给定的 diagonal 是斜边长
Bounds2f Film::GetPhysicalExtent() const {
    Float aspect = (Float)fullResolution.y / (Float)fullResolution.x;
    Float x = std::sqrt(diagonal * diagonal / (1 + aspect * aspect));
    Float y = aspect * x;
    return Bounds2f(Point2f(-x / 2, -y / 2), Point2f(x / 2, y / 2));
}

// 为了适应多线程的优化，这里提供了一个 Tile 的定义
// 参数是传入一个 采样范围（样本范围），返回一个 unique 的 FileTile 指针
// 每一块 FilmTile 指向了一部分的 像素信息，这一部分的信息是 这个对应范围内，对像素的贡献程度（也就是累加结果），而线程对这个 Tile 是 unique 掌控的
// 所以不需要担心和其他线程的竞争关系
// 完成了 Tile 上的工作后，我们会传回 Tile，并且做一次合并
std::unique_ptr<FilmTile> Film::GetFilmTile(const Bounds2i &sampleBounds) {
    // Bound image pixels that samples in _sampleBounds_ contribute to
    Vector2f halfPixel = Vector2f(0.5f, 0.5f);
    Bounds2f floatBounds = (Bounds2f)sampleBounds;
	// 给定了样本的范围，现在我们要计算，样本影响到的像素的对应范围，很显然，就是以样本为中心，扩大滤波器的范围
	// 考量到 离散坐标转连续坐标的影响，这里参考 Sampling_Theory.html 部分，将离散的点，转换为连续的点的做法，为什么大的那部分也是减呢，因为我们做的是一个范围，而不是当个的点
	// 加减滤波器的半径
    Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
    Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
                 Point2i(1, 1);
	// 要和图形边界求交，图像外的像素无需考虑
	// 这个范围，就是，给定的 样本范围，能够影响的 像素范围
    Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
	// 返回一个 FilmTile 指针
    return std::unique_ptr<FilmTile>(new FilmTile(
        tilePixelBounds, filter->radius, filterTable, filterTableWidth,
        maxSampleLuminance));
}

void Film::Clear() {
    for (Point2i p : croppedPixelBounds) {
        Pixel &pixel = GetPixel(p);
        for (int c = 0; c < 3; ++c)
            pixel.splatXYZ[c] = pixel.xyz[c] = 0;
        pixel.filterWeightSum = 0;
    }
}

// 渲染线程调用，合并 FilmTile，把计算结果叠加
// 这里传入一个 FilmTile 的 u_q，所以在传入的时候，FilmTile 的拥有者，就是当前函数了
// 完成了本次调用后，内存会被自动释放，就不可以再使用 FilmTile 并往里面写数据了
void Film::MergeFilmTile(std::unique_ptr<FilmTile> tile) {
    ProfilePhase p(Prof::MergeFilmTile);
    VLOG(1) << "Merging film tile " << tile->pixelBounds;
	// 用一个互斥锁，保证每次只有一次做数据更新
    std::lock_guard<std::mutex> lock(mutex);
	// 只需要遍历，这一块，影响的 像素范围即可
    for (Point2i pixel : tile->GetPixelBounds()) {
        // Merge _pixel_ into _Film::pixels_
        const FilmTilePixel &tilePixel = tile->GetPixel(pixel);
        Pixel &mergePixel = GetPixel(pixel);
        Float xyz[3];
        tilePixel.contribSum.ToXYZ(xyz);
		// 转成 xyz
		// 合并到，我们 Film 类的 Pixel 中，这里存的 Pixel 是 XYZ
        for (int i = 0; i < 3; ++i) mergePixel.xyz[i] += xyz[i];
		// 滤波器权重和相加，这个也是分母
        mergePixel.filterWeightSum += tilePixel.filterWeightSum;
    }
}

// 这个是给整个 Film，一次性提供所有的光谱值
void Film::SetImage(const Spectrum *img) const {
    int nPixels = croppedPixelBounds.Area();
    for (int i = 0; i < nPixels; ++i) {
        Pixel &p = pixels[i];
        img[i].ToXYZ(p.xyz); // 转成 xyz 直接传入
        p.filterWeightSum = 1;
        p.splatXYZ[0] = p.splatXYZ[1] = p.splatXYZ[2] = 0;
    }
}

// 一些光线传入算法，比如 16.3 会出现的  bidirectional path tracing
// 需要一个 "Splat" 贡献度，可能是对任意一个像素的
// 在这里，Splat 不是做加权求和，而是做简单的相加
// 这里大致的想法是，如果一个像素，它旁边的 Splat 越多，他会越亮
// Pixel::splatXYZ 是原子类型的 AtomicFloat，允许多个写出直接通过 AddSplat 相加，不用写锁
void Film::AddSplat(const Point2f &p, Spectrum v) {
    ProfilePhase pp(Prof::SplatFilm);

    if (v.HasNaNs()) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with NaN values "
                                   "at (%f, %f)", p.x, p.y);
        return;
    } else if (v.y() < 0.) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with negative "
                                   "luminance %f at (%f, %f)", v.y(), p.x, p.y);
        return;
    } else if (std::isinf(v.y())) {
        LOG(ERROR) << StringPrintf("Ignoring splatted spectrum with infinite "
                                   "luminance at (%f, %f)", p.x, p.y);
        return;
    }

    Point2i pi = Point2i(Floor(p));
    if (!InsideExclusive(pi, croppedPixelBounds)) return;
    if (v.y() > maxSampleLuminance)
        v *= maxSampleLuminance / v.y();
    Float xyz[3];
    v.ToXYZ(xyz);
    Pixel &pixel = GetPixel(pi);
	// 直接相加即可，不做加权求和
    for (int i = 0; i < 3; ++i) pixel.splatXYZ[i].Add(xyz[i]);
}

// 将图像写入文件中
// 传入一个 splat 缩放因子，这里会在 AddSplact 中用到，详细看 16.4.5 的内容
void Film::WriteImage(Float splatScale) {
    // Convert image to RGB and compute final pixel values
    LOG(INFO) <<
        "Converting image to RGB and computing final weighted pixel values";
	// 首先是分配对应大小，存储 rgb 的图像内容
    std::unique_ptr<Float[]> rgb(new Float[3 * croppedPixelBounds.Area()]);
    int offset = 0;
	// 存储流程，遍历所有的像素点
    for (Point2i p : croppedPixelBounds) {
        // Convert pixel XYZ color to RGB
		// 转换到 RGB，存储到 RGB 空间
        Pixel &pixel = GetPixel(p);
        XYZToRGB(pixel.xyz, &rgb[3 * offset]);

        // Normalize pixel with weight sum
        Float filterWeightSum = pixel.filterWeightSum;
		// 除法计算
        if (filterWeightSum != 0) {
			// 因为可能存在负数，所以要规范到 > 0
            Float invWt = (Float)1 / filterWeightSum;
            rgb[3 * offset] = std::max((Float)0, rgb[3 * offset] * invWt);
            rgb[3 * offset + 1] =
                std::max((Float)0, rgb[3 * offset + 1] * invWt);
            rgb[3 * offset + 2] =
                std::max((Float)0, rgb[3 * offset + 2] * invWt);
        }

        // Add splat value at pixel
        Float splatRGB[3];
        Float splatXYZ[3] = {pixel.splatXYZ[0], pixel.splatXYZ[1],
                             pixel.splatXYZ[2]};
		// 同时计算 Splat，并做缩放计算
        XYZToRGB(splatXYZ, splatRGB);
        rgb[3 * offset] += splatScale * splatRGB[0];
        rgb[3 * offset + 1] += splatScale * splatRGB[1];
        rgb[3 * offset + 2] += splatScale * splatRGB[2];

        // Scale pixel value by _scale_
		// 最终写入的像素，用户指定的系数，一般来说，都会是1????
        rgb[3 * offset] *= scale;
        rgb[3 * offset + 1] *= scale;
        rgb[3 * offset + 2] *= scale;
        ++offset;
    }

    // Write RGB image
    LOG(INFO) << "Writing image " << filename << " with bounds " <<
        croppedPixelBounds;
	// 调用 pbrt 的接口，如果写入的是8位整数格式，在其转换为整数前，会进行 gamma 校正 ????
    pbrt::WriteImage(filename, &rgb[0], croppedPixelBounds, fullResolution);
}

Film *CreateFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
    std::string filename;
    if (PbrtOptions.imageFile != "") {
        filename = PbrtOptions.imageFile;
        std::string paramsFilename = params.FindOneString("filename", "");
        if (paramsFilename != "")
            Warning(
                "Output filename supplied on command line, \"%s\" is overriding "
                "filename provided in scene description file, \"%s\".",
                PbrtOptions.imageFile.c_str(), paramsFilename.c_str());
    } else
        filename = params.FindOneString("filename", "pbrt.exr");

    int xres = params.FindOneInt("xresolution", 1280);
    int yres = params.FindOneInt("yresolution", 720);
    if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
    Bounds2f crop;
    int cwi;
    const Float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
        crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
        crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
        crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
    } else if (cr)
        Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);
    else
        crop = Bounds2f(Point2f(Clamp(PbrtOptions.cropWindow[0][0], 0, 1),
                                Clamp(PbrtOptions.cropWindow[1][0], 0, 1)),
                        Point2f(Clamp(PbrtOptions.cropWindow[0][1], 0, 1),
                                Clamp(PbrtOptions.cropWindow[1][1], 0, 1)));

    Float scale = params.FindOneFloat("scale", 1.);
    Float diagonal = params.FindOneFloat("diagonal", 35.);
    Float maxSampleLuminance = params.FindOneFloat("maxsampleluminance",
                                                   Infinity);
    return new Film(Point2i(xres, yres), crop, std::move(filter), diagonal,
                    filename, scale, maxSampleLuminance);
}

}  // namespace pbrt
