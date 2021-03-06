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


// accelerators/bvh.cpp*
#include "accelerators/bvh.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "parallel.h"
#include <algorithm>

namespace pbrt {

STAT_MEMORY_COUNTER("Memory/BVH tree", treeBytes);
STAT_RATIO("BVH/Primitives per leaf node", totalPrimitives, totalLeafNodes);
STAT_COUNTER("BVH/Interior nodes", interiorNodes);
STAT_COUNTER("BVH/Leaf nodes", leafNodes);

// BVHAccel Local Declarations
struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() {}
    BVHPrimitiveInfo(size_t primitiveNumber, const Bounds3f &bounds)
        : primitiveNumber(primitiveNumber),
          bounds(bounds),
		  // 初始化时，计算重心
          centroid(.5f * bounds.pMin + .5f * bounds.pMax) {}
    size_t primitiveNumber;
    Bounds3f bounds;
    Point3f centroid;
};

// BVH 树节点
struct BVHBuildNode {
    // BVHBuildNode Public Methods
    void InitLeaf(int first, int n, const Bounds3f &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
		// 叶子节点 没有子节点
        children[0] = children[1] = nullptr;
        ++leafNodes;
        ++totalLeafNodes;
        totalPrimitives += n;
    }
	// 当子节点存在时, 父节点非常好实现, 挂两个指针, Union 一次 Bound 即可
    void InitInterior(int axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
        ++interiorNodes;
    }
	// 表示所有子节点的边界 Union做法得到的 bounds
    Bounds3f bounds;
	// 非叶子节点, 记录他们的左右孩子
    BVHBuildNode *children[2];
	// splitAixs 记录将 2个子节点 分开的轴, 提高遍历的效率[在 InsertP 遍历中, 该轴对应的 Dir方向, 决定了先看左节点, 还是先看右节点, 因为我们求的是 hit min
	// firstPrimOffset, nPrimitives 用于描述 叶子节点 包含的 primitive, [firstPrimOffset, firstPrimOffset+nPrimitives) 不包含右边
	// 因此，需要对 prim 数组进行重新排序
    int splitAxis, firstPrimOffset, nPrimitives;
};

struct MortonPrimitive {
	// 它的在 primitiveInfo  下标
    int primitiveIndex;
	// 莫顿码
    uint32_t mortonCode;
};

// 表示一颗子树
struct LBVHTreelet {
	// 起点, 个数
    int startIndex, nPrimitives;
	// 结果节点, 这是指向子树根节点的指针
    BVHBuildNode *buildNodes;
};

// 线性节点的结构, 将2叉树做列表进行保存
struct LinearBVHNode {
    Bounds3f bounds;
    union {
        int primitivesOffset;   // leaf  对于叶子节点, 这个相当于 BVHBuildNode 的 firstPrimOffset
        int secondChildOffset;  // interior 记录中间节点的右节点, 它的偏移量
    };
    uint16_t nPrimitives;  // 0 -> interior node 中间节点是 0, 对于叶子节点, 这个相当于 BVHBuildNode 的 nPrimitives
    uint8_t axis;          // interior node: xyz
    uint8_t pad[1];        // ensure 32 byte total size 为了支持CPU的高速缓存架构, 这里做内存对齐
};

// BVHAccel Utility Functions
inline uint32_t LeftShift3(uint32_t x) {
    CHECK_LE(x, (1 << 10));
    if (x == (1 << 10)) --x; // 如果恰好等于边界, 就减去1
#ifdef PBRT_HAVE_BINARY_CONSTANTS
	// 0b 做为 C++ 的2进制拓展
	// 这一段的目的是, 对 x 做一个维度的 Morton 编码, 目的是减少计算次数
	// http://www.pbr-book.org/3ed-2018/Primitives_and_Intersection_Acceleration/Bounding_Volume_Hierarchies.html#LeftShift3
    x = (x | (x << 16)) & 0b00000011000000000000000011111111;
    // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0b00000011000000001111000000001111;
    // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0b00000011000011000011000011000011;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0b00001001001001001001001001001001;
    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#else
    x = (x | (x << 16)) & 0x30000ff;
    // x = ---- --98 ---- ---- ---- ---- 7654 3210
    x = (x | (x << 8)) & 0x300f00f;
    // x = ---- --98 ---- ---- 7654 ---- ---- 3210
    x = (x | (x << 4)) & 0x30c30c3;
    // x = ---- --98 ---- 76-- --54 ---- 32-- --10
    x = (x | (x << 2)) & 0x9249249;
    // x = ---- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
#endif // PBRT_HAVE_BINARY_CONSTANTS
    return x;
}

inline uint32_t EncodeMorton3(const Vector3f &v) {
    CHECK_GE(v.x, 0);
    CHECK_GE(v.y, 0);
    CHECK_GE(v.z, 0);
	// 结果上 z 左移2次, y 左移1次
    return (LeftShift3(v.z) << 2) | (LeftShift3(v.y) << 1) | LeftShift3(v.x);
}

// 基数排序的写法, 这里是 lowbit 开始的, 写得非常简练
// 对于这里的二进制, 基数排序 会比 快排 和 插入排序做混合的 std::sort 要快
static void RadixSort(std::vector<MortonPrimitive> *v) {
    std::vector<MortonPrimitive> tempVector(v->size());
    PBRT_CONSTEXPR int bitsPerPass = 6; // 我们选取的排序位数是 6 也就是 2^6
    PBRT_CONSTEXPR int nBits = 30; // 总共 30 位数需要排序, 所以我们需要执行 5 次
    static_assert((nBits % bitsPerPass) == 0,
                  "Radix sort bitsPerPass must evenly divide nBits");
    PBRT_CONSTEXPR int nPasses = nBits / bitsPerPass;

    for (int pass = 0; pass < nPasses; ++pass) {
        // Perform one pass of radix sort, sorting _bitsPerPass_ bits
		// 从 lowBit 开始运算
        int lowBit = pass * bitsPerPass;

        // Set in and out vector pointers for radix sort pass
        std::vector<MortonPrimitive> &in = (pass & 1) ? tempVector : *v;
        std::vector<MortonPrimitive> &out = (pass & 1) ? *v : tempVector;

        // Count number of zero bits in array for current radix sort bit
        PBRT_CONSTEXPR int nBuckets = 1 << bitsPerPass;
        int bucketCount[nBuckets] = {0};
        PBRT_CONSTEXPR int bitMask = (1 << bitsPerPass) - 1;
        for (const MortonPrimitive &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
            CHECK_GE(bucket, 0);
            CHECK_LT(bucket, nBuckets);
            ++bucketCount[bucket]; // 记录该基数出现的次数
        }

        // Compute starting index in output array for each bucket
        int outIndex[nBuckets];
        outIndex[0] = 0;
        for (int i = 1; i < nBuckets; ++i)
			// 记录前一项 i-1 的 buckets 总共有多少个, 所以这个也是它的序号
            outIndex[i] = outIndex[i - 1] + bucketCount[i - 1];

        // Store sorted values in output array
        for (const MortonPrimitive &mp : in) {
            int bucket = (mp.mortonCode >> lowBit) & bitMask;
			// 序号+1 记录内容
            out[outIndex[bucket]++] = mp;
        }
    }
    // Copy final result from _tempVector_, if needed
	// 我们要返回的是 v
    if (nPasses & 1) std::swap(*v, tempVector);
}

// BVHAccel Method Definitions
BVHAccel::BVHAccel(std::vector<std::shared_ptr<Primitive>> p,
                   int maxPrimsInNode, SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)),
      splitMethod(splitMethod),
      primitives(std::move(p)) {
    ProfilePhase _(Prof::AccelConstruction);
    if (primitives.empty()) return;
    // Build BVH from _primitives_

    // Initialize _primitiveInfo_ array for primitives
	// 存入 prim 的 index 和 bound
    std::vector<BVHPrimitiveInfo> primitiveInfo(primitives.size());
    for (size_t i = 0; i < primitives.size(); ++i)
		// 这里记录了 i 是在 info 在 vector 中的位置
        primitiveInfo[i] = {i, primitives[i]->WorldBound()};

    // Build BVH tree for primitives using _primitiveInfo_
    MemoryArena arena(1024 * 1024);
    int totalNodes = 0;
    std::vector<std::shared_ptr<Primitive>> orderedPrims;
	// 预先分配内存
    orderedPrims.reserve(primitives.size());
    BVHBuildNode *root;
    if (splitMethod == SplitMethod::HLBVH)
        root = HLBVHBuild(arena, primitiveInfo, &totalNodes, orderedPrims);
    else
		// 返回 BVH 树的根节点, 节点被放在 arena 中的内存中, 总节点数 totalNodes
		// orderedPrims 会带上顺序
        root = recursiveBuild(arena, primitiveInfo, 0, primitives.size(),
                              &totalNodes, orderedPrims);
	
	primitives.swap(orderedPrims);
    primitiveInfo.resize(0);
    LOG(INFO) << StringPrintf("BVH created with %d nodes for %d "
                              "primitives (%.2f MB), arena allocated %.2f MB",
                              totalNodes, (int)primitives.size(),
                              float(totalNodes * sizeof(LinearBVHNode)) /
                              (1024.f * 1024.f),
                              float(arena.TotalAllocated()) /
                              (1024.f * 1024.f));

    // Compute representation of depth-first traversal of BVH tree
    treeBytes += totalNodes * sizeof(LinearBVHNode) + sizeof(*this) +
                 primitives.size() * sizeof(primitives[0]);
    nodes = AllocAligned<LinearBVHNode>(totalNodes); // 这里就是 totalNode 的用处 AllocAligned 是 _aligned_malloc 的高层封装
    int offset = 0;
	// 二叉树压缩成线性
    flattenBVHTree(root, &offset);
    CHECK_EQ(totalNodes, offset);
}

Bounds3f BVHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : Bounds3f();
}

struct BucketInfo {
    int count = 0;
    Bounds3f bounds;
};

// 递归式调用
// start, end 代表 [start, end) 即 primInfo[start] 到 primInfo[end-1] 的  prim 子集
// [start, end) 只有一个值时, 返回一个 叶子节点
// 所以该方法就是 [start, mid) [mid, end) 依次递归
// totalNodes 用于构造更紧凑的 LinerBVHTree
BVHBuildNode *BVHAccel::recursiveBuild(
    MemoryArena &arena, std::vector<BVHPrimitiveInfo> &primitiveInfo, int start,
    int end, int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims) {
    CHECK_NE(start, end);
    BVHBuildNode *node = arena.Alloc<BVHBuildNode>();
    (*totalNodes)++;
    // Compute bounds of all primitives in BVH node
    Bounds3f bounds;
    for (int i = start; i < end; ++i)
        bounds = Union(bounds, primitiveInfo[i].bounds);
    int nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _BVHBuildNode_
		// 每生成一个 叶子节点
		// 我们记录位置，并且把 prim 压入 orderedPrims
		// 这么做的好处: 1. 得到有序的 primsInfo 2. BVH树中只需要记录 offset
        int firstPrimOffset = orderedPrims.size();
        for (int i = start; i < end; ++i) {
            int primNum = primitiveInfo[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;
    } else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        Bounds3f centroidBounds;
        for (int i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, primitiveInfo[i].centroid);
		// dim 是拓展范围最宽的一条轴, 在这里算了之后, 这个子集里面, 就以 dim 进行分类
		// 遍历一次取最大
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        int mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
			// 如果大家都重复在一个地方, 这里明显是重合了, 那么也视为叶子节点
            // Create leaf _BVHBuildNode_
            int firstPrimOffset = orderedPrims.size();
            for (int i = start; i < end; ++i) {
                int primNum = primitiveInfo[i].primitiveNumber;
                orderedPrims.push_back(primitives[primNum]);
            }
            node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
            return node;
        } else {
            // Partition primitives based on _splitMethod_
			// 基于不同的算法做分类
            switch (splitMethod) {
            case SplitMethod::Middle: {
                // Partition primitives through node's midpoint
				// C++ https://en.cppreference.com/w/cpp/algorithm/partition
                Float pmid =
                    (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]) / 2;
                BVHPrimitiveInfo *midPtr = std::partition(
                    &primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                    [dim, pmid](const BVHPrimitiveInfo &pi) {
                        return pi.centroid[dim] < pmid;
                    });
				// C++ 的例子中, 这里是用指针减法的，但是因为这里涉及到 start 这个递归变量
				// 所以要玩指针运算
                mid = midPtr - &primitiveInfo[0];
                // For lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall
                // through
                // to EqualCounts.
				// 这种情况是可能出现的, 当一个超级大的 prim 包裹住所有的 prim 时, 这里就无法做区分了
				// 如果没有 break 掉, 那么会执行下一个 case
                if (mid != start && mid != end) break;
            }
            case SplitMethod::EqualCounts: {
                // Partition primitives into equally-sized subsets
				// 退化成取中间
                mid = (start + end) / 2;
				// nth_element 的时间复杂度是 On, NlogN 次交换，比全排序要好上一点
				// 这里的结果, [start, mid) [mid, end) 这两段里面，不保证有序
				// https://zh.cppreference.com/w/cpp/algorithm/nth_element
                std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
                                 &primitiveInfo[end - 1] + 1,
                                 [dim](const BVHPrimitiveInfo &a,
                                       const BVHPrimitiveInfo &b) {
                                     return a.centroid[dim] < b.centroid[dim];
                                 });
                break;
            }
            case SplitMethod::SAH:
            default: {
                // Partition primitives using approximate SAH
                if (nPrimitives <= 2) {
					// 如果只有两个的情况下, 这里做 SAH 贪心算法的运算是相当不划算的（教程上是4个
					// 直接对半分即可
                    // Partition primitives into equally-sized subsets
                    mid = (start + end) / 2;
                    std::nth_element(&primitiveInfo[start], &primitiveInfo[mid],
                                     &primitiveInfo[end - 1] + 1,
                                     [dim](const BVHPrimitiveInfo &a,
                                           const BVHPrimitiveInfo &b) {
                                         return a.centroid[dim] <
                                                b.centroid[dim];
                                     });
                } else {
                    // Allocate _BucketInfo_ for SAH partition buckets
					// 固定分成 12 个 Bucket
                    PBRT_CONSTEXPR int nBuckets = 12;
                    BucketInfo buckets[nBuckets];

                    // Initialize _BucketInfo_ for SAH partition buckets
                    for (int i = start; i < end; ++i) {
						// 算出这个 primitiveInfo[i] 在第几个 bucket
                        int b = nBuckets *
                                centroidBounds.Offset(
                                    primitiveInfo[i].centroid)[dim];
                        if (b == nBuckets) b = nBuckets - 1;
                        CHECK_GE(b, 0);
                        CHECK_LT(b, nBuckets);
						// 增加数量
                        buckets[b].count++;
						// 拓展 bucket 的 bound，合并在它区间内的 prim
                        buckets[b].bounds =
                            Union(buckets[b].bounds, primitiveInfo[i].bounds);
                    }

                    // Compute costs for splitting after each bucket
					// 因为是 n 个 bucket，所以有 n-1 种 split 方法，我们把 任何一次 RayInsert 计算的代价 看做 1
                    Float cost[nBuckets - 1];
					// 这里的开销是 On^2 但是 n 只有 11 完全可以接受
                    for (int i = 0; i < nBuckets - 1; ++i) {
                        Bounds3f b0, b1;
                        int count0 = 0, count1 = 0;
                        for (int j = 0; j <= i; ++j) {
                            b0 = Union(b0, buckets[j].bounds);
                            count0 += buckets[j].count;
                        }
                        for (int j = i + 1; j < nBuckets; ++j) {
                            b1 = Union(b1, buckets[j].bounds);
                            count1 += buckets[j].count;
                        }
						// 这个 1 是 Ttrav 的代价, 有数学意义, 但是对比较来说没太大用处, 用 0.125,0 都可以
						// 这里用的是 表面积 来表示概率 ???? 原因1: 射线只和表面相交
                        cost[i] = 1 +
                                  (count0 * b0.SurfaceArea() +
                                   count1 * b1.SurfaceArea()) /
                                      bounds.SurfaceArea();
                    }

                    // Find bucket to split at that minimizes SAH metric
                    Float minCost = cost[0];
                    int minCostSplitBucket = 0;
                    for (int i = 1; i < nBuckets - 1; ++i) {
                        if (cost[i] < minCost) {
                            minCost = cost[i];
                            minCostSplitBucket = i;
                        }
                    }

                    // Either create leaf or split primitives at selected SAH
                    // bucket
                    Float leafCost = nPrimitives;
					// 起码得比遍历开销低
                    if (nPrimitives > maxPrimsInNode || minCost < leafCost) {
                        BVHPrimitiveInfo *pmid = std::partition(
                            &primitiveInfo[start], &primitiveInfo[end - 1] + 1,
                            [=](const BVHPrimitiveInfo &pi) {
								// 算出 pi 属于哪个 bucket
                                int b = nBuckets *
                                        centroidBounds.Offset(pi.centroid)[dim];
                                if (b == nBuckets) b = nBuckets - 1;
                                CHECK_GE(b, 0);
                                CHECK_LT(b, nBuckets);
                                return b <= minCostSplitBucket;
                            });
                        mid = pmid - &primitiveInfo[0];
                    } else {
                        // Create leaf _BVHBuildNode_
                        int firstPrimOffset = orderedPrims.size();
                        for (int i = start; i < end; ++i) {
                            int primNum = primitiveInfo[i].primitiveNumber;
                            orderedPrims.push_back(primitives[primNum]);
                        }
                        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
                        return node;
                    }
                }
                break;
            }
            }
			// 递归执行
            node->InitInterior(dim,
                               recursiveBuild(arena, primitiveInfo, start, mid,
                                              totalNodes, orderedPrims),
                               recursiveBuild(arena, primitiveInfo, mid, end,
                                              totalNodes, orderedPrims));
        }
    }
    return node;
}

// 使用 Morton codes
BVHBuildNode *BVHAccel::HLBVHBuild(
    MemoryArena &arena, const std::vector<BVHPrimitiveInfo> &primitiveInfo,
    int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims) const {
    // Compute bounding box of all primitive centroids
	// LBVH 只考虑 Prim 的重心, 不会去关注他们实际的展开大小
    Bounds3f bounds;
    for (const BVHPrimitiveInfo &pi : primitiveInfo)
        bounds = Union(bounds, pi.centroid);

    // Compute Morton indices of primitives
	// 得到总体的范围 Overrall Bounds, 计算每个 Prim 的 Morton Code
    std::vector<MortonPrimitive> mortonPrims(primitiveInfo.size());
	// lambda 用法 & 表示 所有变量可以用引用的方式获取
    ParallelFor([&](int i) {
		// 每组512个 prim
        // Initialize _mortonPrims[i]_ for _i_th primitive
        PBRT_CONSTEXPR int mortonBits = 10;
        PBRT_CONSTEXPR int mortonScale = 1 << mortonBits;
        mortonPrims[i].primitiveIndex = primitiveInfo[i].primitiveNumber;
        Vector3f centroidOffset = bounds.Offset(primitiveInfo[i].centroid); // 这里得到的是 [0,1] 的三维数组
        mortonPrims[i].mortonCode = EncodeMorton3(centroidOffset * mortonScale);
    }, primitiveInfo.size(), 512);

    // Radix sort primitive Morton indices
	// 对 mortonPrim 排序
    RadixSort(&mortonPrims);
	// 这里得到的结果, motronPrims 中的 Prim 已经是根据 他们的 MortonCode 排好序的了

    // Create LBVH treelets at bottom of BVH

    // Find intervals of primitives for each treelet
    std::vector<LBVHTreelet> treeletsToBuild;
    for (int start = 0, end = 1; end <= (int)mortonPrims.size(); ++end) {
#ifdef PBRT_HAVE_BINARY_CONSTANTS
      uint32_t mask = 0b00111111111111000000000000000000;
	  // 进行一次莫顿的分组, 前12个数字, 表示把整个空间分为了 2^12=4096 组, 每个维度是 2^4=16 组
#else
      uint32_t mask = 0x3ffc0000;
#endif
      if (end == (int)mortonPrims.size() ||
            ((mortonPrims[start].mortonCode & mask) !=
             (mortonPrims[end].mortonCode & mask))) {
            // Add entry to _treeletsToBuild_ for this treelet
            int nPrimitives = end - start;
            int maxBVHNodes = 2 * nPrimitives; // BVH树的节点数量, 是叶子节点的两倍 - 1
			// false 表示不需要执行 BVHBuildNode 的初始构造函数, 因为在接下来的并发中, Node会执行一次构造, 这里做并发会消耗大量的性能
            BVHBuildNode *nodes = arena.Alloc<BVHBuildNode>(maxBVHNodes, false); // 在线性内容中分配内存会更好 ????
            treeletsToBuild.push_back({start, nPrimitives, nodes});
            start = end;
        }
    }

    // Create LBVHs for treelets in parallel
	// 并发的构建子树
	// 记录总数, 偏移, 这些需要线程安全的操作
    std::atomic<int> atomicTotal(0), orderedPrimsOffset(0);
    orderedPrims.resize(primitives.size());
    ParallelFor([&](int i) {
        // Generate _i_th LBVH treelet
        int nodesCreated = 0;
        const int firstBitIndex = 29 - 12; // 因为是30位的莫顿码, 分了 12 组, 所以他们的 Index 是 29-12
        LBVHTreelet &tr = treeletsToBuild[i];
        tr.buildNodes =
            emitLBVH(tr.buildNodes, primitiveInfo, &mortonPrims[tr.startIndex],
                     tr.nPrimitives, &nodesCreated, orderedPrims,
                     &orderedPrimsOffset, firstBitIndex);
        atomicTotal += nodesCreated;
    }, treeletsToBuild.size());
    *totalNodes = atomicTotal;

    // Create and return SAH BVH from LBVH treelets
    std::vector<BVHBuildNode *> finishedTreelets;
    finishedTreelets.reserve(treeletsToBuild.size());
    for (LBVHTreelet &treelet : treeletsToBuild)
        finishedTreelets.push_back(treelet.buildNodes);
	// 因为这里只有 4096 棵树, 所以建立 SAH 会很快
    return buildUpperSAH(arena, finishedTreelets, 0, finishedTreelets.size(),
                         totalNodes);
}

// 子树也要继续做分割, 基于莫顿码的特性
BVHBuildNode *BVHAccel::emitLBVH(
    BVHBuildNode *&buildNodes,
    const std::vector<BVHPrimitiveInfo> &primitiveInfo,
    MortonPrimitive *mortonPrims, int nPrimitives, int *totalNodes,
    std::vector<std::shared_ptr<Primitive>> &orderedPrims,
    std::atomic<int> *orderedPrimsOffset, int bitIndex) const {
    CHECK_GT(nPrimitives, 0);
    if (bitIndex == -1 || nPrimitives < maxPrimsInNode) {
        // Create and return leaf node of LBVH treelet
		// 当分到最细时, 或者当前节点数, 少于一个节点可以容纳的最大节点时(一般是4),创建一个叶子节点
        (*totalNodes)++;
        BVHBuildNode *node = buildNodes++; // buildNodes 递进, 因为新增了一个叶子节点
        Bounds3f bounds;
		// 这个变量需要原子操作来执行
        int firstPrimOffset = orderedPrimsOffset->fetch_add(nPrimitives);
        for (int i = 0; i < nPrimitives; ++i) {
            int primitiveIndex = mortonPrims[i].primitiveIndex;
            orderedPrims[firstPrimOffset + i] = primitives[primitiveIndex]; // 这样做是安全的, 因为只会修改到 firstPrimOffset - i
            bounds = Union(bounds, primitiveInfo[primitiveIndex].bounds);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bounds);
        return node;
    } else {
        int mask = 1 << bitIndex;
        // Advance to next subtree level if there's no LBVH split for this bit
        if ((mortonPrims[0].mortonCode & mask) ==
            (mortonPrims[nPrimitives - 1].mortonCode & mask))
			// 如果不可拆分的话, 那就走到下一级
            return emitLBVH(buildNodes, primitiveInfo, mortonPrims, nPrimitives,
                            totalNodes, orderedPrims, orderedPrimsOffset,
                            bitIndex - 1);

        // Find LBVH split point for this dimension
		// 2分查找, 找到 mask & 不一致的位置
        int searchStart = 0, searchEnd = nPrimitives - 1;
        while (searchStart + 1 != searchEnd) {
            CHECK_NE(searchStart, searchEnd);
            int mid = (searchStart + searchEnd) / 2;
            if ((mortonPrims[searchStart].mortonCode & mask) ==
                (mortonPrims[mid].mortonCode & mask))
                searchStart = mid;
            else {
                CHECK_EQ(mortonPrims[mid].mortonCode & mask,
                         mortonPrims[searchEnd].mortonCode & mask);
                searchEnd = mid;
            }
        }
        int splitOffset = searchEnd;
        CHECK_LE(splitOffset, nPrimitives - 1);
        CHECK_NE(mortonPrims[splitOffset - 1].mortonCode & mask,
                 mortonPrims[splitOffset].mortonCode & mask);

        // Create and return interior LBVH node
        (*totalNodes)++;
        BVHBuildNode *node = buildNodes++; // buildNodes 递进, 因为新增了一个中间节点
        BVHBuildNode *lbvh[2] = {
            emitLBVH(buildNodes, primitiveInfo, mortonPrims, splitOffset,
                     totalNodes, orderedPrims, orderedPrimsOffset,
                     bitIndex - 1),
            emitLBVH(buildNodes, primitiveInfo, &mortonPrims[splitOffset],
                     nPrimitives - splitOffset, totalNodes, orderedPrims,
                     orderedPrimsOffset, bitIndex - 1)};
        int axis = bitIndex % 3; // 区分的轴
        node->InitInterior(axis, lbvh[0], lbvh[1]); // 这里也无需排序,  因为 lbvh[0] 的莫顿码 一定比 lbvh[1] 高
        return node;
    }
}

// 对 LBVH 做 SAH 算法
BVHBuildNode *BVHAccel::buildUpperSAH(MemoryArena &arena,
                                      std::vector<BVHBuildNode *> &treeletRoots,
                                      int start, int end,
                                      int *totalNodes) const {
    CHECK_LT(start, end);
    int nNodes = end - start;
    if (nNodes == 1) return treeletRoots[start];
    (*totalNodes)++;
    BVHBuildNode *node = arena.Alloc<BVHBuildNode>();

    // Compute bounds of all nodes under this HLBVH node
    Bounds3f bounds;
    for (int i = start; i < end; ++i)
        bounds = Union(bounds, treeletRoots[i]->bounds);

    // Compute bound of HLBVH node centroids, choose split dimension _dim_
    Bounds3f centroidBounds;
    for (int i = start; i < end; ++i) {
        Point3f centroid =
            (treeletRoots[i]->bounds.pMin + treeletRoots[i]->bounds.pMax) *
            0.5f;
        centroidBounds = Union(centroidBounds, centroid);
    }
    int dim = centroidBounds.MaximumExtent();
    // FIXME: if this hits, what do we need to do?
    // Make sure the SAH split below does something... ?
    CHECK_NE(centroidBounds.pMax[dim], centroidBounds.pMin[dim]);

    // Allocate _BucketInfo_ for SAH partition buckets
    PBRT_CONSTEXPR int nBuckets = 12;
    struct BucketInfo {
        int count = 0;
        Bounds3f bounds;
    };
    BucketInfo buckets[nBuckets];

    // Initialize _BucketInfo_ for HLBVH SAH partition buckets
    for (int i = start; i < end; ++i) {
        Float centroid = (treeletRoots[i]->bounds.pMin[dim] +
                          treeletRoots[i]->bounds.pMax[dim]) *
                         0.5f;
        int b =
            nBuckets * ((centroid - centroidBounds.pMin[dim]) /
                        (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
        if (b == nBuckets) b = nBuckets - 1;
        CHECK_GE(b, 0);
        CHECK_LT(b, nBuckets);
        buckets[b].count++;
        buckets[b].bounds = Union(buckets[b].bounds, treeletRoots[i]->bounds);
    }

    // Compute costs for splitting after each bucket
    Float cost[nBuckets - 1];
    for (int i = 0; i < nBuckets - 1; ++i) {
        Bounds3f b0, b1;
        int count0 = 0, count1 = 0;
        for (int j = 0; j <= i; ++j) {
            b0 = Union(b0, buckets[j].bounds);
            count0 += buckets[j].count;
        }
        for (int j = i + 1; j < nBuckets; ++j) {
            b1 = Union(b1, buckets[j].bounds);
            count1 += buckets[j].count;
        }
        cost[i] = .125f +
                  (count0 * b0.SurfaceArea() + count1 * b1.SurfaceArea()) /
                      bounds.SurfaceArea();
    }

    // Find bucket to split at that minimizes SAH metric
    Float minCost = cost[0];
    int minCostSplitBucket = 0;
    for (int i = 1; i < nBuckets - 1; ++i) {
        if (cost[i] < minCost) {
            minCost = cost[i];
            minCostSplitBucket = i;
        }
    }

    // Split nodes and create interior HLBVH SAH node
    BVHBuildNode **pmid = std::partition(
        &treeletRoots[start], &treeletRoots[end - 1] + 1,
        [=](const BVHBuildNode *node) {
            Float centroid =
                (node->bounds.pMin[dim] + node->bounds.pMax[dim]) * 0.5f;
            int b = nBuckets *
                    ((centroid - centroidBounds.pMin[dim]) /
                     (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
            if (b == nBuckets) b = nBuckets - 1;
            CHECK_GE(b, 0);
            CHECK_LT(b, nBuckets);
            return b <= minCostSplitBucket;
        });
    int mid = pmid - &treeletRoots[0];
    CHECK_GT(mid, start);
    CHECK_LT(mid, end);
    node->InitInterior(
        dim, this->buildUpperSAH(arena, treeletRoots, start, mid, totalNodes),
        this->buildUpperSAH(arena, treeletRoots, mid, end, totalNodes));
    return node;
}

int BVHAccel::flattenBVHTree(BVHBuildNode *node, int *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    int myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        CHECK(!node->children[0] && !node->children[1]);
        CHECK_LT(node->nPrimitives, 65536);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        // Create interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
		// 直接递归生成子树即可
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset =
            flattenBVHTree(node->children[1], offset);
    }
	// 返回这棵树 的 总体偏移  myOffset - origin_offset = 树的长度
    return myOffset;
}

BVHAccel::~BVHAccel() { FreeAligned(nodes);/*释放 nodes 的内存*/ }

bool BVHAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    if (!nodes) return false;
    ProfilePhase p(Prof::AccelIntersect);
    bool hit = false;
	// 预计算一些内容
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    // Follow ray through BVH nodes to find primitive intersections
	// toVisitOffset 栈记录 currentNodeIndex 在node列表中的偏移
    int toVisitOffset = 0, currentNodeIndex = 0;
    int nodesToVisit[64]; // 深度遍历的栈
    while (true) {
		// 开始遍历
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        // Check ray against BVH node
        if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
				// 这里是叶子节点
                // Intersect ray with primitives in leaf BVH node
                for (int i = 0; i < node->nPrimitives; ++i)
					// 如果找到了交点, 很可能不止一个, 所以需要继续遍历
					// 而且我们的 Intersect 会更新 tMax, 可以避免许多的判断
                    if (primitives[node->primitivesOffset + i]->Intersect(
                            ray, isect))
                        hit = true;
                if (toVisitOffset == 0) break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
				// 这里是中间节点
                // Put far BVH node on _nodesToVisit_ stack, advance to near
                // node
				// 这里考虑到有2个子节点, 为了不让子节点进行后续的比较, 这里要从前往后
				// 对于正向的 d, 我们先检测左节点, 后检测右节点, 反之亦然
                if (dirIsNeg[node->axis]) {
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            if (toVisitOffset == 0) break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }
    return hit;
}

bool BVHAccel::IntersectP(const Ray &ray) const {
    if (!nodes) return false;
    ProfilePhase p(Prof::AccelIntersectP);
    Vector3f invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    int dirIsNeg[3] = {invDir.x < 0, invDir.y < 0, invDir.z < 0};
    int nodesToVisit[64];
    int toVisitOffset = 0, currentNodeIndex = 0;
    while (true) {
        const LinearBVHNode *node = &nodes[currentNodeIndex];
        if (node->bounds.IntersectP(ray, invDir, dirIsNeg)) {
            // Process BVH node _node_ for traversal
            if (node->nPrimitives > 0) {
                for (int i = 0; i < node->nPrimitives; ++i) {
                    if (primitives[node->primitivesOffset + i]->IntersectP(
                            ray)) {
                        return true;
                    }
                }
                if (toVisitOffset == 0) break;
                currentNodeIndex = nodesToVisit[--toVisitOffset];
            } else {
                if (dirIsNeg[node->axis]) {
                    /// second child first
                    nodesToVisit[toVisitOffset++] = currentNodeIndex + 1;
                    currentNodeIndex = node->secondChildOffset;
                } else {
                    nodesToVisit[toVisitOffset++] = node->secondChildOffset;
                    currentNodeIndex = currentNodeIndex + 1;
                }
            }
        } else {
            if (toVisitOffset == 0) break;
            currentNodeIndex = nodesToVisit[--toVisitOffset];
        }
    }
    return false;
}

std::shared_ptr<BVHAccel> CreateBVHAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    std::string splitMethodName = ps.FindOneString("splitmethod", "sah");
    BVHAccel::SplitMethod splitMethod;
    if (splitMethodName == "sah")
        splitMethod = BVHAccel::SplitMethod::SAH;
    else if (splitMethodName == "hlbvh")
        splitMethod = BVHAccel::SplitMethod::HLBVH;
    else if (splitMethodName == "middle")
        splitMethod = BVHAccel::SplitMethod::Middle;
    else if (splitMethodName == "equal")
        splitMethod = BVHAccel::SplitMethod::EqualCounts;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                splitMethodName.c_str());
        splitMethod = BVHAccel::SplitMethod::SAH;
    }

    int maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return std::make_shared<BVHAccel>(std::move(prims), maxPrimsInNode, splitMethod);
}

}  // namespace pbrt
