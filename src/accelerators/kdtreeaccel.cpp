
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


// accelerators/kdtreeaccel.cpp*
#include "accelerators/kdtreeaccel.h"
#include "paramset.h"
#include "interaction.h"
#include "stats.h"
#include <algorithm>

namespace pbrt {

// KdTreeAccel Local Declarations
struct KdAccelNode {
    // KdAccelNode Methods
    void InitLeaf(int *primNums, int np, std::vector<int> *primitiveIndices);
    void InitInterior(int axis, int ac, Float s) {
        split = s;
        flags = axis;
        aboveChild |= (ac << 2);
    }
	// 这些函数是用来处理 union 的使用的
    Float SplitPos() const { return split; }
    int nPrimitives() const { return nPrims >> 2; }
    int SplitAxis() const { return flags & 3; }
    bool IsLeaf() const { return (flags & 3) == 3; }
    int AboveChild() const { return aboveChild >> 2; }

	// 重复使用内存块, 让节点大小变得更小, 加上 private 是 8bytes 的节点（一个节点大小
    union {
        Float split;                 // Interior 分轴的位置
        int onePrimitive;            // Leaf 如果只有 1 个 Prim, 记录在 总记录 中的 index, 数据就是 onePrimitive
        int primitiveIndicesOffset;  // Leaf 如果有多个 Prim, 这里记录的在 primitiveIndices 中的偏移, 数据存在 primitiveIndices 中
    };

  private:
    union {
        int flags;       // Both 0:x split-aix, 1:y split-aix, 2:z split-aix, 3:leaf 这里只用了2位 bit
        int nPrims;      // Leaf 叶子节点有几个 Prim, 这里用了 30 bit
        int aboveChild;  // Interior 因为是线性化的存储方式, 所以我们只需要记录一个节点的偏移（另一个节点会在父节点的下一个, 跟 BVH 一致
    };
};

enum class EdgeType { Start, End };
struct BoundEdge {
    // BoundEdge Public Methods
    BoundEdge() {}
    BoundEdge(Float t, int primNum, bool starting) : t(t), primNum(primNum) {
        type = starting ? EdgeType::Start : EdgeType::End;
    }
    Float t;
    int primNum; // 这里会记录它属于哪个 Prim, 这是 Prim 在 AllPrim 中的位置, 因为中途有排序的过程, 所以这里会乱
    EdgeType type;
};

// KdTreeAccel Method Definitions
KdTreeAccel::KdTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
                         int isectCost, int traversalCost, Float emptyBonus,
                         int maxPrims, int maxDepth)
    : isectCost(isectCost),
      traversalCost(traversalCost),
      maxPrims(maxPrims),
      emptyBonus(emptyBonus),
      primitives(std::move(p)) {
    // Build kd-tree for accelerator
    ProfilePhase _(Prof::AccelConstruction);
    nextFreeNode = nAllocedNodes = 0; // 初始化为 0
    if (maxDepth <= 0)
		// 设定最大深度
        maxDepth = std::round(8 + 1.3f * Log2Int(int64_t(primitives.size())));

    // Compute bounds for kd-tree construction
	// 先把 Bounds3f 存一遍（缓存优化
    std::vector<Bounds3f> primBounds;
    primBounds.reserve(primitives.size());
    for (const std::shared_ptr<Primitive> &prim : primitives) {
        Bounds3f b = prim->WorldBound();
        bounds = Union(bounds, b); // 这个是类变量, 扩大一次范围
        primBounds.push_back(b);
    }

    // Allocate working memory for kd-tree construction
	// 边相关的设定, 用于 cost 计算
	// 存储边, 我们每一个 Prim 用 3个轴,  2条边 用于计算
    std::unique_ptr<BoundEdge[]> edges[3];
    for (int i = 0; i < 3; ++i)
        edges[i].reset(new BoundEdge[2 * primitives.size()]);
	// 重复使用的 2个 缓存数组, 用于记录 当前节点，使用到的是 index(在all中的) 列表
    std::unique_ptr<int[]> prims0(new int[primitives.size()]);
    std::unique_ptr<int[]> prims1(new int[(maxDepth + 1) * primitives.size()]);

    // Initialize _primNums_ for kd-tree construction
	// 根节点的表示, 因为它包含了所有的节点, 所以现在 primNums 对应 所有的 primitives, 而且是 一一对应
    std::unique_ptr<int[]> primNums(new int[primitives.size()]);
    for (size_t i = 0; i < primitives.size(); ++i) primNums[i] = i;

    // Start recursive construction of kd-tree
	// 每一个树节点都会调用一次 buildTress
    buildTree(0, bounds, primBounds, primNums.get(), primitives.size(),
              maxDepth, edges, prims0.get(), prims1.get());
}

void KdAccelNode::InitLeaf(int *primNums, int np,
                           std::vector<int> *primitiveIndices) {
    flags = 3; // 这是一个叶子节点
	// 前30位 来表示 Prim 的数量
    nPrims |= (np << 2);
    // Store primitive ids for leaf node
    if (np == 0)
		// 没有 Prim, 直接是 0
        onePrimitive = 0;
    else if (np == 1)
		// 如果只有一个节点, 也不需要添加内存 因为记录总位置就可以了
        onePrimitive = primNums[0];
    else {
		// 记录多位置, 我们把多位置存储到了 primitiveIndices 中, 叶子节点通过->Indices中偏移+np数量->总记录位置->一串Prim 结构很像 VAO VBO 的设定
        primitiveIndicesOffset = primitiveIndices->size();
        for (int i = 0; i < np; ++i) primitiveIndices->push_back(primNums[i]);
    }
}

KdTreeAccel::~KdTreeAccel() { FreeAligned(nodes); }

// nodeNum :  the offset into the array of KdAccelNodes to use for the node that it creates 创建的该节点, 在 KdAccelNodes 中的位移
// nodeBounds ： 整体范围
// allPrimBounds : 所有 Prim 的 Bound，这个是总数据
// primNums : 包含的 Prim 的 在allPrimBounds数组的下标 的数组，当前节点有多少 Prim
// nPrimitives : 包含的 Prim 的个数
// depth： 当前深度
// edges： 边长（3个轴向，前后两条边
void KdTreeAccel::buildTree(int nodeNum, const Bounds3f &nodeBounds,
                            const std::vector<Bounds3f> &allPrimBounds,
                            int *primNums, int nPrimitives, int depth,
                            const std::unique_ptr<BoundEdge[]> edges[3],
                            int *prims0, int *prims1, int badRefines) {
    CHECK_EQ(nodeNum, nextFreeNode);
    // Get next free node from _nodes_ array
    if (nextFreeNode == nAllocedNodes) {
		// 如果节点分配完了, 或者是初始化的时候,  我们分配 2 * nAllocedNodes 个节点
        int nNewAllocNodes = std::max(2 * nAllocedNodes, 512);
		// 开辟内存空间
        KdAccelNode *n = AllocAligned<KdAccelNode>(nNewAllocNodes);
        if (nAllocedNodes > 0) {
			// 复制拷贝
            memcpy(n, nodes, nAllocedNodes * sizeof(KdAccelNode));
			// 释放原来的节点
            FreeAligned(nodes);
        }
        nodes = n;
        nAllocedNodes = nNewAllocNodes;
    }
    ++nextFreeNode; // 表示这个节点被使用了

    // Initialize leaf node if termination criteria met
	// 如果达到了叶子节点的程度, 这里生成一个 叶子节点并返回
    if (nPrimitives <= maxPrims || depth == 0) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Initialize interior node and continue recursion

    // Choose split axis position for interior node
	// 为了得到最小的开销, 我们需要计算出 2 个数, bestAixs 最合适的轴, bestOffset 最合适的偏移
    int bestAxis = -1, bestOffset = -1;
    Float bestCost = Infinity;
    Float oldCost = isectCost * Float(nPrimitives);
    Float totalSA = nodeBounds.SurfaceArea();
    Float invTotalSA = 1 / totalSA; // 总的表面积, 我们先预计算出一个 倒数
    Vector3f d = nodeBounds.pMax - nodeBounds.pMin; // nodeBounds 的 斜向量

    // Choose which axis to split along
    int axis = nodeBounds.MaximumExtent(); // 先尝试最长的那条边, 这个会让子空间更倾向于正方体，而不是长方体，看上去就蛮合理的
    int retries = 0;
retrySplit:

    // Initialize edges for _axis_
    for (int i = 0; i < nPrimitives; ++i) {
        int pn = primNums[i];
        const Bounds3f &bounds = allPrimBounds[pn]; // 找到对应的 Prim 的 Bound
        edges[axis][2 * i] = BoundEdge(bounds.pMin[axis], pn, true);
        edges[axis][2 * i + 1] = BoundEdge(bounds.pMax[axis], pn, false);
    }

    // Sort _edges_ for _axis_
	// 使用 STD 库的 sort
    std::sort(&edges[axis][0], &edges[axis][2 * nPrimitives],
              [](const BoundEdge &e0, const BoundEdge &e1) -> bool {
                  if (e0.t == e1.t)
                      return (int)e0.type < (int)e1.type;
                  else
                      return e0.t < e1.t;
              });

    // Compute cost of all splits for _axis_ to find best
	// nBelow,bAbove 区分 Split 上下的边界
    int nBelow = 0, nAbove = nPrimitives;
    for (int i = 0; i < 2 * nPrimitives; ++i) {
		// 遇到一个 End 的边, 那么 nAbove 就减去 1
        if (edges[axis][i].type == EdgeType::End) --nAbove;
        Float edgeT = edges[axis][i].t;
		// 这里的 = 没有必要算, 因为这样等于没有做切割
        if (edgeT > nodeBounds.pMin[axis] && edgeT < nodeBounds.pMax[axis]) {
            // Compute cost for split at _i_th edge

            // Compute child surface areas for split at _edgeT_
            int otherAxis0 = (axis + 1) % 3, otherAxis1 = (axis + 2) % 3;
			// 这里指的是，假设分成了两块，这两大块的 表面积 计算
            Float belowSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                 (edgeT - nodeBounds.pMin[axis]) *
                                     (d[otherAxis0] + d[otherAxis1]));
            Float aboveSA = 2 * (d[otherAxis0] * d[otherAxis1] +
                                 (nodeBounds.pMax[axis] - edgeT) *
                                     (d[otherAxis0] + d[otherAxis1]));
            Float pBelow = belowSA * invTotalSA;
            Float pAbove = aboveSA * invTotalSA;
            Float eb = (nAbove == 0 || nBelow == 0) ? emptyBonus : 0;
            Float cost =
                traversalCost +
                isectCost * (1 - eb) * (pBelow * nBelow + pAbove * nAbove);
				// 个人对这个 eb, 也就是 emptybounds 的理解, 可能存在一种情况, 有一边分割会很空, 这样分割确实是有利的, 所以有 eb

            // Update best split if this is lowest cost so far
            if (cost < bestCost) {
                bestCost = cost;
                bestAxis = axis;
                bestOffset = i;
            }
        }
        // 如果遇到一个 Start 的边, 那么 nBelow +1
		if (edges[axis][i].type == EdgeType::Start) ++nBelow;
    }
    CHECK(nBelow == nPrimitives && nAbove == 0);

    // Create leaf if no good splits were found
    if (bestAxis == -1 && retries < 2) {
		// 可能出现的情况
		//	做分割时，所有的 Prim 的 Bound 都包含了这个 NodeBound 的范围, 所以根本没有进循环逻辑 -> 单独做一个叶子节点
        ++retries;
        axis = (axis + 1) % 3;
        goto retrySplit;
    }
	// 可能出现的情况
	// 计算出来的 cost 比不拆分 还大，存在情况
	//    1. 事实就是如此, 这个时候 Prim  的个数已经可以接受了, 就  -> 单独做一个叶子节点
	//    2. 可能这一步的算法有问题，如果继续拆分, 会有比较好的结果 -> 默认是尝试 3 次的继续拆分
    if (bestCost > oldCost) ++badRefines;
    if ((bestCost > 4 * oldCost && nPrimitives < 16) || bestAxis == -1 ||
        badRefines == 3) {
        nodes[nodeNum].InitLeaf(primNums, nPrimitives, &primitiveIndices);
        return;
    }

    // Classify primitives with respect to split
	// 注意,  这里的区分是用 Edge 所以会存在有 Prim 在两边都存在的情况
	// 这里很巧妙, 首先是 利用边的 start - end  特性, 把 处于左边的 Prim 放入 prims0 把右边的 Prim 放入 prims1 这里指的是下标记录
    int n0 = 0, n1 = 0;
    for (int i = 0; i < bestOffset; ++i)
        if (edges[bestAxis][i].type == EdgeType::Start)
            prims0[n0++] = edges[bestAxis][i].primNum;
    for (int i = bestOffset + 1; i < 2 * nPrimitives; ++i)
        if (edges[bestAxis][i].type == EdgeType::End)
            prims1[n1++] = edges[bestAxis][i].primNum;

    // Recursively initialize children nodes
    Float tSplit = edges[bestAxis][bestOffset].t;
    Bounds3f bounds0 = nodeBounds, bounds1 = nodeBounds;
    bounds0.pMax[bestAxis] = bounds1.pMin[bestAxis] = tSplit;
	// 因为左边是立刻调用的,  所以我们直接传入 prims0 就算它被覆盖了, 也没有关系, 因为每一次递归会重新生成, 而这一次使用是没有风险的
    buildTree(nodeNum + 1, bounds0, allPrimBounds, prims0, n0, depth - 1, edges,
              prims0, prims1 + nPrimitives, badRefines);
    int aboveChild = nextFreeNode;
    nodes[nodeNum].InitInterior(bestAxis, aboveChild, tSplit);
	// prims1 需要特殊处理, 因为我们传入的时候, 前面会修改 prims1 里面的值, 但是, 对于当前而言 prims1[0~nPrims-1] 这个范围是我们记录了这次递归数据的
	// 所以传入下一次递归的数据 直接冲 prims1 + nP 开始记录即可 (极端情况可能晒满，但一般来说会 空一半
	// 所以这也是 prims1 初始化大小是 (maxDepth + 1) * primitives.size() 的原因
    buildTree(aboveChild, bounds1, allPrimBounds, prims1, n1, depth - 1, edges,
              prims0, prims1 + nPrimitives, badRefines);
}

bool KdTreeAccel::Intersect(const Ray &ray, SurfaceInteraction *isect) const {
    ProfilePhase p(Prof::AccelIntersect);
    // Compute initial parametric range of ray inside kd-tree extent
	// 和最初始的 根节点 做检测
	// tMin 和 tMax 表示当前的 Intersect 范围
    Float tMin, tMax;
    if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }

    // Prepare to traverse kd-tree for ray
	// 预处理
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    PBRT_CONSTEXPR int maxTodo = 64;
	// 存储要做判断的节点
    KdToDo todo[maxTodo]; // TODO 的堆栈
    int todoPos = 0;

    // Traverse kd-tree nodes in order for ray
    bool hit = false;
    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {
        // Bail out if we found a hit closer than the current node
		// 这种情况是可能发生的, 因为有叶子节点中的 Prim，超出了当前节点的方位，导致交点在外面
		// 这种情况需要做一次叶子节点遍历之后, 再做退出
        if (ray.tMax < tMin) break;
        if (!node->IsLeaf()) {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
			// 算出分割平面
            int axis = node->SplitAxis();
            Float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
			// 这里要计算出, 射线是从哪个方向进的, 是先 below 还是先 above, 我们计算相交, 一定要先算前面的
            int belowFirst =
                (ray.o[axis] < node->SplitPos()) ||
                (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tPlane > tMax || tPlane <= 0)
				// 直接跳到下一个子节点, 一个是在内测, 或者是方向相反
                node = firstChild;
            else if (tPlane < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in todo list
				// _secondChild_ 需要入栈
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;
            }
        } else {
            // Check for intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                    primitives[node->onePrimitive];
                // Check one primitive inside leaf node
                if (p->Intersect(ray, isect)) hit = true;
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int index =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &p = primitives[index];
                    // Check one primitive inside leaf node
                    if (p->Intersect(ray, isect)) hit = true;
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
        }
    }
    return hit;
}

bool KdTreeAccel::IntersectP(const Ray &ray) const {
    ProfilePhase p(Prof::AccelIntersectP);
    // Compute initial parametric range of ray inside kd-tree extent
    Float tMin, tMax;
    if (!bounds.IntersectP(ray, &tMin, &tMax)) {
        return false;
    }

    // Prepare to traverse kd-tree for ray
    Vector3f invDir(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
    PBRT_CONSTEXPR int maxTodo = 64;
    KdToDo todo[maxTodo];
    int todoPos = 0;
    const KdAccelNode *node = &nodes[0];
    while (node != nullptr) {
        if (node->IsLeaf()) {
            // Check for shadow ray intersections inside leaf node
            int nPrimitives = node->nPrimitives();
            if (nPrimitives == 1) {
                const std::shared_ptr<Primitive> &p =
                    primitives[node->onePrimitive];
                if (p->IntersectP(ray)) {
                    return true;
                }
            } else {
                for (int i = 0; i < nPrimitives; ++i) {
                    int primitiveIndex =
                        primitiveIndices[node->primitiveIndicesOffset + i];
                    const std::shared_ptr<Primitive> &prim =
                        primitives[primitiveIndex];
                    if (prim->IntersectP(ray)) {
                        return true;
                    }
                }
            }

            // Grab next node to process from todo list
            if (todoPos > 0) {
                --todoPos;
                node = todo[todoPos].node;
                tMin = todo[todoPos].tMin;
                tMax = todo[todoPos].tMax;
            } else
                break;
        } else {
            // Process kd-tree interior node

            // Compute parametric distance along ray to split plane
            int axis = node->SplitAxis();
            Float tPlane = (node->SplitPos() - ray.o[axis]) * invDir[axis];

            // Get node children pointers for ray
            const KdAccelNode *firstChild, *secondChild;
            int belowFirst =
                (ray.o[axis] < node->SplitPos()) ||
                (ray.o[axis] == node->SplitPos() && ray.d[axis] <= 0);
            if (belowFirst) {
                firstChild = node + 1;
                secondChild = &nodes[node->AboveChild()];
            } else {
                firstChild = &nodes[node->AboveChild()];
                secondChild = node + 1;
            }

            // Advance to next child node, possibly enqueue other child
            if (tPlane > tMax || tPlane <= 0)
                node = firstChild;
            else if (tPlane < tMin)
                node = secondChild;
            else {
                // Enqueue _secondChild_ in todo list
                todo[todoPos].node = secondChild;
                todo[todoPos].tMin = tPlane;
                todo[todoPos].tMax = tMax;
                ++todoPos;
                node = firstChild;
                tMax = tPlane;
            }
        }
    }
    return false;
}

std::shared_ptr<KdTreeAccel> CreateKdTreeAccelerator(
    std::vector<std::shared_ptr<Primitive>> prims, const ParamSet &ps) {
    int isectCost = ps.FindOneInt("intersectcost", 80);
    int travCost = ps.FindOneInt("traversalcost", 1); // 这里的比值是 80:1
    Float emptyBonus = ps.FindOneFloat("emptybonus", 0.5f); // 如果切割时，有一边
    int maxPrims = ps.FindOneInt("maxprims", 1);
    int maxDepth = ps.FindOneInt("maxdepth", -1);
    return std::make_shared<KdTreeAccel>(std::move(prims), isectCost, travCost, emptyBonus,
                                         maxPrims, maxDepth);
}

}  // namespace pbrt
