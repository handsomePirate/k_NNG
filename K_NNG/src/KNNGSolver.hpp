#pragma once
#include "Vector.hpp"
#include "PriorityQueueContainer.hpp"
#include "MortonOrdering.hpp"

#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>

#include <chrono>
#include <limits>
#include <vector>

#define MortonMethod 1
#if (MortonMethod == 1)
#define MortonClass Morton
#else
#define MortonClass FloatMortonLess
#endif

#define KNNGTemplate template<typename T, int K>
#define KNNGTemplate2 template<int N, typename T, int K>


KNNGTemplate class KNNGNode;

KNNGTemplate struct KNNGEdge
{
	T dist;
	int next;
	
	KNNGEdge(T dist, int next)
		: dist(dist), next(next) {}
};


KNNGTemplate struct KNNGEdgeComparator
{
	bool operator()(const KNNGEdge<T, K>& e1, const KNNGEdge<T, K>& e2)
	{
		return e1.dist < e2.dist;
	}
};


KNNGTemplate class KNNGNode
{
public:
	KNNGNode() = default;
	void InsertNeighbour(T dist, int neighbour)
	{
		neighbours_.emplace(dist, neighbour);
		if (neighbours_.size() > K)
		{
			neighbours_.pop();
		}
	}
	const std::vector<KNNGEdge<T, K>>& GetNeighbours() const
	{
		return Container(neighbours_);
	}
private:
	std::priority_queue<KNNGEdge<T, K>, std::vector<KNNGEdge<T, K>>, KNNGEdgeComparator<T, K>> neighbours_;
};

using namespace AppliedGeometry;

class KNNGSolver
{
public:
	/*<N, T, K>*/
	KNNGTemplate2 static std::vector<KNNGNode<T, K>> SolveNaive(const std::vector<VectorN<N, T>>& points)
	{
		std::vector<KNNGNode<T, K>> result;
		for (int j = 0; j < points.size(); ++j)
		{
			result.emplace_back();
			KNNGNode<T, K>& node = result[result.size() - 1];
			for (int i = 0; i < points.size(); ++i)
			{
				if (i != j)
				{
					T dist = (points[j] - points[i]).Size();
					node.InsertNeighbour(dist, i);
				}
			}
		}
		return result;
	}

#define PARALLEL
#undef PARALLEL

	/// Solves the K NNG problem using a parallel algorithm that makes use of Morton ordering
	/// In order to be able to modify the points, they need to be passed in as an rvalue,
	/// otherwise a copy will need to be made.
	KNNGTemplate2/*<N, T, K>*/ static std::vector<KNNGNode<T, K>> SolveMortonParallel(
		std::vector<VectorN<N, T>>&& points, bool precise = true)
	{
		std::vector<T> mins;
		mins.resize(N);
		for (T& min : mins)
		{
			min = std::numeric_limits<T>::max();
		}

		struct reduce_struct {
			std::vector<T> operator()(const std::vector<T>& localMins1, const std::vector<T>& localMins2) const
			{
				std::vector<T> localMins = localMins2;
				for (int d = 0; d < localMins1.size(); ++d)
				{
					if (localMins1[d] < localMins[d])
					{
						localMins[d] = localMins1[d];
					}
				}
				return localMins;
			}
		};

		mins = tbb::parallel_reduce(
			tbb::blocked_range<int>(0, points.size()),
			mins,
			[&](tbb::blocked_range<int> r, const std::vector<T>& localMins)
			{
				std::vector<T> localMinsRes = localMins;
				for (int i = r.begin(); i < r.end(); ++i)
				{
					for (int d = 0; d < N; ++d)
					{
						if (points[i].GetComponent(d) < localMinsRes[d])
						{
							localMinsRes[d] = points[i].GetComponent(d) - 0.5f;
						}
					}
				}

				return localMinsRes;
			},
			reduce_struct());

#if (MortonMethod == 1)
		MortonClass<N, T> morton(mins.data(), mins.size());
#else
		MortonClass<N, T> morton;
#endif

		tbb::parallel_sort(points.begin(), points.end(), morton);

		std::vector<KNNGNode<float, K>> edges;
		edges.resize(points.size());

		const float c = 1.f;

		static const int ck = int(c * K);
#ifdef PARALLEL
		tbb::parallel_for(tbb::blocked_range<int>(0, points.size()), [&](tbb::blocked_range<int> r) {
			for (int i = r.begin(); i < r.end(); ++i)
#else
		for (int i = 0; i < points.size(); ++i)
#endif
		{
			for (int j = -ck; j <= ck; ++j)
			{
				if (j != 0 && i + j >= 0 && i + j < points.size())
				{
					T dist = points[i].DistanceSqr(points[i + j]);
					edges[i].InsertNeighbour(dist, i + j);
				}
			}
			if (precise)
			{
				T radius = ceil(FarthestPoint(points[i], edges[i]));
				VectorN<2, T> lower = points[i] - radius;
				VectorN<2, T> upper = points[i] + radius;

				int u;
				int l;
				if (i + ck < points.size() && morton.operator()(upper, points[i + ck]))
				{
					u = i;
				}
				else
				{
					int I = 2;
					for (; i + I < points.size() && morton.operator()(upper, points[i + I]); I <<= 1) {}
					u = i + I;
					if (points.size() <= i + I)
					{
						u = points.size() - 1;
					}
				}
				if (i - ck >= 0 && morton.operator()(points[i - ck], lower))
				{
					l = i;
				}
				else
				{
					int I = 2;
					for (; i - I >= 0 && morton.operator()(points[i - I], lower); I <<= 1) {}
					l = i - I;
					if (0 > i - I)
					{
						l = 0;
					}
				}
				if (u != l)
				{
					//CSearch(edges, points, i, l, u, morton, lower, upper);
				}
			}
		}
#ifdef PARALLEL
		});
#endif

		return edges;
	}

	KNNGTemplate2/*<N, T, K>*/ static std::vector<KNNGNode<T, K>> SolveMortonParallel(
		const std::vector<VectorN<N, T>>& points, bool precise = true)
	{
		std::vector<VectorN<N, T>> pointsCopy = points;
		return SolveMortonParallel<N, T, K>(std::move(pointsCopy), precise);
	}

private:
	template<typename T, int N> static T BoxDist(const VectorN<N, T>& point,
		const VectorN<N, T>& lower, const VectorN<N, T>& upper)
	{
		T dx = (std::max)((std::max)(lower.GetComponent(0) - point.GetComponent(0), T(0)),
			point.GetComponent(0) - upper.GetComponent(0));
		T dy = (std::max)((std::max)(lower.GetComponent(1) - point.GetComponent(1), T(0)),
			point.GetComponent(1) - upper.GetComponent(1));
		return sqrt(dx * dx + dy * dy);
	}

	KNNGTemplate2/*<N, T, K>*/ static T FarthestPoint(const VectorN<N, T>& point,
		const KNNGNode<T, K>& edges)
	{
		const auto& neighbours = edges.GetNeighbours();
		T radius = 0.f;
		for (int j = 0; j < neighbours.size(); ++j)
		{
			if (radius < neighbours[j].dist)
			{
				radius = neighbours[j].dist;
			}
		}
		return radius;
	}

	KNNGTemplate2/*<N, T, K>*/ static void CSearch(std::vector<KNNGNode<T, K>>& edges,
		const std::vector<VectorN<N, T>>& points, int index, int l, int u,
		const MortonClass<N, T>& morton, const VectorN<N, T>& lower, const VectorN<N, T>& upper)
	{
		const int closingConstant = K * 2;
		if (u - l < closingConstant)
		{
			for (int k = l; k <= u; ++k)
			{
				bool present = false;
				for (const KNNGEdge<T, K>& edge : edges[index].GetNeighbours())
				{
					if (edge.next == k)
					{
						present = true;
						break;
					}
				}
				if (k != index && ! present && k >= 0 && k < points.size())
				{
					T dist = points[index].DistanceSqr(points[k]);
					edges[index].InsertNeighbour(dist, k);
				}
			}
			return;
		}

		int m = (u + l) / 2;
		if (m != index)
		{
			T dist = points[index].DistanceSqr(points[m]);
			edges[index].InsertNeighbour(dist, m);
		}

		if (BoxDist(points[index], points[l], points[u]) >= FarthestPoint(points[index], edges[index]))
		{
			return;
		}

		if (morton.operator()(points[index], points[m]))
		{
			CSearch(edges, points, index, l, m - 1, morton, lower, upper);
			if (morton.operator()(points[m], upper))
			{
				CSearch(edges, points, index, m + 1, u, morton, lower, upper);
			}
		}
		else
		{
			CSearch(edges, points, index, m + 1, u, morton, lower, upper);
			if (morton.operator()(lower, points[m]))
			{
				CSearch(edges, points, index, l, m - 1, morton, lower, upper);
			}
		}
	}
};
