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

template<typename T>
class KNNGNode;

template<typename T>
struct KNNGEdge
{
	T dist;
	int next;
	
	KNNGEdge(T dist, int next)
		: dist(dist), next(next) {}
};

template<typename T>
struct KNNGEdgeComparator
{
	bool operator()(const KNNGEdge<T>& e1, const KNNGEdge<T>& e2)
	{
		return e1.dist < e2.dist;
	}
};


template<typename T>
class KNNGNode
{
public:
	KNNGNode(int K)
		: K_(K) {}

	void InsertNeighbour(T dist, int neighbour)
	{
		neighbours_.emplace(dist, neighbour);
		if (neighbours_.size() > K_)
		{
			neighbours_.pop();
		}
	}
	const std::vector<KNNGEdge<T>>& GetNeighbours() const
	{
		return Container(neighbours_);
	}
private:
	int K_;
	std::priority_queue<KNNGEdge<T>, std::vector<KNNGEdge<T>>, KNNGEdgeComparator<T>> neighbours_;
};

using namespace AppliedGeometry;

class KNNGSolver
{
public:
	template<typename T, int N>
	static std::vector<KNNGNode<T>> SolveNaive(const std::vector<VectorN<N, T>>& points)
	{
		std::vector<KNNGNode<T>> result;
		for (int j = 0; j < points.size(); ++j)
		{
			result.emplace_back();
			KNNGNode<T>& node = result[result.size() - 1];
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
//#undef PARALLEL

	/// Solves the K NNG problem using a parallel algorithm that makes use of Morton ordering
	/// In order to be able to modify the points, they need to be passed in as an rvalue,
	/// otherwise a copy will need to be made.
	template<typename T, int N>
	static std::vector<std::unique_ptr<KNNGNode<T>>> SolveMortonParallel(std::vector<VectorN<N, T>>&& points, int K)
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
							localMinsRes[d] = points[i].GetComponent(d);
						}
					}
				}

				return localMinsRes;
			},
			reduce_struct());

		Morton<N, T> morton(mins.data(), mins.size());

		tbb::parallel_sort(points.begin(), points.end(), morton);

		std::vector<std::unique_ptr<KNNGNode<float>>> edges;
		edges.resize(points.size());
		for (int d = 0; d < edges.size(); ++d)
		{
			edges[d] = std::make_unique<KNNGNode<float>>(K);
		}

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
					float dist = points[i].DistanceSqr(points[i + j]);
					edges[i]->InsertNeighbour(dist, i + j);
				}
			}
		}
#ifdef PARALLEL
		});
#endif

		return edges;
	}

	template<typename T, int N>
	static std::vector<std::unique_ptr<KNNGNode<T>>> SolveMortonParallel(const std::vector<VectorN<N, T>>& points, int K)
	{
		std::vector<VectorN<N, T>> pointsCopy = points;
		return SolveMortonParallel<T, N>(std::move(pointsCopy), K);
	}
};
