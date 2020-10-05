#pragma once
#include "PriorityQueueContainer.hpp"
#include "Vector.hpp"
#include <vector>

#define KNNGTemplate template<typename T, int K, typename std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
#define KNNGTemplate2 template<typename T, int N, int K, typename std::enable_if_t<std::is_floating_point<T>::value, int> = 0>


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

	KNNGTemplate2 static std::vector<KNNGNode<T, K>> SolveDivideAndConquer(const std::vector<VectorN<N, T>>& points)
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

private:

	//static 
};
