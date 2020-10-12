#include "KNNGSolver.hpp"
#include "MortonOrdering.hpp"
#include <iostream>
#include <fstream>
#include <chrono>
#include <tbb/parallel_sort.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_for.h>
#include <limits>

#define TEST 1

#if (TEST == 0)
using type = int;
#elif (TEST == 1)
using type = float;
#else
using type = double;
#endif

int main(int argc, char* argv[])
{
	const int K = 1;
	const int pointCount = 1024;

	std::vector<AppliedGeometry::VectorN<2, type>> points;
	
	const float pmax = 1000.f;
	for (int p = 0; p < pointCount; ++p)
	{
		points.emplace_back(rand() % int(pmax) / pmax - .5f, rand() % int(pmax) / pmax - .5f);
	}
	/*
	for (int y = 0; y < sqrt(pointCount); ++y)
	{
		for (int x = 0; x < sqrt(pointCount); ++x)
		{
#if (TEST > 0)
			points.emplace_back(type(x) * type(.5f) - type(sqrt(pointCount)) * type(.25f), 
				type(y) * type(.5f) - type(sqrt(pointCount)) * type(.25f));
#else
			points.emplace_back(type(x), type(y));
#endif
		}
	}
	*/
	std::vector<type> mins;
	mins.resize(points[0].Count());
	for (type& min : mins)
	{
		min = std::numeric_limits<type>::max();
	}

	struct reduce_struct {
		std::vector<type> operator()(const std::vector<type>& localMins1, const std::vector<type>& localMins2) const
		{
			std::vector<type> localMins = localMins2;
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

	auto start = std::chrono::high_resolution_clock::now();

	mins = tbb::parallel_reduce(
		tbb::blocked_range<int>(0, points.size()),
		mins,
		[&](tbb::blocked_range<int> r, const std::vector<type>& localMins)
		{
			std::vector<type> localMinsRes = localMins;
			for (int i = r.begin(); i < r.end(); ++i)
			{
				for (int d = 0; d < points[i].Count(); ++d)
				{
					if (points[i][d] < localMinsRes[d])
					{
						localMinsRes[d] = points[i][d];
					}
				}
			}

			return localMinsRes;
		}, 
		reduce_struct());
	
	Morton<2, type> morton(mins.data(), mins.size());

	auto mid1 = std::chrono::high_resolution_clock::now();

	// TODO: Possibly use the GPU.
	tbb::parallel_sort(points.begin(), points.end(), morton);
	//std::sort(points.begin(), points.end(), morton);
	
	auto mid2 = std::chrono::high_resolution_clock::now();
	
	//std::cout << "x,y" << std::endl;
	//for (const auto& point : points)
	//{
	//	std::cout << point.GetComponent(0) << "," << point.GetComponent(1) << std::endl;
	//}
	//return 0;

	std::vector<KNNGNode<float, K>> edges;
	edges.resize(points.size());

	const float c = 1.f;
	/*
	static const int ck = int(c * K);
	for (int i = 0; i < points.size(); ++i)
	{
		for (int j = -ck; j <= ck; ++j)
		{
			if (j != 0 && i + j >= 0 && i + j < points.size())
			{
				float dist = points[i].DistanceSqr(points[i + j]);
				edges[i].InsertNeighbour(dist, i + j);
			}
		}
	}
	*/
	tbb::parallel_for(tbb::blocked_range<int>(0, points.size()), [&](tbb::blocked_range<int> r)
		{
			static const int ck = int(c * K);
			for (int i = r.begin(); i < r.end(); ++i)
			{
				for (int j = -ck; j <= ck; ++j)
				{
					if (j != 0 && i + j >= 0 && i + j < points.size())
					{
						float dist = points[i].DistanceSqr(points[i + j]);
						edges[i].InsertNeighbour(dist, i + j);
					}
				}
			}
		});
	
	auto end = std::chrono::high_resolution_clock::now();

	//std::cout << std::endl;
	//std::cout << std::chrono::duration_cast<std::chrono::microseconds>(mid1 - start).count() / 1000.f << " ms." << std::endl;
	//std::cout << std::chrono::duration_cast<std::chrono::microseconds>(mid2 - mid1).count() / 1000.f << " ms." << std::endl;
	//std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - mid2).count() / 1000.f << " ms." << std::endl;

	//return 0;

	//auto start = std::chrono::high_resolution_clock::now();

	//auto result = KNNGSolver::SolveNaive<type, 2, K>(points);

	//auto end = std::chrono::high_resolution_clock::now();

	std::ofstream ofsPoints("points.csv");
	ofsPoints << "x,y" << std::endl;
	for (int i = 0; i < points.size(); ++i)
	{
		ofsPoints << points[i].GetComponent(0) << "," << points[i].GetComponent(1) << std::endl;
	}
	ofsPoints.close();

	std::ofstream ofsEdges("edges.csv");
	ofsEdges << "e1,e2" << std::endl;
	for (int i = 0; i < edges.size(); ++i)
	{
		auto neighbours = edges[i].GetNeighbours();
		for (auto&& node : neighbours)
		{
			ofsEdges << i << "," << node.next << std::endl;
		}
	}
	ofsEdges.close();
	//for (const auto& node : edges)
	//{
	//	auto neighbours = node.GetNeighbours();
	//	std::sort(neighbours.begin(), neighbours.end(), KNNGEdgeComparator<float, K>());
	//	const int neighbourCount = neighbours.size() > K ? K : neighbours.size();
	//	for (int k = 0; k < neighbourCount; ++k)
	//	{
	//		std::cout << neighbours[k].next << '(' << neighbours[k].dist << ')';
	//		if (k < K - 1)
	//		{
	//			std::cout << ' ';
	//		}
	//	}
	//	std::cout << std::endl;
	//}
	

	//std::cout << std::endl;
	//std::cout << "K = " << K << " nearest neighbours in ";
	//std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.f << " ms." << std::endl;

	return 0;
	
}