#include "KNNGSolver.hpp"
#include "MortonOrdering.hpp"
#include <iostream>
#include <chrono>

#define TEST 2

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
	
	//const float pmax = 1000.f;
	//for (int p = 0; p < pointCount; ++p)
	//{
	//	points.emplace_back(rand() % int(pmax) / pmax - .5f, rand() % int(pmax) / pmax - .5f);
	//}

	for (int y = 0; y < sqrt(pointCount); ++y)
	{
		for (int x = 0; x < sqrt(pointCount); ++x)
		{
			points.emplace_back(type(x) * type(.125f)/* - .2f*/, type(y) * type(.125f)/* - .2f*/);
			//points.emplace_back(type(x), type(y));
		}
	}
	
	std::vector<type> mins;
	mins.resize(points[0].Count());
	for (int i = 0; i < points.size(); ++i)
	{
		for (int d = 0; d < points[i].Count(); ++d)
		{
			mins[d] = points[i][d] > mins[d] ? mins[d] : points[i][d];
		}
	}
	
	Morton<2, type> morton(mins.data(), mins.size());
	std::sort(points.begin(), points.end(), morton);

	std::cout << "x,y" << std::endl;
	for (const auto& point : points)
	{
		std::cout << point.GetComponent(0) << "," << point.GetComponent(1) << std::endl;
	}
	return 0;

	auto start = std::chrono::high_resolution_clock::now();

	//auto result = KNNGSolver::SolveNaive<float, 2, K>(points);

	auto end = std::chrono::high_resolution_clock::now();

	/*
	for (const auto& node : result)
	{
		auto neighbours = node.GetNeighbours();
		std::sort(neighbours.begin(), neighbours.end(), KNNGEdgeComparator<float, K>());
		const int neighbourCount = neighbours.size() > K ? K : neighbours.size();
		for (int k = 0; k < neighbourCount; ++k)
		{
			std::cout << neighbours[k].next << '(' << neighbours[k].dist << ')';
			if (k < K - 1)
			{
				std::cout << ' ';
			}
		}
		std::cout << std::endl;
	}
	*/

	std::cout << std::endl;
	std::cout << "K = " << K << " nearest neighbours in ";
	std::cout << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.f << " ms." << std::endl;

	return 0;
}