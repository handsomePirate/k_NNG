#include "KNNGSolver.hpp"
#include "MortonOrdering.hpp"
#include <iostream>
#include <chrono>

int main(int argc, char* argv[])
{
	const int K = 1;
	const int pointCount = 8;

	std::vector<AppliedGeometry::VectorN<2, float>> points;
	std::vector<AppliedGeometry::VectorN<2, int>> pointsint;

#define POINTS points
	
	//for (int p = 0; p < pointCount; ++p)
	//{
	//	points.emplace_back(rand() % 100 / 100.f, rand() % 100 / 100.f);
	//	pointsint.emplace_back(rand() % 100, rand() % 100);
	//}

	for (int y = 0; y < 8; ++y)
	{
		for (int x = 0; x < 8; ++x)
		{
			POINTS.emplace_back(float(x)/* * 0.1f*/, float(y)/* * 0.1f*/);
		}
	}

	std::sort(POINTS.begin(), POINTS.end(), Morton());

	for (const auto& point : POINTS)
	{
		std::cout << "(" << point.GetComponent(0) << ", " << point.GetComponent(1) << ")" << std::endl;
	}
	return 0;

	auto start = std::chrono::high_resolution_clock::now();

	auto result = KNNGSolver::SolveNaive<float, 2, K>(points);

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