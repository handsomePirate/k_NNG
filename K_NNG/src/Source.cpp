#include "KNNGSolver.hpp"
#include "Vector.hpp"
#include <iostream>
#include <fstream>

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
	constexpr int K = 8;
	const int pointCount = 1024;

	std::vector<AppliedGeometry::VectorN<2, type>> points;
	
	const float pmax = 1000.f;
	for (int p = 0; p < pointCount; ++p)
	{
		float x = rand() % int(pmax) / pmax - .5f;
		float y = rand() % int(pmax) / pmax - .5f;
		points.emplace_back(x, y);
	}

	auto result = KNNGSolver::SolveMortonParallel<2, type, K>(std::move(points), false);

	std::ofstream ofsPoints("points.csv");
	ofsPoints << "x,y" << std::endl;
	for (int i = 0; i < points.size(); ++i)
	{
		ofsPoints << points[i].GetComponent(0) << "," << points[i].GetComponent(1) << std::endl;
	}
	ofsPoints.close();

	std::ofstream ofsEdges("edges.csv");
	ofsEdges << "e1,e2" << std::endl;
	for (int i = 0; i < result.size(); ++i)
	{
		auto neighbours = result[i].GetNeighbours();
		for (auto&& node : neighbours)
		{
			ofsEdges << i << "," << node.next << std::endl;
		}
	}
	ofsEdges.close();

	return 0;
	
}