#include "KNNGSolver.hpp"
#include "Vector.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>

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
	const int pointCount = 64;

	std::vector<AppliedGeometry::VectorN<2, type>> points;
	
	const float pmax = 1000.f;
	
	//for (int p = 0; p < pointCount; ++p)
	//{
	//	float x = rand() % int(pmax) / pmax - .5f;
	//	float y = rand() % int(pmax) / pmax - .5f;
	//	points.emplace_back(x, y);
	//}
	
	for (int y = 0; y < int(sqrt(pointCount)); ++y)
	{
		for (int x = 0; x < int(sqrt(pointCount)); ++x)
		{
			float fx = x;// (x - .5f) * .25f;
			float fy = y;// (y - .5f) * .25f;
			if (rand() % 2 > 0)
			{
				points.emplace_back(fx, fy);
			}
		}
	}

	//float x = .5f;
	//float y = .5f;
	//points.emplace_back(x, y);

	auto result = KNNGSolver::SolveMortonParallel<2, type, K>(std::move(points), false);
	//auto result2 = KNNGSolver::SolveMortonParallel<2, type, K>(std::move(points), true);
	//auto result3 = KNNGSolver::SolveNaive<2, type, K>(points);

	//int errorCount = 0;
	//for (int i = 0; i < result.size(); ++i)
	//{
	//	for (int k = 0; k < K; ++k)
	//	{
	//		bool found = false;
	//		for (int k2 = 0; k2 < K; ++k2)
	//		{
	//			if (result[i].GetNeighbours()[k].next == result3[i].GetNeighbours()[k2].next)
	//			{
	//				found = true;
	//			}
	//		}
	//		++errorCount;
	//	}
	//}

	//std::cout << errorCount << std::endl;

	std::ofstream ofsPoints("points.csv");
	ofsPoints << "x,y" << std::endl;
	ofsPoints << std::setprecision(16);
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