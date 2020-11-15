#include "KNNGSolver.hpp"
#include "Vector.hpp"
#include <iostream>
#include <fstream>
#include <string>

std::vector<std::string> Split(const std::string& str, const char delimiter = ' ')
{
	std::vector<std::string> result;
	std::size_t previous = 0;
	std::size_t current = str.find(delimiter);
	while (current != std::string::npos)
	{
		result.push_back(str.substr(previous, current - previous));
		previous = current + 1;
		current = str.find(delimiter, previous);
	}
	result.push_back(str.substr(previous, current - previous));

	return result;
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		std::cout << "Usage: <exe> <txt-with-points>" << std::endl;
		return 0;
	}

	std::string filename = argv[1];
	int dim = 0;
	std::ifstream ifs(filename);
	if (!ifs || ifs.bad())
	{
		std::cout << "Cannot open specified file." << std::endl;
		return 0;
	}
	std::string line;
	std::getline(ifs, line);
	int K;
	try
	{
		K = stoi(line);
	}
	catch (...)
	{
		std::cout << "Wrong file format." << std::endl;
		return 0;
	}
	std::vector<AppliedGeometry::VectorN<2, float>> vec2d;
	std::vector<AppliedGeometry::VectorN<3, float>> vec3d;

	while (ifs.good() && !ifs.eof())
	{
		std::getline(ifs, line);
		auto splitLine = Split(line);
		if (dim == 0)
		{
			dim = splitLine.size();
			if (dim != 2 && dim != 3)
			{
				std::cout << "Wrong file format." << std::endl;
				return 0;
			}
		}
		else
		{
			if (dim != splitLine.size())
			{
				std::cout << "Wrong file format." << std::endl;
				return 0;
			}
		}
		if (dim == 2)
		{
			float x = std::stof(splitLine[0]);
			float y = std::stof(splitLine[1]);
			vec2d.emplace_back(x, y);
		}
		else if (dim == 3)
		{
			float x = std::stof(splitLine[0]);
			float y = std::stof(splitLine[1]);
			float z = std::stof(splitLine[2]);
			vec3d.emplace_back(x, y, z);
		}
	}


	std::vector<std::unique_ptr<KNNGNode<float>>> result;
	if (dim == 2)
	{
		result = KNNGSolver::SolveMortonParallel<float, 2>(std::move(vec2d), K);
	}
	else if (dim == 2)
	{
		result = KNNGSolver::SolveMortonParallel<float, 3>(std::move(vec3d), K);
	}

	std::ofstream ofsPoints("points.csv");
	ofsPoints << "x,y" << std::endl;
	for (int i = 0; i < vec2d.size(); ++i)
	{
		if (dim == 2)
		{
			ofsPoints << vec2d[i].GetComponent(0) << "," << vec2d[i].GetComponent(1) << std::endl;
		}
		else if (dim == 3)
		{
			ofsPoints << vec3d[i].GetComponent(0) << "," << vec3d[i].GetComponent(1) << "," << vec3d[i].GetComponent(2) << std::endl;
		}
	}
	ofsPoints.close();

	std::ofstream ofsEdges("edges.csv");
	ofsEdges << "e1,e2" << std::endl;
	for (int i = 0; i < result.size(); ++i)
	{
		auto neighbours = result[i]->GetNeighbours();
		for (auto&& node : neighbours)
		{
			ofsEdges << i << "," << node.next << std::endl;
		}
	}
	ofsEdges.close();

	return 0;
	
}