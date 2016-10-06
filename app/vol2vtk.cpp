/*
 * vol2vtk.cpp
 *
 *  Created on: May 26, 2015
 *      Author: cotrik
 */

#include "include.h"

using namespace std;

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "Usage: vol2vtk vol_file tet_file" << std::endl;
		return -1;
	}

	ifstream file(argv[1]);

	std::vector <Vertex> V;
	std::vector <Cell> C;

	std::string str(1024, 0);
	// read cells
	unsigned long cellsNum = 0;
	while (file.getline((char *)str.data(), 1024))
	{
		if (str.find("volumeelements") != str.npos)
		{
			file.getline((char *)str.data(), 1024);
			stringstream strstream(str.c_str());
			strstream >> cellsNum;
			std::cout << "cells : " << cellsNum << std::endl;
			for (unsigned long i = 0; i < cellsNum; i++)
			{
				file.getline((char *)str.data(), 1024);
				stringstream cell_strstream(str.c_str());
				unsigned long num[6];
				cell_strstream >> num[0] >> num[1] >> num[2] >> num[3] >> num[4] >> num[5];
				Cell cell(4, 0);
				cell.at(0) = num[3] - 1;
				cell.at(1) = num[4] - 1;
				cell.at(2) = num[5] - 1;
				cell.at(3) = num[2] - 1;
				C.push_back(cell);
			}
			C.resize(C.size());
			break;
		}
	}

	// read vertices
	unsigned long verticesNum = 0;
	while (file.getline((char *)str.data(), 1024))
	{
		if (str.find("points") != str.npos)
		{
			file.getline((char *)str.data(), 1024);
			stringstream strstream(str.c_str());
			strstream >> verticesNum;
			std::cout << "points : " << verticesNum << std::endl;
			for (unsigned long i = 0; i < verticesNum; i++)
			{
				file.getline((char *)str.data(), 1024);
				stringstream vertex_strstream(str.c_str());
				Vertex v;
				vertex_strstream >> v.x >> v.y >> v.z;
				V.push_back(v);
			}
			V.resize(V.size());
			break;
		}
	}

	MeshFileWriter writer(V, C, argv[2], TETRAHEDRA);
	writer.WriteFile();
	file.close();
}


