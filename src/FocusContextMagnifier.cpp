/*
 * FocusContextMagnifier.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: cotrik
 */

#include "FocusContextMagnifier.h"
#ifdef _WIN32
#include "Eigen/Core"
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Sparse"
//#include "Eigen/MPRealSupport"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "Eigen/SparseCore"
//#include "Eigen/src/MPRealSupport/mpfr/mpfr.h"
#else
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Eigen"
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
//#include "eigen3/Eigen/MPRealSupport"
#include "eigen3/Eigen/Cholesky"
#include "eigen3/Eigen/LU"
#include "eigen3/Eigen/SVD"
#include "eigen3/Eigen/SparseCore"
#endif

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "include.h"

const double W = 1.0;
const double g_delta = 1.0;
FocusContextMagnifier::FocusContextMagnifier(const Mesh& mesh)
: m_mesh((Mesh& )mesh)
{
	// TODO Auto-generated constructor stub

}

FocusContextMagnifier::~FocusContextMagnifier()
{
	// TODO Auto-generated destructor stub
}

const size_t FocusContextMagnifier::GetHexVerticesNum() const
{
	return m_mesh.m_vecPolycubeHexVertex.size();
}

const size_t FocusContextMagnifier::GetHexCellsNum() const
{
	return m_mesh.m_vecPolycubeHexCell.size();
}

const size_t FocusContextMagnifier::GetHexEdgesNum()// const
{
	const size_t verticesNum = GetHexVerticesNum();
//	std::vector<Edge> vecEdge;
	m_vecHexEdge.clear();
	const size_t cellsNum = GetHexCellsNum();
	for (unsigned long i = 0; i < cellsNum; i++)
	{
		const Cell& cell = m_mesh.m_vecPolycubeHexCell.at(i);
		if (cell.size() > 8)
		{
			continue;
		}
		for (unsigned long j = 0; j < 12; j++)
		{
			const unsigned long& p0 = cell.at(HexEdge[j][0]);
			const unsigned long& p1 = cell.at(HexEdge[j][1]);
//			if (p0 < p1)
//			{
				const Edge e(p0, p1);
				if (std::find(m_vecHexEdge.begin(), m_vecHexEdge.end(), e) == m_vecHexEdge.end())
				m_vecHexEdge.push_back(e);
//			}
//			else
//			{
//				const Edge e(p1, p0);
//				if (std::find(m_vecHexEdge.begin(), m_vecHexEdge.end(), e) == m_vecHexEdge.end())
//				m_vecHexEdge.push_back(e);
//			}
		}
	}
	std::sort(m_vecHexEdge.begin(), m_vecHexEdge.end());
	std::vector<Edge>::iterator iter = std::unique(m_vecHexEdge.begin(), m_vecHexEdge.end());
	unsigned long numOfEdges = std::distance(m_vecHexEdge.begin(), iter);
	m_vecHexEdge.resize(numOfEdges);
	return numOfEdges;
}

void FocusContextMagnifier::UpdateVerticesPositions()
{
	ConstructLinearEquation();
	SolveLinearEquation();
}

//void FocusContextMagnifier::SolveLinearEquation()
//{
//	int i = 1;
//	while (i-- > 0){
//	const size_t verticesNum = GetHexVerticesNum();
//	const size_t cellsNum = GetHexCellsNum();
//	const size_t edgesNum = GetHexEdgesNum();
//	std::cout << "number of hex vertices is " << verticesNum << std::endl;
//	std::cout << "number of hex cells is " << cellsNum << std::endl;
//	std::cout << "number of hex edges is " << edgesNum << std::endl;
//
//	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(edgesNum, edgesNum);
//	Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(edgesNum, verticesNum);
//	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(verticesNum, verticesNum);
//	Eigen::MatrixXd I = Eigen::MatrixXd::Zero(verticesNum, verticesNum);
//	Eigen::MatrixXd AA = Eigen::MatrixXd::Zero(2*verticesNum + edgesNum, verticesNum);
//
//	Eigen::VectorXd Hx = Eigen::VectorXd::Zero(verticesNum);
//	Eigen::VectorXd Hy = Eigen::VectorXd::Zero(verticesNum);
//	Eigen::VectorXd Hz = Eigen::VectorXd::Zero(verticesNum);
//	Eigen::VectorXd Vx = Eigen::VectorXd::Zero(verticesNum);
//	Eigen::VectorXd Vy = Eigen::VectorXd::Zero(verticesNum);
//	Eigen::VectorXd Vz = Eigen::VectorXd::Zero(verticesNum);
//
//	Eigen::VectorXd bx = Eigen::VectorXd::Zero(2*verticesNum + edgesNum);
//	Eigen::VectorXd by = Eigen::VectorXd::Zero(2*verticesNum + edgesNum);
//	Eigen::VectorXd bz = Eigen::VectorXd::Zero(2*verticesNum + edgesNum);
//	// S
//	for (unsigned long u = 0; u < edgesNum/2; u++) {
//		for (unsigned long v = 0; v < edgesNum/2; v++) {
//			if (v == u)	{
//				S(u,v) = 1.0;
//			}
////			else {
////				S(u,v) = 0;
////			}
//		}
//	}
//
//	for (unsigned long u = edgesNum/2; u < edgesNum; u++) {
//		for (unsigned long v = edgesNum/2; v < edgesNum; v++) {
//			if (v == u)	{
//				S(u,v) = 1.0;
//			}
////			else {
////				S(u,v) = 0;
////			}
//		}
//	}
//	// Q
//	for (unsigned long u = 0; u < edgesNum; u++) {
//		const Edge& e = m_vecHexEdge.at(u);
//		//const unsigned long& p = u;
//		const unsigned long& i = e.p0;
//		const unsigned long& j = e.p1;
//		for (unsigned long v = 0; v < verticesNum; v++)	{
//			if (/*u == p && */v == i) {
//				Q(u,v) = 1;
//			}
//			else if (/*u == p && */v == j)	{
//				Q(u,v) = -1;
//			}
//			else{
//				Q(u,v) = 0;
//			}
//		}
//	}
/////	std::cout << "Q = \n" << Q << std::endl;
//	// I
//	for (unsigned long u = 0; u < verticesNum; u++)	{
//		for (unsigned long v = 0; v < verticesNum; v++)	{
//			if (v == u)	{
//				I(u,v) = W;
//			}
//			else {
//				I(u,v) = 0;
//			}
//		}
//	}
//	// R
//	for (unsigned long u = 0; u < verticesNum; u++) {
//		for (unsigned long v = 0; v < verticesNum; v++)	{
//			if (v == u)	{
//				R(u, v) = g_delta;
//			}
//			else {
//				R(u, v) = 0;
//			}
//		}
//	}
//	// V
//	for (unsigned long i = 0; i < verticesNum; i++)	{
//		const Vertex& vertex = m_mesh.m_vecPolycubeHexVertex.at(i);
//		Vx(i) = vertex.x;
//		Vy(i) = vertex.y;
//		Vz(i) = vertex.z;
//	}
//	// H
//	for (unsigned long i = 0; i < verticesNum; i++)	{
//		if (Vx(i) < m_mesh.m_minVertex.x)
//		{
//			Hx(i) = g_delta*m_mesh.m_minVertex.x;
//		}
//		else if (Vx(i) > m_mesh.m_maxVertex.x)
//		{
//			Hx(i) = g_delta*m_mesh.m_maxVertex.x;
//		}
//
//		if (Vy(i) < m_mesh.m_minVertex.y)
//		{
//			Hy(i) = g_delta*m_mesh.m_minVertex.y;
//		}
//		else if (Vy(i) > m_mesh.m_maxVertex.y)
//		{
//			Hy(i) = g_delta*m_mesh.m_maxVertex.y;
//		}
//
//		if (Vz(i) < m_mesh.m_minVertex.z)
//		{
//			Hz(i) = g_delta*m_mesh.m_minVertex.z;
//		}
//		else if (Vz(i) > m_mesh.m_maxVertex.z)
//		{
//			Hz(i) = g_delta*m_mesh.m_maxVertex.z;
//		}
//	}
//
//	// A
//	for (unsigned long u = 0; u < edgesNum; u++) {
//		for (unsigned long v = 0; v < verticesNum; v++)	{
//			AA(u, v) = Q(u, v);
//		}
//	}
//	for (unsigned long u = edgesNum; u < edgesNum + verticesNum; u++) {
//		for (unsigned long v = 0; v < verticesNum; v++)	{
//			AA(u, v) = I(u - edgesNum, v);
//		}
//	}
//	for (unsigned long u = edgesNum + verticesNum; u < edgesNum + 2*verticesNum; u++) {
//		for (unsigned long v = 0; v < verticesNum; v++)	{
//			AA(u, v) = R(u - edgesNum - verticesNum, v);
//		}
//	}
//
//	Eigen::MatrixXd AT = AA.transpose();
//	Eigen::MatrixXd A = AT*AA;
////	std::cout << "A = \n" << A << std::endl;
//	// b
//	Eigen::VectorXd SQPx = Vx.transpose()*Q.transpose()*S.transpose();
//	for (unsigned long i = 0; i < edgesNum; i++) {
//		bx(i) = SQPx(i);
//	}
//	for (unsigned long i = edgesNum; i < edgesNum + verticesNum; i++) {
//		bx(i) = Vx(i - edgesNum);
//	}
//	for (unsigned long i = edgesNum + verticesNum; i < edgesNum + 2*verticesNum; i++) {
//		bx(i) = Hx(i - edgesNum - verticesNum);
//	}
//
//	Eigen::VectorXd SQPy = Vy.transpose()*Q.transpose()*S.transpose();
////	std::cout << "SQPy = \n" << SQPy << std::endl;
//	for (unsigned long i = 0; i < edgesNum; i++) {
//		by(i) = SQPy(i);
//	}
//	for (unsigned long i = edgesNum; i < edgesNum + verticesNum; i++) {
//		by(i) = Vy(i - edgesNum);
//	}
//	for (unsigned long i = edgesNum + verticesNum; i < edgesNum + 2*verticesNum; i++) {
//		by(i) = Hy(i - edgesNum - verticesNum);
//	}
//
//	Eigen::VectorXd SQPz = Vz.transpose()*Q.transpose()*S.transpose();
//	for (unsigned long i = 0; i < edgesNum; i++) {
//		bz(i) = SQPz(i);
//	}
//	for (unsigned long i = edgesNum; i < edgesNum + verticesNum; i++) {
//		bz(i) = Vz(i - edgesNum);
//	}
//	for (unsigned long i = edgesNum + verticesNum; i < edgesNum + 2*verticesNum; i++) {
//		bz(i) = Hz(i - edgesNum - verticesNum);
//	}
//
//	Eigen::VectorXd Bx = AT*bx;
//	Eigen::VectorXd By = AT*by;
//	Eigen::VectorXd Bz = AT*bz;
//
//	Eigen::VectorXd x = A.llt().solve(Bx);
//	Eigen::VectorXd y = A.llt().solve(By);
//	Eigen::VectorXd z = A.llt().solve(Bz);
//
//	for (unsigned long i = 0; i < verticesNum; i++) {
//		Vertex& vertex = m_mesh.m_vecPolycubeHexVertex.at(i);
//		vertex.x = x(i);
//		vertex.y = y(i);
//		vertex.z = z(i);
//	}
//	}
//}
static int iter_time = 0;
// kitty
//const unsigned int enlargeCubeIndices[] =
//{
//    //11, 18, 41, 49
//	428, 597, 415, 584, // right ear
//	427, 596, 414, 583, // right ear
//
//	1413, 1426, 1439,   // left ear
//	1582, 1595, 1608,   // left ear
//
//	1412, 1425, 1438,   // left ear
//	1581, 1594, 1607,   // left ear
//
//	1411, 1424, 1437,   // left ear
//	1580, 1593, 1606,   // left ear
//	1410,
//
//	1245, 1258, 1271,
//	1246, 1259, 1272,
//
//	426, 595, 401, 570, 608,  // right ear
//
//	1244, 1258, 1271,
//};
// rockerArm_1
//const unsigned int enlargeCubeIndices[] =
//{
//	// 4 times larger
//	753, 766,
//	752, 765,
//	922, 935, 948,
//	921, 934, 947,
//	// 3 times larger
//	920, 933, 946,
//	1091, 1104, 1117,
//	1090, 1103, 1116,
//	1089, 1102, 1115,
//	// 2 times larger
//    1260, 1273, 1286,
//    1259, 1272, 1285,
//    1258, 1271, 1284
//};

// bunny
/*const unsigned int enlargeCubeIndices[] =
{
	// 4 times larger
	// left ear
	1338, 1156, 1325,
	1337, 1155, 1324,
	1142, 1311,
	1141, 1310, 1323,
	// right ear
	1481, 1650, 1468, 1637,
	1480, 1649, 1647, 1636,

	// 3 times larger
    1129, 1298, 1116, 1285,
    1128, 1297,
	// 2 times larger
    1115, 1284
};*/

// canwt
const unsigned int enlargeCubeIndices[] =
{
	// 4 times larger
	// right front toe
	1454, 1441, 1428,
	1285, 1272, 1259,
	1116, 1103, 1090,

	1453, 1440, 1427,
	1284, 1271, 1258,
	1115, 1102, 1089,

	// left front toe
	752,
	596, 583,

	427, 414, 401,
	764, 751, 738,
	595, 582, 569,
	426, 413, 400,
	242, 231,

	// right ear
	1118, 1287, 1131, 1300,

	// left ear
	638, 807, 651, 820,
	639, 808, 652, 821,
	// 3 times larger
	// right ear
    1119, 1288, 1132, 1301
	// 2 times larger
};

const IndexScale indexScale[] =
{
	{38,   1.00},
	{51,   1.00},
	{64,   1.00},
	{77,   1.00},
	{90,   1.00},
	{103,  1.00},
	{116,  1.00},
	{129,  1.00},
	{207,  1.00},
	{220,  1.00},
	{233,  1.10},    ///
	{246,  1.10},    ///
	{259,  1.10},    ///
	{272,  1.10},    ///
	{285,  1.00},
	{298,  1.00},
	{376,  1.00},
	{389,  1.00},
	{402,  1.10},    ///
	{415,  10.20},    //
	{428,  10.20},    //
	{441,  1.10},    ///
	{454,  1.00},
	{467,  1.00},
	{545,  1.00},
	{558,  1.00},
	{571,  1.10},    ///
	{584,  10.20},    //
	{597,  10.20},    //
	{610,  1.10},    ///
	{623,  1.00},
	{636,  1.00},
	{714,  1.00},
	{727,  1.00},
	{740,  1.10},    ///
	{753,  1.10},    ///
	{766,  1.10},    ///
	{779,  1.10},    ///
	{792,  1.00},
	{805,  1.00},
	{883,  1.00},
	{896,  1.00},
	{909,  1.00},
	{922,  1.00},
	{935,  1.00},
	{948,  1.00},
	{961,  1.00},
	{974,  1.00},
	{1052, 1.00},
	{1065, 1.00},
	{1078, 1.00},
	{1091, 1.00},
	{1104, 1.00},
	{1117, 1.00},
	{1130, 1.00},
	{1143, 1.00},
	{1221, 1.00},
	{1234, 1.00},
	{1247, 1.00},
	{1260, 1.00},
	{1273, 1.00},
	{1286, 1.00},
	{1299, 1.00},
	{1312, 1.00},

    {37,   1.00},
    {50,   1.00},
    {63,   1.00},
    {76,   1.00},
    {89,   1.00},
    {102,  1.00},
    {115,  1.00},
    {128,  1.00},
    {206,  1.00},
    {219,  1.00},
    {232,  1.00},
    {245,  1.00},
    {258,  1.00},
    {271,  1.00},
    {284,  1.00},
    {297,  1.00},
    {375,  1.00},
    {388,  1.00},
    {401,  1.20},    ///
    {414,  10.20},    //
    {427,  10.20},    //
    {440,  1.20},    ///
    {453,  1.00},
    {466,  1.00},
    {544,  1.00},
    {557,  1.00},
    {570,  1.20},    ///
    {583,  10.20},    //
    {596,  10.20},    //
    {609,  1.20},    ///
    {622,  1.00},
    {635,  1.00},
    {713,  1.00},
    {726,  1.00},
    {739,  1.00},
    {752,  1.20},   ///
    {765,  1.20},   ///
    {778,  1.00},
    {791,  1.00},
    {804,  1.00},
    {882,  1.00},
    {895,  1.00},
    {908,  1.00},
    {921,  1.00},
    {934,  1.00},
    {947,  1.00},
    {960,  1.00},
    {973,  1.00},
    {1051, 1.00},
    {1064, 1.00},
    {1077, 1.00},
    {1090, 1.00},
    {1103, 1.00},
    {1116, 1.00},
    {1129, 1.00},
    {1142, 1.00},
    {1220, 1.00},
    {1233, 1.00},
    {1246, 1.00},
    {1259, 1.00},
    {1272, 1.00},
    {1285, 1.00},
    {1298, 1.00},
    {1311, 1.00},

    {36,   1.00},
    {49,   1.00},
    {62,   1.00},
    {75,   1.00},
    {88,   1.00},
    {101,  1.00},
    {114,  1.00},
    {127,  1.00},
    {205,  1.00},
    {218,  1.00},
    {231,  1.00},
    {244,  1.00},
    {257,  1.00},
    {270,  1.00},
    {283,  1.00},
    {296,  1.00},
    {374,  1.00},
    {387,  1.00},
    {400,  1.00},
    {413,  1.00},
    {426,  1.10},     ///
    {439,  1.00},
    {452,  1.00},
    {465,  1.00},
    {543,  1.00},
    {556,  1.00},
    {569,  1.10},
    {582,  1.10},
    {595,  1.10},     ///
    {608,  1.00},
    {621,  1.00},
    {634,  1.00},
    {712,  1.00},
    {725,  1.00},
    {738,  1.10},
    {751,  1.10},
    {764,  1.10},
    {777,  1.00},
    {790,  1.00},
    {803,  1.00},
    {881,  1.00},
    {894,  1.00},
    {907,  1.00},
    {920,  1.00},
    {933,  1.00},
    {946,  1.00},
    {959,  1.00},
    {972,  1.00},
    {1050, 1.00},
    {1063, 1.00},
    {1076, 1.00},
    {1089, 1.00},
    {1102, 1.00},
    {1115, 1.00},
    {1128, 1.00},
    {1141, 1.00},
    {1219, 1.00},
    {1232, 1.00},
    {1245, 1.00},
    {1258, 1.00},
    {1271, 1.00},
    {1284, 1.00},
    {1297, 1.00},
    {1310, 1.00},

    {35,   1.00},
    {48,   1.00},
    {61,   1.00},
    {74,   1.00},
    {87,   1.00},
    {100,  1.00},
    {113,  1.00},
    {126,  1.00},
    {204,  1.00},
    {217,  1.00},
    {230,  1.00},
    {243,  1.00},
    {256,  1.00},
    {269,  1.00},
    {282,  1.00},
    {295,  1.00},
    {373,  1.00},
    {386,  1.00},
    {399,  1.00},
    {412,  1.00},
    {425,  1.00},
    {438,  1.00},
    {451,  1.00},
    {464,  1.00},
    {542,  1.00},
    {555,  1.00},
    {568,  1.00},
    {581,  1.00},
    {594,  1.00},
    {607,  1.00},
    {620,  1.00},
    {633,  1.00},
    {711,  1.00},
    {724,  1.00},
    {737,  1.00},
    {750,  1.00},
    {763,  1.00},
    {776,  1.00},
    {789,  1.00},
    {802,  1.00},
    {880,  1.00},
    {893,  1.00},
    {906,  1.00},
    {919,  1.00},
    {932,  1.00},
    {945,  1.00},
    {958,  1.00},
    {971,  1.00},
    {1049, 1.00},
    {1062, 1.00},
    {1075, 1.00},
    {1088, 1.00},
    {1101, 1.00},
    {1114, 1.00},
    {1127, 1.00},
    {1140, 1.00},
    {1218, 1.00},
    {1231, 1.00},
    {1244, 1.00},
    {1257, 1.00},
    {1270, 1.00},
    {1283, 1.00},
    {1296, 1.00},
    {1309, 1.00}
};
std::string strScaleParas;
#include "MeshDisplayer.h"
void ReadIndexScale(std::vector<IndexScale> &indexScale)
{
	std::ifstream scaleFileStream(strScaleParas.c_str());
	IndexScale index_scale;
	while (scaleFileStream >> index_scale.index)
	{
		scaleFileStream >> index_scale.scale;
		indexScale.push_back(index_scale);
	}
	scaleFileStream.close();

	if (indexScale.empty())
	{
		if (SELECTED_CUBE_INDEX == 0)
			return;
		// layer 1
		if (SELECTED_CUBE_INDEX != 0 && SELECTED_CUBE_INDEX < 2197 && SELECTED_CUBE_INDEX > 0)
		{
			/*
			 * ----------------------
			 * |      |      |      |
			 * |c - 14|c - 13|c - 12|
			 * |      |      |      |
			 * ----------------------
			 * |      |      |      |
			 * |c - 1 |   c  | c + 1|
			 * |      |      |      |
			 * ----------------------
			 * |      |      |      |
			 * |c + 12|c + 13|c + 14|
			 * |      |      |      |
			 * ----------------------
			 */
			int c = SELECTED_CUBE_INDEX - 169;
			int index = c;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 14;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 13;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 12;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 1;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 1;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 12;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 13;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 14;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}

			c = SELECTED_CUBE_INDEX;
			index = c;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 14;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 13;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 12;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 1;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 1;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 12;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 13;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 14;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}

			c = SELECTED_CUBE_INDEX + 169;
			index = c;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 14;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 13;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 12;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c - 1;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 1;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 12;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 13;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
			index = c + 14;
			if (index >= 0 && index < 2197){
			IndexScale i_s(index, 1.1);
			indexScale.push_back(i_s);}
		}
		if (SELECTED_CUBE_INDEX != 0 && SELECTED_CUBE_INDEX < 2197 && SELECTED_CUBE_INDEX > 0)
		{
			IndexScale i_s(SELECTED_CUBE_INDEX, 4);
			indexScale.push_back(i_s);
		}
	}
}

void FocusContextMagnifier::SolveLinearEquation()
{
	const size_t verticesNum = GetHexVerticesNum();
	const size_t cellsNum = GetHexCellsNum();
	const size_t edgesNum = GetHexEdgesNum();
	std::cout << "number of hex vertices is " << verticesNum << std::endl;
	std::cout << "number of hex cells is " << cellsNum << std::endl;
	std::cout << "number of hex edges is " << edgesNum << std::endl;

	const unsigned int n = round(pow(verticesNum, 1.0/3.0)) - 1;
	const unsigned int N = 12*n*n*n;
	const unsigned int M = verticesNum;

	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(N, N);
	Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(N, M);
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(M, M);
	Eigen::MatrixXd I = Eigen::MatrixXd::Zero(M, M);
	Eigen::MatrixXd AA = Eigen::MatrixXd::Zero(N + 2*M, M);

	Eigen::VectorXd Hx = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Hy = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Hz = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Vx = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Vy = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Vz = Eigen::VectorXd::Zero(M);

	Eigen::VectorXd bx = Eigen::VectorXd::Zero(N + 2*M);
	Eigen::VectorXd by = Eigen::VectorXd::Zero(N + 2*M);
	Eigen::VectorXd bz = Eigen::VectorXd::Zero(N + 2*M);
	// S
//	for (unsigned long u = 0; u < N; u++) {
//		for (unsigned long v = 0; v < N; v++) {
//			if (v == u)	{
//				S(u,v) = 1.0;
//			}
//			else {
//				S(u,v) = 0;
//			}
//		}
//	}

//	std::vector<Vertex> vec_vertex;
//	for (int i = 0; i < sizeof(enlargeCubeIndices)/sizeof(enlargeCubeIndices[0]); i++)
//	{
//		const Cell& cell = m_mesh.m_vecPolycubeHexCell.at(enlargeCubeIndices[i]);
//		for (int c = 0; c < 8; c++)
//		{
//			vec_vertex.push_back(m_mesh.m_vecPolycubeHexVertex[cell[c]]);
//		}
//	}
//	Vertex v_sum;
//	for (int i = 0; i < vec_vertex.size(); i++)
//	{
//		v_sum.x += vec_vertex.at(i).x;
//		v_sum.y += vec_vertex.at(i).y;
//		v_sum.z += vec_vertex.at(i).z;
//	}
//	Vertex v_center(v_sum.x/vec_vertex.size(), v_sum.y/vec_vertex.size(), v_sum.z/vec_vertex.size());
//
//	std::vector<float> vec_scale;
//	for (int i = 0; i < m_mesh.m_vecPolycubeHexCell.size()/*sizeof(enlargeCubeIndices)/sizeof(enlargeCubeIndices[0])*/; i++)
//	{
//		const Cell& cell = m_mesh.m_vecPolycubeHexCell.at(i/*enlargeCubeIndices[i]*/);
//		glm::vec3 centerPoint(0.0f, 0.0f, 0.0f);
//		size_t cellSize = cell.size();
//		for (unsigned int c = 0; c < cellSize; c++)
//		{
//			centerPoint.x += m_mesh.m_vecPolycubeHexVertex[cell.at(c)].x;
//			centerPoint.y += m_mesh.m_vecPolycubeHexVertex[cell.at(c)].y;
//			centerPoint.z += m_mesh.m_vecPolycubeHexVertex[cell.at(c)].z;
//		}
//		centerPoint.x = centerPoint.x / cellSize;
//		centerPoint.y = centerPoint.y / cellSize;
//		centerPoint.z = centerPoint.z / cellSize;
//		glm::vec3 targetPoint(v_center.x, v_center.y, v_center.z);
//		glm::vec3 dir = centerPoint - targetPoint;
//		float max = fabs(dir.z);
//		if (fabs(dir.y) > max)
//			max = fabs(dir.y);
//		if (fabs(dir.z) > max)
//			max = fabs(dir.z);
//		int numberOfGrids = max/0.08;
//		float scale = 1.5 - numberOfGrids*0.06;
//		vec_scale.push_back(scale<1.0f ? 1.0f : scale);
////		float distance = glm::length(dir);
////		if (distance < 1.0/10)
////		{
////			vec_scale.push_back(1.5);
////		}
////		else
////			vec_scale.push_back(2.5/(10*distance));
//	}
//
//	for (int i = 0; i < m_mesh.m_vecPolycubeHexCell.size()/*sizeof(enlargeCubeIndices)/sizeof(enlargeCubeIndices[0])*/; i++)
//	{
//		const unsigned int& index = i;
//		for (unsigned long u = index*12; u < index*12 + 12; u++) {
//			for (unsigned long v = index*12; v < index*12 + 12; v++) {
//				if (v == u)	{
//					S(u,v) = vec_scale[i];
//				}
//			}
//		}
//	}

//	for (int i = 0; i < sizeof(enlargeCubeIndices)/sizeof(enlargeCubeIndices[0]); i++)
//	{
//		const unsigned int& index = enlargeCubeIndices[i];
//		for (unsigned long u = index*12; u < index*12 + 12; u++) {
//			for (unsigned long v = index*12; v < index*12 + 12; v++) {
//				if (v == u)	{
//					S(u,v) = 4.0;
//				}
//			}
//		}
//	}
//	for (int i = sizeof(enlargeCubeIndices)/sizeof(enlargeCubeIndices[0]) - 4; i < sizeof(enlargeCubeIndices)/sizeof(enlargeCubeIndices[0]); i++)
//	{
//		const unsigned int& index = enlargeCubeIndices[i];
//		for (unsigned long u = index*12; u < index*12 + 12; u++) {
//			for (unsigned long v = index*12; v < index*12 + 12; v++) {
//				if (v == u)	{
//					S(u,v) = 3.0;
//				}
//			}
//		}
//	}
/*
 	*/
//	for (unsigned long u = 0; u < N/2; u++) {
//		for (unsigned long v = 0; v < N/2; v++) {
//			if (v == u)	{
//				S(u,v) = 1.0;
//			}
//			else {
//				S(u,v) = 0;
//			}
//		}
//	}
//

//	for (int i = 0; i < sizeof(indexScale)/sizeof(indexScale[0]); i++)
//	{
//		const unsigned int& index = indexScale[i].index;
//		for (unsigned long u = index*12; u < index*12 + 12; u++) {
//			for (unsigned long v = index*12; v < index*12 + 12; v++) {
//				if (v == u)	{
//					S(u,v) = indexScale[i].scale;//1.0;//
//				}
//			}
//		}
//	}

	// S
//	for (unsigned long u = 0; u < N; u++) {
//		for (unsigned long v = 0; v < N; v++) {
//			if (v == u)	{
//				S(u,v) = 20.0;
//			}
//		}
//	}
//	for (unsigned long u = N/2; u < N; u++) {
//		for (unsigned long v = N/2; v < N; v++) {
//			if (v == u)	{
//				S(u,v) = 2.0;
//			}
//			else {
//				S(u,v) = 0;
//			}
//		}
//	}
//
//	for (unsigned long u = N*(1 + iter_time)/(2 + iter_time); u < N; u++) {
//		for (unsigned long v = N*(1 + iter_time)/(2 + iter_time); v < N; v++) {
//			if (v == u)	{
//				S(u,v) = 1.05 + 0.02*iter_time;
//			}
//			else {
//				S(u,v) = 0;
//			}
//		}
//	}
//	iter_time++;

	// S
	for (unsigned long u = 0; u < N; u++) {
		for (unsigned long v = 0; v < N; v++) {
			if (v == u)	{
				S(u,v) = 1.0;
			}
		}
	}
	std::vector<IndexScale> indexScale;
	ReadIndexScale(indexScale);
	for (int i = 0; i < indexScale.size(); i++)
	{
		const unsigned int& index = indexScale[i].index;
		for (unsigned long u = index * 12; u < index * 12 + 12; u++)
		{
			for (unsigned long v = index * 12; v < index * 12 + 12; v++)
			{
				if (v == u)
				{
					S(u, v) = indexScale[i].scale;
				}
			}
		}
	}
	// Q
	for (unsigned long u = 0; u < N/12; u++) {
		const Cell& cell = m_mesh.m_vecPolycubeHexCell.at(u);
		if (cell.size() > 8)
			continue;
		for (unsigned int eIndex = 0; eIndex < 12; eIndex++) {
			const unsigned long& i = cell.at(HexEdge[eIndex][0]);
			const unsigned long& j = cell.at(HexEdge[eIndex][1]);
		    Q(u*12 + eIndex, i) = 1;
		    Q(u*12 + eIndex, j) = -1;
		}
	}
    //	std::cout << "Q = \n" << Q << std::endl;
	// I
	for (unsigned long u = 0; u < M; u++)	{
		for (unsigned long v = 0; v < M; v++)	{
			if (v == u)	{
				I(u,v) = W;
			}
//			else {
//				I(u,v) = 0;
//			}
		}
	}
	// R
//	for (unsigned long u = 0; u < M; u++) {
//		for (unsigned long v = 0; v < M; v++)	{
//			if (v == u)	{
//				R(u, v) = 0;
//			}
//			else {
//				R(u, v) = 0;
//			}
//		}
//	}
	// V
	for (unsigned long i = 0; i < M; i++)	{
		const Vertex& vertex = m_mesh.m_vecPolycubeHexVertex.at(i);
		Vx(i) = vertex.x;
		Vy(i) = vertex.y;
		Vz(i) = vertex.z;
	}
	// H
	for (unsigned long i = 0; i < M; i++)	{
		if (Vx(i) < m_mesh.m_minVertex.x)
		{
			Hx(i) = g_delta*m_mesh.m_minVertex.x;
		}
		else if (Vx(i) > m_mesh.m_maxVertex.x)
		{
			Hx(i) = g_delta*m_mesh.m_maxVertex.x;
		}

		if (Vy(i) < m_mesh.m_minVertex.y)
		{
			Hy(i) = g_delta*m_mesh.m_minVertex.y;
		}
		else if (Vy(i) > m_mesh.m_maxVertex.y)
		{
			Hy(i) = g_delta*m_mesh.m_maxVertex.y;
		}

		if (Vz(i) < m_mesh.m_minVertex.z)
		{
			Hz(i) = g_delta*m_mesh.m_minVertex.z;
		}
		else if (Vz(i) > m_mesh.m_maxVertex.z)
		{
			Hz(i) = g_delta*m_mesh.m_maxVertex.z;
		}
	}

	// A
	for (unsigned long u = 0; u < N; u++) {
		for (unsigned long v = 0; v < M; v++)	{
			AA(u, v) = Q(u, v);
		}
	}
	for (unsigned long u = N; u < N + M; u++) {
		for (unsigned long v = 0; v < M; v++)	{
			AA(u, v) = I(u - N, v);
		}
	}
	for (unsigned long u = N + M; u < N + 2*M; u++) {
		for (unsigned long v = 0; v < M; v++)	{
			AA(u, v) = R(u - N - M, v);
		}
	}

	Eigen::MatrixXd AT = AA.transpose();
	Eigen::MatrixXd A = AT*AA;
//	std::cout << "A = \n" << A << std::endl;
	// b
	Eigen::VectorXd SQPx = Vx.transpose()*Q.transpose()*S.transpose();
	for (unsigned long i = 0; i < N; i++) {
		bx(i) = SQPx(i);
	}
	for (unsigned long i = N; i < N + M; i++) {
		bx(i) = Vx(i - N);
	}
	for (unsigned long i = N + M; i < N + 2*M; i++) {
		bx(i) = Hx(i - N - M);
	}

	Eigen::VectorXd SQPy = Vy.transpose()*Q.transpose()*S.transpose();
//	std::cout << "SQPy = \n" << SQPy << std::endl;
	for (unsigned long i = 0; i < N; i++) {
		by(i) = SQPy(i);
	}
	for (unsigned long i = N; i < N + M; i++) {
		by(i) = Vy(i - N);
	}
	for (unsigned long i = N + M; i < N + 2*M; i++) {
		by(i) = Hy(i - N - M);
	}

	Eigen::VectorXd SQPz = Vz.transpose()*Q.transpose()*S.transpose();
	for (unsigned long i = 0; i < N; i++) {
		bz(i) = SQPz(i);
	}
	for (unsigned long i = N; i < N + M; i++) {
		bz(i) = Vz(i - N);
	}
	for (unsigned long i = N + M; i < N + 2*M; i++) {
		bz(i) = Hz(i - N - M);
	}

	Eigen::VectorXd Bx = AT*bx;
	Eigen::VectorXd By = AT*by;
	Eigen::VectorXd Bz = AT*bz;

	Eigen::VectorXd x = A.llt().solve(Bx);
	Eigen::VectorXd y = A.llt().solve(By);
	Eigen::VectorXd z = A.llt().solve(Bz);

	for (unsigned long i = 0; i < M; i++) {
		Vertex& vertex = m_mesh.m_vecPolycubeHexVertex.at(i);
		vertex.x = x(i);
		vertex.y = y(i);
		vertex.z = z(i);
	}
}

void FocusContextMagnifier::SolveLinearEquation_TestLLT_PR_LU()
{
	const size_t verticesNum = GetHexVerticesNum();
	const size_t cellsNum = GetHexCellsNum();
	const size_t edgesNum = GetHexEdgesNum();
	std::cout << "number of hex vertices is " << verticesNum << std::endl;
	std::cout << "number of hex cells is " << cellsNum << std::endl;
	std::cout << "number of hex edges is " << edgesNum << std::endl;

	const unsigned int n = round(pow(verticesNum, 1.0/3.0)) - 1;
	const unsigned int N = 12*n*n*n;
	const unsigned int M = verticesNum;

	Eigen::MatrixXd S = Eigen::MatrixXd::Zero(N, N);
	Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(N, M);
	Eigen::MatrixXd R = Eigen::MatrixXd::Zero(M, M);
	Eigen::MatrixXd I = Eigen::MatrixXd::Zero(M, M);
	Eigen::MatrixXd AA = Eigen::MatrixXd::Zero(N + 2*M, M);

	Eigen::VectorXd Hx = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Hy = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Hz = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Vx = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Vy = Eigen::VectorXd::Zero(M);
	Eigen::VectorXd Vz = Eigen::VectorXd::Zero(M);

	Eigen::VectorXd bx = Eigen::VectorXd::Zero(N + 2*M);
	Eigen::VectorXd by = Eigen::VectorXd::Zero(N + 2*M);
	Eigen::VectorXd bz = Eigen::VectorXd::Zero(N + 2*M);

	// S
	for (unsigned long u = 0; u < N; u++) {
		for (unsigned long v = 0; v < N; v++) {
			if (v == u)	{
				S(u,v) = 1.0;
			}
		}
	}
	std::vector<IndexScale> indexScale;
	ReadIndexScale(indexScale);
	for (int i = 0; i < indexScale.size(); i++)
	{
		const unsigned int& index = indexScale[i].index;
		for (unsigned long u = index * 12; u < index * 12 + 12; u++)
		{
			for (unsigned long v = index * 12; v < index * 12 + 12; v++)
			{
				if (v == u)
				{
					S(u, v) = indexScale[i].scale;
				}
			}
		}
	}
	// Q
	for (unsigned long u = 0; u < N/12; u++) {
		const Cell& cell = m_mesh.m_vecPolycubeHexCell.at(u);
		if (cell.size() > 8)
			continue;
		for (unsigned int eIndex = 0; eIndex < 12; eIndex++) {
			const unsigned long& i = cell.at(HexEdge[eIndex][0]);
			const unsigned long& j = cell.at(HexEdge[eIndex][1]);
		    Q(u*12 + eIndex, i) = 1;
		    Q(u*12 + eIndex, j) = -1;
		}
	}
    //	std::cout << "Q = \n" << Q << std::endl;
	// I
	for (unsigned long u = 0; u < M; u++)	{
		for (unsigned long v = 0; v < M; v++)	{
			if (v == u)	{
				I(u,v) = W;
			}
		}
	}
	// V
	for (unsigned long i = 0; i < M; i++)	{
		const Vertex& vertex = m_mesh.m_vecPolycubeHexVertex.at(i);
		Vx(i) = vertex.x;
		Vy(i) = vertex.y;
		Vz(i) = vertex.z;
	}
	// H
	for (unsigned long i = 0; i < M; i++)	{
		if (Vx(i) < m_mesh.m_minVertex.x)
		{
			Hx(i) = g_delta*m_mesh.m_minVertex.x;
		}
		else if (Vx(i) > m_mesh.m_maxVertex.x)
		{
			Hx(i) = g_delta*m_mesh.m_maxVertex.x;
		}

		if (Vy(i) < m_mesh.m_minVertex.y)
		{
			Hy(i) = g_delta*m_mesh.m_minVertex.y;
		}
		else if (Vy(i) > m_mesh.m_maxVertex.y)
		{
			Hy(i) = g_delta*m_mesh.m_maxVertex.y;
		}

		if (Vz(i) < m_mesh.m_minVertex.z)
		{
			Hz(i) = g_delta*m_mesh.m_minVertex.z;
		}
		else if (Vz(i) > m_mesh.m_maxVertex.z)
		{
			Hz(i) = g_delta*m_mesh.m_maxVertex.z;
		}
	}

	// A
	for (unsigned long u = 0; u < N; u++) {
		for (unsigned long v = 0; v < M; v++)	{
			AA(u, v) = Q(u, v);
		}
	}
	for (unsigned long u = N; u < N + M; u++) {
		for (unsigned long v = 0; v < M; v++)	{
			AA(u, v) = I(u - N, v);
		}
	}
	for (unsigned long u = N + M; u < N + 2*M; u++) {
		for (unsigned long v = 0; v < M; v++)	{
			AA(u, v) = R(u - N - M, v);
		}
	}

	Eigen::MatrixXd AT = AA.transpose();
	Eigen::MatrixXd A = AT*AA;
//	std::cout << "A = \n" << A << std::endl;
	// b
	Eigen::VectorXd SQPx = Vx.transpose()*Q.transpose()*S.transpose();
	for (unsigned long i = 0; i < N; i++) {
		bx(i) = SQPx(i);
	}
	for (unsigned long i = N; i < N + M; i++) {
		bx(i) = Vx(i - N);
	}
	for (unsigned long i = N + M; i < N + 2*M; i++) {
		bx(i) = Hx(i - N - M);
	}

	Eigen::VectorXd SQPy = Vy.transpose()*Q.transpose()*S.transpose();
//	std::cout << "SQPy = \n" << SQPy << std::endl;
	for (unsigned long i = 0; i < N; i++) {
		by(i) = SQPy(i);
	}
	for (unsigned long i = N; i < N + M; i++) {
		by(i) = Vy(i - N);
	}
	for (unsigned long i = N + M; i < N + 2*M; i++) {
		by(i) = Hy(i - N - M);
	}

	Eigen::VectorXd SQPz = Vz.transpose()*Q.transpose()*S.transpose();
	for (unsigned long i = 0; i < N; i++) {
		bz(i) = SQPz(i);
	}
	for (unsigned long i = N; i < N + M; i++) {
		bz(i) = Vz(i - N);
	}
	for (unsigned long i = N + M; i < N + 2*M; i++) {
		bz(i) = Hz(i - N - M);
	}

	Eigen::VectorXd Bx = AT*bx;
	Eigen::VectorXd By = AT*by;
	Eigen::VectorXd Bz = AT*bz;

	std::cout << "#############################################" << std::endl;
	std::cout << "A: " << A.rows() << " X " << A.cols() << std::endl;
	Timer timerllt = Timer();
	timerllt.start();
	Eigen::VectorXd x = A.llt().solve(Bx);
	Eigen::VectorXd y = A.llt().solve(By);
	Eigen::VectorXd z = A.llt().solve(Bz);
	double durationllt = timerllt.stop();
	std::cout << "LLT SOLVER: ";
	timerllt.printTime(durationllt/3);

	Timer timerqr = Timer();
	timerqr.start();
	Eigen::VectorXd xqr = A.fullPivHouseholderQr().solve(Bx);
	Eigen::VectorXd yqr = A.fullPivHouseholderQr().solve(By);
	Eigen::VectorXd zqr = A.fullPivHouseholderQr().solve(Bz);
	double durationqr = timerqr.stop();
	std::cout << "QR SOLVER: ";
	timerqr.printTime(durationqr/3);

	Timer timerlu = Timer();
	timerlu.start();
	Eigen::VectorXd xlu = A.lu().solve(Bx);
	Eigen::VectorXd ylu = A.lu().solve(By);
	Eigen::VectorXd zlu = A.lu().solve(Bz);
	double durationlu = timerlu.stop();
	std::cout << "LU SOLVER: ";
	timerlu.printTime(durationlu/3);

	Timer timerFullPivLU = Timer();
	timerFullPivLU.start();
	Eigen::VectorXd xFullPivLU = A.fullPivLu().solve(Bx);
	Eigen::VectorXd yFullPivLU = A.fullPivLu().solve(By);
	Eigen::VectorXd zFullPivLU = A.fullPivLu().solve(Bz);
	double durationFullPivLU = timerFullPivLU.stop();
	std::cout << "LU SOLVER: ";
	timerFullPivLU.printTime(durationFullPivLU/3);

	double errorllt = 0;
	double errorqr = 0;
	double errorlu = 0;
	for (unsigned long i = 0; i < M; i++) {
		errorllt += (x(i) - xFullPivLU(i))*(x(i) - xFullPivLU(i)) + (y(i) - yFullPivLU(i))*(y(i) - yFullPivLU(i)) + (z(i) - zFullPivLU(i))*(z(i) - zFullPivLU(i));
		errorlu += (xlu(i) - xFullPivLU(i))*(xlu(i) - xFullPivLU(i)) + (ylu(i) - yFullPivLU(i))*(ylu(i) - yFullPivLU(i)) + (zlu(i) - zFullPivLU(i))*(zlu(i) - zFullPivLU(i));
		errorqr += (xqr(i) - xFullPivLU(i))*(xqr(i) - xFullPivLU(i)) + (yqr(i) - yFullPivLU(i))*(yqr(i) - yFullPivLU(i)) + (zqr(i) - zFullPivLU(i))*(zqr(i) - zFullPivLU(i));
	}

	std::cout << "LLT SOLVER ERROR: " << errorllt << std::endl;
	std::cout << "QR SOLVER ERROR: " << errorqr << std::endl;
	std::cout << "LU SOLVER ERROR: " << errorlu << std::endl;
	std::cout << "#############################################" << std::endl;
	for (unsigned long i = 0; i < M; i++) {
		Vertex& vertex = m_mesh.m_vecPolycubeHexVertex.at(i);
		vertex.x = x(i);
		vertex.y = y(i);
		vertex.z = z(i);
	}
}

void FocusContextMagnifier::ConstructLinearEquation()
{
//	Construct_Q();
//	Construct_Rx();
//	Construct_Ry();
//	Construct_Rz();
//
//	Construct_S();
//	Construct_Vx();
//	Construct_Vy();
//	Construct_Vz();
//	Construct_Hx();
//	Construct_Hy();
//	Construct_Hz();
//
//	Construct_A();
//	Construct_x();
//	Construct_b();
}

void FocusContextMagnifier::Construct_A()
{

}

void FocusContextMagnifier::Construct_b()
{

}

void FocusContextMagnifier::Construct_Q()
{

}

