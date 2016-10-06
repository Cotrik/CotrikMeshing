/*
 * gen_rot_meshes.cpp
 *
 *  Created on: Nov 26, 2015
 *      Author: cotrik
 */

#include "include.h"

const char* out_x_filenames[35] =
{
	"_x_5.off",
	"_x_10.off",
	"_x_15.off",
	"_x_20.off",
	"_x_25.off",
	"_x_30.off",
	"_x_35.off",
	"_x_40.off",
	"_x_45.off",
	"_x_50.off",
	"_x_55.off",
	"_x_60.off",
	"_x_65.off",
	"_x_70.off",
	"_x_75.off",
	"_x_80.off",
	"_x_85.off",
	"_x_90.off",
	"_x_95.off",
	"_x_100.off",
	"_x_105.off",
	"_x_110.off",
	"_x_115.off",
	"_x_120.off",
	"_x_125.off",
	"_x_130.off",
	"_x_135.off",
	"_x_140.off",
	"_x_145.off",
	"_x_150.off",
	"_x_155.off",
	"_x_160.off",
	"_x_165.off",
	"_x_170.off",
	"_x_175.off"
};

const char* out_y_filenames[35] =
{
	"_y_5.off",
	"_y_10.off",
	"_y_15.off",
	"_y_20.off",
	"_y_25.off",
	"_y_30.off",
	"_y_35.off",
	"_y_40.off",
	"_y_45.off",
	"_y_50.off",
	"_y_55.off",
	"_y_60.off",
	"_y_65.off",
	"_y_70.off",
	"_y_75.off",
	"_y_80.off",
	"_y_85.off",
	"_y_90.off",
	"_y_95.off",
	"_y_100.off",
	"_y_105.off",
	"_y_110.off",
	"_y_115.off",
	"_y_120.off",
	"_y_125.off",
	"_y_130.off",
	"_y_135.off",
	"_y_140.off",
	"_y_145.off",
	"_y_150.off",
	"_y_155.off",
	"_y_160.off",
	"_y_165.off",
	"_y_170.off",
	"_y_175.off"
};

const char* out_z_filenames[35] =
{
	"_z_5.off",
	"_z_10.off",
	"_z_15.off",
	"_z_20.off",
	"_z_25.off",
	"_z_30.off",
	"_z_35.off",
	"_z_40.off",
	"_z_45.off",
	"_z_50.off",
	"_z_55.off",
	"_z_60.off",
	"_z_65.off",
	"_z_70.off",
	"_z_75.off",
	"_z_80.off",
	"_z_85.off",
	"_z_90.off",
	"_z_95.off",
	"_z_100.off",
	"_z_105.off",
	"_z_110.off",
	"_z_115.off",
	"_z_120.off",
	"_z_125.off",
	"_z_130.off",
	"_z_135.off",
	"_z_140.off",
	"_z_145.off",
	"_z_150.off",
	"_z_155.off",
	"_z_160.off",
	"_z_165.off",
	"_z_170.off",
	"_z_175.off"
};

void SetRotMatrix(double rx[3][3], double ry[3][3], double rz[3][3], double theta_x, double theta_y, double theta_z)
{
  rx[0][0] = 1;
  rx[1][1] = cos(theta_x);
  rx[1][2] = sin(theta_x);
  rx[2][2] = cos(theta_x);
  rx[2][1] = -sin(theta_x);

  ry[1][1] = 1;
  ry[0][0] = cos(theta_y);
  ry[0][2] = -sin(theta_y);
  ry[2][2] = cos(theta_y);
  ry[2][0] = sin(theta_y);

  rz[2][2] = 1;
  rz[0][0] = cos(theta_z);
  rz[0][1] = sin(theta_z);
  rz[1][1] = cos(theta_z);
  rz[1][0] = -sin(theta_z);
}
const double PI = 3.1415926536;
int main(int argc, char* argv[])
{
	int num = 10;
	if (argc == 3)
	{
		num = atoi(argv[2]);
	}
	else if (argc != 2)
	{
		std::cout << "Usage: gen_rot_meshes <triangle_off_file> <num>" << std::endl;
		return -1;
	}
	MeshFileReader vtk(argv[1]);
	const Mesh& mesh = vtk.GetMesh();
	const std::vector<Vertex>& V = mesh.V;
	const std::vector<Cell>& C = mesh.C;

	double rx[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	double ry[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	double rz[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	for (int i = 0; i < num; i++)
	{
		float theta_x = 2.0 * PI * i / num;
		for (int j = 0; j < num; j++)
		{
			float theta_y = 2.0 * PI * j/ num;
			for (int k = 0; k < num; k++)
			{
				float theta_z = 2.0 * PI * k/ num;
				SetRotMatrix(rx, ry, rz, theta_x, theta_y, theta_z);
				std::vector<Vertex> Vx(V);
				for (size_t l = 0; l < Vx.size(); l++)
				{
				    Vx[l].x = rx[0][0]*V[l].x + rx[0][1]*V[l].y + rx[0][2]*V[l].z;
				    Vx[l].y = rx[1][0]*V[l].x + rx[1][1]*V[l].y + rx[1][2]*V[l].z;
				    Vx[l].z = rx[2][0]*V[l].x + rx[2][1]*V[l].y + rx[2][2]*V[l].z;

					Vx[l].x = ry[0][0]*V[l].x + ry[0][1]*V[l].y + ry[0][2]*V[l].z;
					Vx[l].y = ry[1][0]*V[l].x + ry[1][1]*V[l].y + ry[1][2]*V[l].z;
					Vx[l].z = ry[2][0]*V[l].x + ry[2][1]*V[l].y + ry[2][2]*V[l].z;

					Vx[l].x = rz[0][0]*V[l].x + rz[0][1]*V[l].y + rz[0][2]*V[l].z;
					Vx[l].y = rz[1][0]*V[l].x + rz[1][1]*V[l].y + rz[1][2]*V[l].z;
					Vx[l].z = rz[2][0]*V[l].x + rz[2][1]*V[l].y + rz[2][2]*V[l].z;
				}
				std::string strx(32, 0);// = itoa(int(360*i/num));
				std::string stry(32, 0);// = itoa(int(360*j/num));
				std::string strz(32, 0);// = itoa(int(360*k/num));
				sprintf((char*)strx.data(), "%d",  360*i/num);
				sprintf((char*)stry.data(), "%d",  360*j/num);
				sprintf((char*)strz.data(), "%d",  360*k/num);
				const std::string str = std::string(argv[1]).substr(0, strlen(argv[1]) - 4) + std::string("_") + std::string(strx.c_str()) + std::string("_") + std::string(stry.c_str()) + std::string("_") + std::string(strz.c_str()) + std::string(".off");
				MeshFileWriter meshWriter(Vx, C, str.c_str());
				meshWriter.WriteMeshFile();

			}
		}
	}
	/*
	for (int i = 0; i < 35; i++)
	{
		std::vector<Vertex> Vx(V);
		for (size_t j = 0; j < Vx.size(); j++)
		{
			Vx[j].x = X_RM[i][0][0]*V[j].x + X_RM[i][0][1]*V[j].y + X_RM[i][0][2]*V[j].z;
			Vx[j].y = X_RM[i][1][0]*V[j].x + X_RM[i][1][1]*V[j].y + X_RM[i][1][2]*V[j].z;
			Vx[j].z = X_RM[i][2][0]*V[j].x + X_RM[i][2][1]*V[j].y + X_RM[i][2][2]*V[j].z;
		}
		const std::string str = std::string(argv[1]).substr(0, strlen(argv[1]) - 4) + std::string(out_x_filenames[i]);
		MeshFileWriter meshWriter(Vx, C, str.c_str(), MESH_TYPE_TRIANGLE_OFF);
		meshWriter.WriteMeshFile(str.c_str(), MESH_TYPE_TRIANGLE_OFF);
	}

	for (int i = 0; i < 35; i++)
	{
		std::vector<Vertex> Vy(V);
		for (size_t j = 0; j < Vy.size(); j++)
		{
			Vy[j].x = Y_RM[i][0][0]*V[j].x + Y_RM[i][0][1]*V[j].y + Y_RM[i][0][2]*V[j].z;
			Vy[j].y = Y_RM[i][1][0]*V[j].x + Y_RM[i][1][1]*V[j].y + Y_RM[i][1][2]*V[j].z;
			Vy[j].z = Y_RM[i][2][0]*V[j].x + Y_RM[i][2][1]*V[j].y + Y_RM[i][2][2]*V[j].z;
		}
		const std::string str = std::string(argv[1]).substr(0, strlen(argv[1]) - 4) + std::string(out_y_filenames[i]);
		MeshFileWriter meshWriter(Vy, C, str.c_str(), MESH_TYPE_TRIANGLE_OFF);
		meshWriter.WriteMeshFile(str.c_str(), MESH_TYPE_TRIANGLE_OFF);
	}

	for (int i = 0; i < 35; i++)
	{
		std::vector<Vertex> Vz(V);
		for (size_t j = 0; j < Vz.size(); j++)
		{
			Vz[j].x = Z_RM[i][0][0]*V[j].x + Z_RM[i][0][1]*V[j].y + Z_RM[i][0][2]*V[j].z;
			Vz[j].y = Z_RM[i][1][0]*V[j].x + Z_RM[i][1][1]*V[j].y + Z_RM[i][1][2]*V[j].z;
			Vz[j].z = Z_RM[i][2][0]*V[j].x + Z_RM[i][2][1]*V[j].y + Z_RM[i][2][2]*V[j].z;
		}
		const std::string str = std::string(argv[1]).substr(0, strlen(argv[1]) - 4) + std::string(out_z_filenames[i]);
		MeshFileWriter meshWriter(Vz, C, str.c_str(), MESH_TYPE_TRIANGLE_OFF);
		meshWriter.WriteMeshFile(str.c_str(), MESH_TYPE_TRIANGLE_OFF);
	}
	*/
	return 0;
}
