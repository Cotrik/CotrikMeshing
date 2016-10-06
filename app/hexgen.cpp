/*
 * hexgen.cpp
 *
 *  Created on: Jun 3, 2015
 *      Author: cotrik
 */

#include "include.h"
std::vector<bool> paras_flag;

void GetParameter(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p);
void MapbackToOrigTet(const Mesh& tetMesh, const Mesh& hexMesh, const TetHexParas& p, std::vector<Vertex>& hexV);
int main(int argc, char* argv[])
{
	std::string strInputOrigTetFileName;
	std::string strInputPolyTetFileName;
	std::string strInputPolyHexFileName;
	std::string strOutputHexFileName;

	try
	{
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message.")
			("thread,j", po::value<unsigned int>(&THREAD_NUM)->default_value(4), "threads number")
			// ########## input #########
			// -i
			("input-tet-file,i", po::value<std::string>(&strInputOrigTetFileName)->default_value("orig.tet.vtk"),
				"name of the orig input tetrahedral vtk file")
			// -t
			("polycube-tet-file,t", po::value<std::string>(&strInputPolyTetFileName)->default_value("polycube.tet.vtk"),
				"name of the polycube input tetrahedral vtk file")
			// -p
			("polycube-hex-file,p", po::value<std::string>(&strInputPolyHexFileName)->default_value("polycube.hex.vtk"),
				"name of the polycube input hexahedral vtk file")
			// ########## output #########
			// -o
			("output-hex-file,o", po::value<std::string>(&strOutputHexFileName)->default_value("hex.vtk"),
				"name of the output hexahedral vtk file")
			;

		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);

		if (vm.count("help"))
		{
			std::cout << desc << "\n";
			return 1;
		}

		if (vm.count("input-tet-file"))		std::cout << "strInputOrigTetFileName = " << strInputOrigTetFileName << std::endl;
		if (vm.count("polycube-tet-file"))		std::cout << "strInputPolyTetFileName = " << strInputPolyTetFileName << std::endl;
		if (vm.count("polycube-hex-file"))		std::cout << "strInputPolyHexFileName = " << strInputPolyHexFileName << std::endl;
		if (vm.count("output-hex-file"))		std::cout << "strOutputHexFileName = " << strOutputHexFileName << std::endl;
#else
		strInputOrigTetFileName = "orig.tet.vtk";
		strInputPolyTetFileName = "polycube.tet.vtk";
		strInputPolyHexFileName = "polycube.hex.vtk";
		strOutputHexFileName = "hex.vtk";
#endif
	} catch (std::exception& e)
	{
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	} catch (...)
	{
		std::cerr << "Exception of unknown type!\n";
	}

	MeshFileReader orgTetMeshFileReader(strInputOrigTetFileName.c_str());
	Mesh origTetMesh(orgTetMeshFileReader.GetMesh());

	MeshFileReader polycubeTetMeshFileReader(strInputPolyTetFileName.c_str());
	PolycubeMesh polycubeTetMesh(polycubeTetMeshFileReader.GetMesh());
	////////
	polycubeTetMesh.GetMaxMinCoordinates();
	polycubeTetMesh.ExtractSurface();
	polycubeTetMesh.GetFaceType();
	polycubeTetMesh.GetFaceAndNeighborFaces();
	polycubeTetMesh.GetFacePatches_N();
	polycubeTetMesh.GetNormalOfSurfaceFaces();
	polycubeTetMesh.GetNormalOfSurfaceVertices();
	polycubeTetMesh.GetFacePatches();
	polycubeTetMesh.GetEdgePatches();
	polycubeTetMesh.GetVertexPatches();
	//////////////////////////////
	polycubeTetMesh.GetPatches();
	polycubeTetMesh.SortPatches();
	polycubeTetMesh.ModifyPatchesPosition();

	polycubeTetMesh.GetMaxMinCoordinates_Patch();
	polycubeTetMesh.AdjustPolycubeSurfaceToHoldAllSmallCubesVolume1(0.02f, polycubeTetMesh.m_minVertex);
	polycubeTetMesh.AddCompensation();

	polycubeTetMesh.GetPatches();
	polycubeTetMesh.SortPatches();
	polycubeTetMesh.ModifyPatchesPosition();

	MeshFileWriter polycubeTetMeshFileWriter(polycubeTetMesh.V, polycubeTetMesh.C, "polycube.tet.n.a.vtk", TETRAHEDRA);
	polycubeTetMeshFileWriter.WriteFile();
	////////

	MeshFileReader polycubeHexMeshFileReader(strInputPolyHexFileName.c_str());
	PolycubeMesh polycubeHexMesh(polycubeHexMeshFileReader.GetMesh());

	TetHexParas p;
	GetParameter(polycubeTetMesh, polycubeHexMesh, p);
	// GetMagnifiedParas(polycubeTetMesh, polycubeHexMesh, p);

	std::vector<Vertex> hexV;
	MapbackToOrigTet(origTetMesh, polycubeHexMesh, p, hexV);

	MeshFileWriter origHexMeshFileWriter(hexV, polycubeHexMesh.C, strOutputHexFileName.c_str(), HEXAHEDRA);
	origHexMeshFileWriter.WriteFile();

	return 0;
}

void GetParameter(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p)
{
    p.paras.resize(hexMesh.V.size());
    p.tetIndex.resize(hexMesh.V.size());
    p.flag.resize(hexMesh.V.size());
    for (int i = 0; i < hexMesh.V.size(); i++)
    {
        bool flag = false;
        const Vertex& v = hexMesh.V.at(i);
        for (int j = 0; j < tetMesh.C.size(); j++)
        {
            const Cell& tet = tetMesh.C.at(j);
            const Vertex& p0 = tetMesh.V[tet.at(0)];
            const Vertex& p1 = tetMesh.V[tet.at(1)];
            const Vertex& p2 = tetMesh.V[tet.at(2)];
            const Vertex& p3 = tetMesh.V[tet.at(3)];

            glm::vec3 v03((p0.x - p3.x), (p0.y - p3.y), (p0.z - p3.z));
            glm::vec3 v13((p1.x - p3.x), (p1.y - p3.y), (p1.z - p3.z));
            glm::vec3 v23((p2.x - p3.x), (p2.y - p3.y), (p2.z - p3.z));

            glm::vec3 r(v.x, v.y, v.z);
            glm::vec3 r4(p3.x, p3.y, p3.z);

            glm::vec3 rr4 = r - r4;

            glm::mat3x3 T(v03, v13, v23);
            glm::mat3x3 T_inverse = glm::inverse(T);
            glm::vec3 rambda = T_inverse * (r - r4);

            if (rambda.x > -1e-3 && rambda.y > -1e-3 && rambda.z>-1e-3 && (rambda.x + rambda.y + rambda.z) < 1.001){
                p.paras.at(i) = rambda;
                p.tetIndex.at(i) = j;
                p.flag.at(i) = true;
                flag = true;
                break;
            }
        }
        if (!flag)
        {
            ///////////////////////
            const glm::vec3 vv(v.x, v.y, v.z);
            double closeCellIndex = 0;
            double dis = 100000000;
            for (int j = 0; j < tetMesh.C.size(); j++)
            {
                const Cell& tet = tetMesh.C.at(j);
                const Vertex& p0 = tetMesh.V[tet.at(0)];
                const Vertex& p1 = tetMesh.V[tet.at(1)];
                const Vertex& p2 = tetMesh.V[tet.at(2)];
                const Vertex& p3 = tetMesh.V[tet.at(3)];

                const glm::vec3 v0(p0.x, p0.y, p0.z);
                const glm::vec3 v1(p1.x, p1.y, p1.z);
                const glm::vec3 v2(p2.x, p2.y, p2.z);
                const glm::vec3 v3(p3.x, p3.y, p3.z);

                const glm::vec3 d0(v0 - vv);
                const glm::vec3 d1(v1 - vv);
                const glm::vec3 d2(v2 - vv);
                const glm::vec3 d3(v3 - vv);

                double d = glm::length(d0) + glm::length(d1) + glm::length(d2) + glm::length(d3);
                if (d < dis)
                {
                    dis = d;
                    closeCellIndex = j;
                }
            }
            const Cell& tet = tetMesh.C.at(closeCellIndex);
            const Vertex& p0 = tetMesh.V[tet.at(0)];
            const Vertex& p1 = tetMesh.V[tet.at(1)];
            const Vertex& p2 = tetMesh.V[tet.at(2)];
            const Vertex& p3 = tetMesh.V[tet.at(3)];

            const glm::vec3 v0(p0.x, p0.y, p0.z);
            const glm::vec3 v1(p1.x, p1.y, p1.z);
            const glm::vec3 v2(p2.x, p2.y, p2.z);
            const glm::vec3 v3(p3.x, p3.y, p3.z);

            const glm::vec3 d0(v0 - vv);
            const glm::vec3 d1(v1 - vv);
            const glm::vec3 d2(v2 - vv);
            const glm::vec3 d3(v3 - vv);
            const double d = glm::length(d0) + glm::length(d1) + glm::length(d2) + glm::length(d3);
            //rambda.x = d0/d; rambda.y = d1/d; rambda.z = d2/d;
            ///////////////////////
            //glm::vec3 rambda(0.25, 0.25, 0.5);
            glm::vec3 rambda(glm::length(d0)/d, glm::length(d1)/d, glm::length(d2)/d);
            p.paras.at(i) = rambda;
            p.tetIndex.at(i) = closeCellIndex;
            p.flag.at(i) = false;

            //errorPointIndices.push_back(i);
            std::cout << "i = " << i << " Fail to parametrization!" << std::endl;
            paras_flag.push_back(false);
        }
    }
}

void MapbackToOrigTet(const Mesh& tetMesh, const Mesh& hexMesh, const TetHexParas& p, std::vector<Vertex>& hexV)
{
    //parameterization of new OrigTet to new Cube Mesh
    std::cout << "MapbackToOrigTet" << std::endl;
    hexV.resize(hexMesh.V.size());
    for (int i = 0; i < hexMesh.V.size(); i++)
    {
        const unsigned long tetIndex = p.tetIndex.at(i);
        const Cell& tet = tetMesh.C.at(tetIndex);

        const Vertex& origT0 = tetMesh.V[tet.at(0)];
        const Vertex& origT1 = tetMesh.V[tet.at(1)];
        const Vertex& origT2 = tetMesh.V[tet.at(2)];
        const Vertex& origT3 = tetMesh.V[tet.at(3)];

        glm::vec3 p03((origT0.x - origT3.x), (origT0.y - origT3.y), (origT0.z - origT3.z));
        glm::vec3 p13((origT1.x - origT3.x), (origT1.y - origT3.y), (origT1.z - origT3.z));
        glm::vec3 p23((origT2.x - origT3.x), (origT2.y - origT3.y), (origT2.z - origT3.z));
        glm::mat3x3 OrigT(p03, p13, p23);
        glm::vec3 new_r4(origT3.x, origT3.y, origT3.z);

        const glm::vec3& rambda = p.paras.at(i);
        glm::vec3 new_r = OrigT * rambda;
        new_r += new_r4;
        Vertex baryCenter(new_r.x, new_r.y, new_r.z);

        Vertex& v = hexV.at(i);
        v.x = baryCenter.x; v.y = baryCenter.y; v.z = baryCenter.z;
    }
}
