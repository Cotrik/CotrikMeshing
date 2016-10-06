/*
 * hexmapback.cpp
 *
 *  Created on: May 14, 2015
 *      Author: cotrik
 */

#include "include.h"

std::vector<Cell> tets;
std::vector<glm::vec3> paras;
std::vector<bool> paras_flag;

void Project(Mesh& hexMesh, const Mesh& tetMesh)
{
    hexMesh.GetSurfaceVertexIndices();
    for (int i = 0; i < hexMesh.surfaceVertexIndices.size(); i++)
    {
        unsigned long index = hexMesh.surfaceVertexIndices.at(i);
        Vertex& v = hexMesh.V.at(index);
        glm::vec3 orig(v.x, v.y, v.z);
        glm::vec3 dir(hexMesh.N_V.at(index));
        glm::vec3 dir2(-dir.x, -dir.y, -dir.z);
        std::vector<Position_Distance> vec_p_d;
        for (int j = 0; j < tetMesh.surface.size(); j++)
        {
            const Face& face = tetMesh.surface.at(j);
            const Vertex& v_0 = tetMesh.V.at(face.at(0));
            const Vertex& v_1 = tetMesh.V.at(face.at(1));
            const Vertex& v_2 = tetMesh.V.at(face.at(2));

            glm::vec3 v0(v_0.x, v_0.y, v_0.z);
            glm::vec3 v1(v_1.x, v_1.y, v_1.z);
            glm::vec3 v2(v_2.x, v_2.y, v_2.z);

            glm::vec3 n = glm::cross(v1 - v0, v2 - v0);
            double d = -glm::dot(n, v0);
            double length = 1.0 / n.length();
            n *= length;
            d *= length;

            double a = n.x;
            double b = n.y;
            double c = n.z;

            double t = (glm::dot(orig, n) + d)/(a*a + b*b + c*c);

            glm::vec3 p(orig.x - a*t, orig.y - b*t, orig.z - c*t);
            if (GeoUtil::IsPointInTriangle(p, v0, v1, v2, i, j))
            {
                Position_Distance pd;
                pd.position = p;
                pd.distance = glm::length(orig - p);
                vec_p_d.push_back(pd);
            }
        }
        if (!vec_p_d.empty())
        {
            std::sort(vec_p_d.begin(), vec_p_d.end());
            if (vec_p_d.at(0).distance < 0.01)
            {
            v.x = vec_p_d.at(0).position.x;
            v.y = vec_p_d.at(0).position.y;
            v.z = vec_p_d.at(0).position.z;
            }
            else
            {
                std::cout << index << " Fail projection, distance too big : " <<  vec_p_d.at(0).distance << std::endl;
            }
        }
        else
        {
            std::cout << index << " Fail projection" << std::endl;
        }
    }
}


static void GetMagnifiedParas2(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p)
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

//                if (i == 2775 || i == 6776 || i == 1875)
//                {
//                    if(fabs(rambda.x) + fabs(rambda.y) + fabs(rambda.z) < 1.3)
//                    std::cout << "i = " << i << " " << "j = " << j << " " << rambda.x << ", " << rambda.y << ", " << rambda.z  << ", " << std::endl;
//                }

                if (rambda.x > -5*1e-2 && rambda.y > -5*1e-2 && rambda.z>-5*1e-2 && (rambda.x + rambda.y + rambda.z) < 1.1){
                p.paras.at(i) = rambda;
                p.tetIndex.at(i) = j;
                p.flag.at(i) = true;
                flag = true;
                break;
                }

        }
        if (!flag)
        {
            glm::vec3 rambda(0.25, 0.25, 0.5);

            p.paras.at(i) = rambda;
            p.tetIndex.at(i) = 0;
            p.flag.at(i) = false;
            if (i == 2775)
            {
                p.paras.at(i) = glm::vec3(-0.0405468, 0.93677, 0.271008);
                p.tetIndex.at(i) = 10315;
            }
            if (i == 6776)
            {
                p.paras.at(i) = glm::vec3(0.397704, 0.0734921, -0.499427);
                p.tetIndex.at(i) = 3161;
            }

 //           errorPointIndices.push_back(i);
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
//      if (p.flag.at(i)){
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
        //std::cout << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << " " << rambda.x << " " << rambda.y << " " << rambda.z << std::endl;
//      }
//      else
//      {
//          hexV.at(i) = hexMesh_u.V.at(i);
//      }
    }
}

int HexMapBack(int argc, char* argv[])
{
    std::string strInputOrigTetFileName;
    std::string strInputMagnifiedTetFileName;
    std::string strInputMagnifiedHexFileName;
//  std::string strInputMagnifiedHexFileName_u;
    std::string strOutputOrigHexFileName;
    try
    {
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help,h", "produce help message.")
            // ########## input #########
            // -i
            ("orig-tet,i", po::value<std::string>(&strInputOrigTetFileName)->default_value("orig.tet.vtk"), "name of the orig input tetrahedral vtk file")
            // -u
            ("magnified-tet,u", po::value<std::string>(&strInputMagnifiedTetFileName)->default_value("magnified.tet.vtk"), "name of the orig input tetrahedral vtk file")
            // -v
            ("magnified-hex,v", po::value<std::string>(&strInputMagnifiedHexFileName)->default_value("magnified.hex_optimized.vtk"), "name of the magnified input hexahedral vtk file")
            // -w
//          ("magnified-hex-u,w", po::value<std::string>(&strInputMagnifiedHexFileName_u)->default_value("magnified.hex.vtk"), "name of the magnified input hexahedral vtk file")
            // ########## output #########
            // -o
            ("orig-hex,o", po::value<std::string>(&strOutputOrigHexFileName)->default_value("orig.hex.vtk"),    "name of the orig output hexahedral vtk file")
            ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        if (vm.count("help"))       {           std::cout << desc << "\n";          return 1;       }
        if (vm.count("orig-tet"))       std::cout << "strInputOrigTetFileName = "      << strInputOrigTetFileName      << std::endl;
        if (vm.count("magnified-tet"))  std::cout << "strInputMagnifiedTetFileName = " << strInputMagnifiedTetFileName << std::endl;
        if (vm.count("magnified-hex"))  std::cout << "strInputMagnifiedHexFileName = " << strInputMagnifiedHexFileName << std::endl;
//      if (vm.count("magnified-hex-u"))    std::cout << "strInputMagnifiedHexFileName_u = " << strInputMagnifiedHexFileName_u << std::endl;
        if (vm.count("orig-hex"))       std::cout << "strOutputOrigHexFileName = "     << strOutputOrigHexFileName     << std::endl;
#else
        strInputOrigTetFileName      = "orig.tet.vtk";
        strInputMagnifiedTetFileName = "magnified.tet.vtk";
        strInputMagnifiedHexFileName = "magnified.hex.vtk";
        strOutputOrigHexFileName     = "orig.hex.vtk";
#endif
    } catch (std::exception& e)
    {       std::cerr << "error: " << e.what() << "\n";     return 1;
    } catch (...)   {
        std::cerr << "Exception of unknown type!\n";
    }

    //////////////////////////////////////////
    // Read Input Tet Mesh
    //////////////////////////////////////////
    MeshFileReader origTetMeshFileReader(strInputOrigTetFileName.c_str());
    Mesh origTetMesh(origTetMeshFileReader.GetMesh());

    MeshFileReader magnifiedTetMeshFileReader(strInputMagnifiedTetFileName.c_str());
    Mesh magnifiedTetMesh(magnifiedTetMeshFileReader.GetMesh());

    MeshFileReader magnifiedHexMeshFileReader(strInputMagnifiedHexFileName.c_str());
    Mesh magnifiedHexMesh(magnifiedHexMeshFileReader.GetMesh());

    MeshFileWriter origHexMeshFileWriter2(magnifiedHexMesh.V, magnifiedHexMesh.C, "optimized.hex.vtk");
    origHexMeshFileWriter2.WriteMeshFile();

    magnifiedTetMesh.ExtractSurface();
    magnifiedHexMesh.ExtractSurface();
    magnifiedHexMesh.GetNormalOfSurfaceFaces();
    magnifiedHexMesh.GetNormalOfSurfaceVertices();
    Project(magnifiedHexMesh, magnifiedTetMesh);

    MeshFileWriter origHexMeshFileWriter1(magnifiedHexMesh.V, magnifiedHexMesh.C, "project.hex.vtk");
    origHexMeshFileWriter1.WriteMeshFile();

    TetHexParas p;
    GetMagnifiedParas2(magnifiedTetMesh, magnifiedHexMesh, p);

    std::vector<Vertex> hexV;
    MapbackToOrigTet(origTetMesh, magnifiedHexMesh/*, magnifiedHexMesh_u*/, p, hexV);

    MeshFileWriter origHexMeshFileWriter(hexV, magnifiedHexMesh.C, strOutputOrigHexFileName.c_str());
    origHexMeshFileWriter.WriteMeshFile();

    MeshFileWriter origHexMeshFileWriter3(hexV, magnifiedHexMesh.C, "mapback.hex.mesh");
    origHexMeshFileWriter3.WriteMeshFile();

    Mesh hexMesh(hexV, magnifiedHexMesh.C);

    //errorPointIndices.push_back(0);
    ///////////////////////////////////////////////
//    pTetMesh = (Mesh*) &origTetMesh;
//    pTetMagnifiedMesh = (Mesh*) &magnifiedTetMesh;
//    pHexMagnifiedMesh = (Mesh*) &magnifiedHexMesh;
//    pHexMesh = (Mesh*) &hexMesh;
//    DisplayMesh(argc, argv, magnifiedHexMesh);

    return 0;
}


int main(int argc, char* argv[])
{
	return HexMapBack(argc, argv);
}


