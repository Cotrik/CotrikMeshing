/*
 * check_overlap.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: cotrik
 */

#include "include.h"

std::vector<bool> paras_flag;
std::vector<unsigned long> degenerated;
void GetParameter(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p);
void GetParameter1(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p);
void MapbackToOrigTet(const Mesh& tetMesh, const Mesh& hexMesh, const TetHexParas& p, std::vector<Vertex>& hexV);
void CheckDegeneracy(const Mesh& tetMesh);

int main(int argc, char* argv[])
{
    if (argc != 3)
    {
        std::cout << "Usage: check_overlap <input_tet_mesh> <input_hex_mesh>" << std::endl;
        return -1;
    }
    MeshFileReader tetReader(argv[1]);
    MeshFileReader hexReader(argv[2]);
//    Mesh tetMesh(tetReader.GetMesh());
//    Mesh hexMesh(hexReader.GetMesh());
    PolycubeMesh polycubeTetMesh(tetReader.GetMesh());
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

    CheckDegeneracy(polycubeTetMesh);
    TetHexParas p;
    GetParameter(polycubeTetMesh, hexReader.GetMesh(), p);

    return 0;
}

void CheckDegeneracy(const Mesh& tetMesh)
{
    for (int j = 0; j < tetMesh.C.size(); j++)
    {
        const Cell& tet = tetMesh.C.at(j);
        if (GeoUtil::IsTetrahedronDegenerated(tet, tetMesh.V))
        {
            std::cout << "Tet " << j << " is Degenerated\n";
            degenerated.push_back(j);
        }
    }
    degenerated.resize(degenerated.size());
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
        int count = 0;
        std::vector<glm::vec3> l;
        std::vector<int> t;
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

//            if (rambda.x >= 0 && rambda.y >= 0 && rambda.z>=0 && (rambda.x + rambda.y + rambda.z) <= 1.00000){
            if (rambda.x > -1e-3 && rambda.y > -1e-3 && rambda.z>-1e-3 && (rambda.x + rambda.y + rambda.z) < 1.001){
                p.paras.at(i) = rambda;
                p.tetIndex.at(i) = j;
                p.flag.at(i) = true;
                flag = true;
                l.push_back(rambda);
                t.push_back(j);
                //break;
//                if (1 - rambda.x < 1e-6 || 1 - rambda.y < 1e-6 || 1 - rambda.z < 1e-6 || (rambda.x + rambda.y + rambda.z) < 1e-6);
//                else
//                    if (std::find(degenerated.begin(), degenerated.end(), j) != degenerated.end())
//                        count++;
                count++;
            }
        }
        if (count > 1)
        {
            std::cout << "polycube.hex.vtk vertex " << i << " overlaps, count = " << count
                    << " (" << p.paras.at(i).x << ", " << p.paras.at(i).y << ", " << p.paras.at(i).z << ")" << std::endl;
//            for (int k = 0; k < l.size(); k++)
//                std::cout << " (" << l.at(k).x << ", " << l.at(k).y << ", " << l.at(k).z << ")" << std::endl;
            std::cout << " ( " ;
            for (int k = 0; k < l.size(); k++)
                std::cout << t.at(k) << " " ;
            std::cout << ")" << std::endl;
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

void GetParameter1(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p)
{
    p.paras.resize(hexMesh.V.size());
    p.tetIndex.resize(hexMesh.V.size());
    p.flag.resize(hexMesh.V.size());
    glm::vec3 lambda;
    for (int i = 0; i < hexMesh.V.size(); i++)
    {
        bool flag = false;
        const Vertex& v = hexMesh.V.at(i);
        int count = 0;
        for (int j = 0; j < tetMesh.C.size(); j++)
        {
            const Cell& tet = tetMesh.C.at(j);
            if (GeoUtil::IsVertexInsideTetrahedron(v, tet, tetMesh.V, lambda))
            {
                p.paras.at(i) = lambda;
                flag = true;
                count++;
            }
        }
        if (count > 1)
        {
            std::cout << "polycube.hex.vtk vertex " << i << " overlaps, count = " << count
                    << " (" << lambda.x << ", " << lambda.y << ", " << lambda.z << ")" << std::endl;
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

