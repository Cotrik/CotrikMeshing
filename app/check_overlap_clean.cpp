/*
 * check_overlap.cpp
 *
 *  Created on: Jul 27, 2016
 *      Author: cotrik
 */

#include "include.h"

std::vector<bool> paras_flag;
std::vector<unsigned long> degenerated;
struct HexV_TetC_Paras
{
    unsigned long hexVId;
    std::vector<unsigned long> tetCIds;
    std::vector<glm::vec3> paras;
    std::vector<glm::vec3> locs;
    std::vector<unsigned long> closestSufaceV;
    std::vector<double> geoDis;
    double min_distance_tetCId;
    double max_distance_tetCId;
    unsigned long min_distance_paras_index;
    unsigned long max_distance_paras_index;
};

//std::vector<HexV_TetC_Paras> hexV_tetC_paras;
std::vector<unsigned long> neigboringHexV;

enum FarNear
{
    FAR,
    NEAR
};

std::vector<FarNear> neigboring_FarNear;
////////////////////////////////////////////////
// GeodesicDistance members
vtkPolyDataReader* reader;
vtkPolyDataNormals* normals;
vtkPolyData* pd;
vtkPolyDataMapper* mapper;
vtkActor* actor;
vtkRenderer* ren;
vtkRenderWindowInteractor* iren;
vtkRenderWindow* renWin;
vtkContourWidget* contourWidget;
vtkPolygonalSurfacePointPlacer* pointPlacer;
vtkPolygonalSurfaceContourLineInterpolator2* interpolator;
////////////////////////////////////////////////
void GetParameter(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p,
        std::vector<unsigned long>& overlapV, std::vector<HexV_TetC_Paras>& hexV_tetC_paras);
void GetParameter1(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p);
void MapbackToOrigTet(const Mesh& tetMesh, const Mesh& hexMesh, const TetHexParas& p,
        std::vector<Vertex>& hexV);
void CheckDegeneracy(const Mesh& tetMesh);
void CheckOverlapCells(const Mesh& hexMesh, const std::vector<unsigned long>& overlapV,
        std::vector<unsigned long>& overlapC);

std::vector<unsigned long> hexvlist;
void GetOverlapHtp(const Mesh& hexMesh, const std::vector<unsigned long>& overlapC,
        std::vector<HexV_TetC_Paras>& hexV_tetC_paras)
{
//    std::vector<unsigned long> hexvlist;

    for (size_t i = 0; i < overlapC.size(); i++)
    {
        const Cell& hex = hexMesh.C.at(overlapC.at(i));
        for (size_t j = 0; j < hex.size(); j++)
            hexvlist.push_back(hex.at(j));
    }

    std::sort(hexvlist.begin(), hexvlist.end());
    std::vector<unsigned long>::iterator iter = std::unique(hexvlist.begin(), hexvlist.end());
    hexvlist.resize(std::distance(hexvlist.begin(), iter));

    ///////////////////
    std::cout << "--------------- hexvlist --------------\n";
    for (size_t i = 0; i < hexvlist.size(); i++)
        std::cout << hexvlist[i] << " ";
    std::cout << endl;

    for (std::vector<HexV_TetC_Paras>::iterator iterhtp = hexV_tetC_paras.begin(); iterhtp != hexV_tetC_paras.end();)
    {
        const unsigned long hexVId = iterhtp->hexVId;
        if (std::find(hexvlist.begin(), hexvlist.end(), hexVId) != hexvlist.end())
        {
            ++iterhtp;
        }
        else
        {
            iterhtp = hexV_tetC_paras.erase(iterhtp);
        }
    }
}
void GetOrigSpaceLocation(const Mesh& tetMesh, const unsigned long tetIndex, const glm::vec3& lambda, glm::vec3& new_r)
{
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

    new_r = OrigT * lambda;
    new_r += new_r4;
}

unsigned long GetClosestSurfacePoint(const Mesh& tetMesh, const glm::vec3& p)
{
    double min_distance = 1000000000000;
    unsigned long vId = 0;
    for (unsigned long i = 0; i < tetMesh.V.size(); i++)
    {
        if (!tetMesh.V.at(i).vinfo.bSurface)
            continue;
        Vertex v(p);
        const double distance = tetMesh.V.at(i).Distance(v);
        if (distance < min_distance)
        {
            vId = i;
            min_distance = distance;
        }
    }

    return vId;
}

void InitGeodesicDistance(const char* polydataFilename)
{
    reader = vtkPolyDataReader::New();
    reader->SetFileName(polydataFilename);

    normals = vtkPolyDataNormals::New();

    const int geodesicMethod = 0;
    const int interpolationOrder = 0;
    const double distanceOffset = 0;

    // We need to ensure that the dataset has normals if a distance offset was
    // specified.
    if (fabs(distanceOffset) > 1e-6)
    {
        normals->SetInputConnection(reader->GetOutputPort());
        normals->SplittingOff();

        // vtkPolygonalSurfacePointPlacer needs cell normals
        // vtkPolygonalSurfaceContourLineInterpolator needs vertex normals
        normals->ComputeCellNormalsOn();
        normals->ComputePointNormalsOn();
        normals->Update();
    }

    vtkPolyData *pd = (fabs(distanceOffset) > 1e-6) ? normals->GetOutput() : reader->GetOutput();

    mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection(fabs(distanceOffset) > 1e-6 ? normals->GetOutputPort() : reader->GetOutputPort());

    actor = vtkActor::New();
    actor->SetMapper(mapper);

    ren = vtkRenderer::New();
    renWin = vtkRenderWindow::New();
    renWin->AddRenderer(ren);
    iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);

    ren->AddActor(actor);
    ren->GetActiveCamera()->SetPosition(-3.68, .447, 1.676);
    ren->GetActiveCamera()->Roll(150);
    ren->ResetCamera();
    ren->ResetCameraClippingRange();

    ren->AddActor(actor);

    // Create a contour widget to interactively trace the path

    contourWidget = vtkContourWidget::New();
    contourWidget->SetInteractor(iren);
    vtkOrientedGlyphContourRepresentation *rep = vtkOrientedGlyphContourRepresentation::SafeDownCast(contourWidget->GetRepresentation());
    rep->GetLinesProperty()->SetColor(1, 0.2, 0);
    rep->GetLinesProperty()->SetLineWidth(5.0);

    pointPlacer = vtkPolygonalSurfacePointPlacer::New();
    pointPlacer->AddProp(actor);
    pointPlacer->GetPolys()->AddItem(pd);
    rep->SetPointPlacer(pointPlacer);

    // Snap the contour nodes to the closest vertices on the mesh
    pointPlacer->SnapToClosestPointOn();

    interpolator = vtkPolygonalSurfaceContourLineInterpolator2::New();
    interpolator->GetPolys()->AddItem(pd);
    interpolator->SetGeodesicMethod(geodesicMethod);
    interpolator->SetInterpolationOrder(interpolationOrder);
    rep->SetLineInterpolator(interpolator);
    if (fabs(distanceOffset) > 1e-6)
    {
        pointPlacer->SetDistanceOffset(distanceOffset);
        interpolator->SetDistanceOffset(distanceOffset);
    }

//    renWin->Render();
//    iren->Initialize();
    contourWidget->EnabledOn();
    //interpolator->InterpolateLine(reader->GetOutput(), 14483, 14574);
}

double ComputeGeodesicDistance(const Mesh& origTet, vtkPolyData *polyData, vtkIdType beginVertId, vtkIdType endVertId)
{
  std::vector<vtkIdType> ids;
  interpolator->InterpolateLine(reader->GetOutput(), beginVertId, endVertId,
          ids);

  double geodesicDistance = 0.0;
  int lastIndex = ids.at(0);
  for (int i = 1; i < ids.size(); i++)
  {
      int curIndex = ids.at(i);
      geodesicDistance += origTet.V.at(curIndex).Distance(origTet.V.at(lastIndex));
      lastIndex = curIndex;
  }
  return geodesicDistance;
}

double GetGeoDis(const Mesh& origTet, unsigned long vId1, unsigned long vId2)
{
    return ComputeGeodesicDistance(origTet, pd, vId1, vId2);
}

void GetNeiboringHexV(const Mesh& hexMesh, const std::vector<HexV_TetC_Paras>& hexV_tetC_paras)
{
    for (size_t i = 0; i < hexV_tetC_paras.size(); i++)
    {
        unsigned long hexVId = hexV_tetC_paras.at(i).hexVId;
        std::pair<std::multimap<unsigned long, unsigned long>::const_iterator,
                  std::multimap<unsigned long, unsigned long>::const_iterator> ret = hexMesh.VI_VI.equal_range(hexVId);
        for (std::multimap<unsigned long, unsigned long>::const_iterator iter = ret.first; iter != ret.second; ++iter)
        {
            neigboringHexV.push_back(iter->second);
        }
    }

    std::sort(neigboringHexV.begin(), neigboringHexV.end());
    std::vector<unsigned long>::iterator iter = std::unique(neigboringHexV.begin(), neigboringHexV.end());
    neigboringHexV.resize(std::distance(neigboringHexV.begin(), iter));


    for (std::vector<unsigned long>::iterator it = neigboringHexV.begin(); it != neigboringHexV.end();)
    {
        const unsigned long id = *it;
        if (std::find(hexvlist.begin(), hexvlist.end(), id) != hexvlist.end())
            it = neigboringHexV.erase(it);
        else
            ++it;
    }

    ///////////////////
    std::cout << "--------------- neigboringHexV --------------\n";
    for (size_t i = 0; i < neigboringHexV.size(); i++)
        std::cout << neigboringHexV[i] << " ";
    std::cout << endl;
}

unsigned long GetTargetSurfaceV(const Mesh& origTet, const TetHexParas& p)
{
    const unsigned long hexVId = neigboringHexV.at(0);
    glm::vec3 P;
    GetOrigSpaceLocation(origTet, p.tetIndex.at(hexVId), p.paras.at(hexVId), P);
    return GetClosestSurfacePoint(origTet, P);
}

void GetParas(const Mesh& origTet, const Mesh& polyTet, const Mesh& hexMesh, const TetHexParas& p, std::vector<HexV_TetC_Paras>& hexV_tetC_paras)
{
    unsigned long targetSurfaceVId = GetTargetSurfaceV(origTet, p);
    std::cout << "#### targetSurfaceVId = " << targetSurfaceVId << "\n";
    std::cout << "----------- tet index ----------- " << std::endl;
    double max_dis = -10000000000000;
    for (size_t i = 0; i < hexV_tetC_paras.size(); i++)
    {
        HexV_TetC_Paras& htp = hexV_tetC_paras.at(i);
        htp.locs.resize(htp.paras.size());
        htp.closestSufaceV.resize(htp.paras.size());
        htp.geoDis.resize(htp.paras.size());
        std::vector<unsigned long> closestSurfaceV;
        double min_distance = 10000000000000;
        double max_distance = -10000000000000;
        for (int j = 0; j < htp.tetCIds.size(); j++)
        {
            GetOrigSpaceLocation(origTet, htp.tetCIds.at(j), htp.paras.at(j), htp.locs.at(j));
            htp.closestSufaceV.at(j) = GetClosestSurfacePoint(origTet, htp.locs.at(j));
            htp.geoDis.at(j) = GetGeoDis(origTet, htp.closestSufaceV.at(j), targetSurfaceVId);
            if (htp.geoDis.at(j) > max_distance)
            {
                max_distance = htp.geoDis.at(j);
                htp.max_distance_tetCId = htp.tetCIds.at(j);
                htp.max_distance_paras_index = j;
                if (max_distance > max_dis)
                    max_dis = max_distance;
            }
            if (htp.geoDis.at(j) < min_distance)
            {
                min_distance = htp.geoDis.at(j);
                htp.min_distance_tetCId = htp.tetCIds.at(j);
                htp.min_distance_paras_index = j;
            }
        }
        std::cout << htp.min_distance_tetCId << " ";
    }
    std::cout << std::endl;
    std::cout << "----------- far tet index ----------- " << std::endl;
    for (size_t i = 0; i < hexV_tetC_paras.size(); i++)
    {
        HexV_TetC_Paras& htp = hexV_tetC_paras.at(i);
        std::cout << htp.max_distance_tetCId << " ";
    }
    std::cout << std::endl;
    std::cout << "----------- neigboringHexV ----------- " << std::endl;
    for (size_t i = 0; i < neigboringHexV.size(); i++)
    {
        unsigned long hexVid = neigboringHexV.at(i);
        std::cout << hexVid << " ";
        glm::vec3 loc;
        GetOrigSpaceLocation(origTet, p.tetIndex.at(hexVid), p.paras.at(hexVid), loc);
        unsigned long surfaceId = GetClosestSurfacePoint(origTet, loc);
        if (GetGeoDis(origTet, surfaceId, targetSurfaceVId) > max_dis/5)
            neigboring_FarNear.push_back(FAR);
        else
            neigboring_FarNear.push_back(NEAR);
    }
    std::cout << std::endl;

    std::cout << "----------- Far neigboringHexV ----------- " << std::endl;
    for (size_t i = 0; i < neigboring_FarNear.size(); i++)
    {
        if (neigboring_FarNear.at(i) == FAR)
            std::cout << neigboringHexV.at(i) << " ";
    }
    std::cout << std::endl;

    std::cout << "----------- Near neigboringHexV ----------- " << std::endl;
    for (size_t i = 0; i < neigboring_FarNear.size(); i++)
    {
        if (neigboring_FarNear.at(i) == NEAR)
            std::cout << neigboringHexV.at(i) << " ";
    }
    std::cout << std::endl;
}
struct Old_New
{
    unsigned long oldId;
    unsigned long newId;
};
std::vector<Old_New> old_news;
//////////////////////////////////////////////////////////
std::vector<Vertex> newHexV;
std::vector<Cell> newHexC;
void ModifyParas(Mesh& hexMesh, const std::vector<HexV_TetC_Paras>& hexV_tetC_paras, TetHexParas& p, const std::vector<unsigned long>& overlapC)
{
    for (size_t i = 0; i < hexV_tetC_paras.size(); i++)
    {
        const HexV_TetC_Paras& htp = hexV_tetC_paras.at(i);
        unsigned long id = htp.hexVId;
        p.paras.at(id) = htp.paras.at(htp.min_distance_paras_index);
        p.tetIndex.at(id) = htp.min_distance_tetCId;
    }

    for (size_t i = 0; i < hexV_tetC_paras.size(); i++)
    {
        const HexV_TetC_Paras& htp = hexV_tetC_paras.at(i);
        unsigned long id = p.paras.size();
        p.paras.push_back(htp.paras.at(htp.max_distance_paras_index));
        p.tetIndex.push_back(htp.max_distance_tetCId);
        Old_New on;
        on.oldId = htp.hexVId;
        on.newId = id;
        old_news.push_back(on);
    }

    ///////////////////
    std::cout << "--------------- old Id --------------\n";
    for (size_t i = 0; i < old_news.size(); i++)
        std::cout << old_news[i].oldId << " ";
    std::cout << endl;
    std::cout << "--------------- new Id --------------\n";
    for (size_t i = 0; i < old_news.size(); i++)
        std::cout << old_news[i].newId << " ";
    std::cout << endl;

    for (int i = 0; i < overlapC.size(); i++)
    {
        //get new cell id
        Cell cell(8);
        for (int j = 0; j < 8; j++)
        {
            unsigned long oldId = hexMesh.C.at(overlapC.at(i)).at(j);
            for (int k = 0; k < old_news.size(); k++)
                if (oldId == old_news.at(k).oldId){
                    cell.at(j) = old_news.at(k).newId;
                    break;
                }
        }
        newHexC.push_back(cell);
    }
    //////////////////////////////////////////////////////////
    //Find out those cell in which some vertices in overlap c and some vertices are in the far side respect to target Surface Point
    std::cout << "--------------- Modify Cell --------------\n";
    for (int i = 0; i < hexMesh.C.size(); i++)
    {
        Cell& cell = hexMesh.C.at(i);
        bool containNeignbor = false;
        bool containOverlap = false;
        for (int j = 0; j < cell.size(); j++)
        {
            unsigned long& id = cell.at(j);
            std::vector<unsigned long>::iterator iter = std::find(neigboringHexV.begin(), neigboringHexV.end(), id);
            int dis = std::distance(neigboringHexV.begin(), iter);
            if (iter != neigboringHexV.end() && neigboring_FarNear.at(dis) == FAR)
                containNeignbor = true;
            if (std::find(hexvlist.begin(), hexvlist.end(), id) != hexvlist.end())
                containOverlap = true;
        }

        if (containNeignbor && containOverlap)
        {
            std::cout << "--------------- Cell Index " <<  i << " --------------\n";
            for (size_t j = 0; j < cell.size(); j++)
                std::cout << cell[j] << " ";
            std::cout << endl;
            for (int j = 0; j < cell.size(); j++)
            {
                unsigned long& id = cell.at(j);
                if (std::find(hexvlist.begin(), hexvlist.end(), id) != hexvlist.end())
                    for (int k = 0; k < old_news.size(); k++)
                        if (id == old_news.at(k).oldId){
                            id = old_news.at(k).newId;
                            break;
                        }
            }
            std::cout << "--------------- After modification --------------\n";
            for (size_t j = 0; j < cell.size(); j++)
                std::cout << cell[j] << " ";
            std::cout << endl;
        }
    }
}

int main(int argc, char* argv[])
{
    if (argc != 4)
    {
        std::cout << "Usage: check_overlap <input_tet_mesh> <input_poly_mesh> <input_hex_mesh>" << std::endl;
        return -1;
    }
    MeshFileReader tetReader(argv[1]);
    MeshFileReader polyTetReader(argv[2]);
    MeshFileReader hexReader(argv[3]);
    Mesh tetMesh(tetReader.GetMesh());
    tetMesh.ExtractSurface();
    MeshFileWriter triMeshWriter(tetMesh.V, tetMesh.surface, "surface.tri.vtk", TRIANGLE);
    triMeshWriter.WriteVtkPolyDataFile();

    PolycubeMesh polycubeTetMesh(polyTetReader.GetMesh());
    Mesh hexMesh(hexReader.GetMesh());
    hexMesh.GetVI_VI();
    CheckDegeneracy(polycubeTetMesh);
    TetHexParas p;
    std::vector<unsigned long> overlapV;
    std::vector<HexV_TetC_Paras> hexV_tetC_paras;
    GetParameter(polycubeTetMesh, hexReader.GetMesh(), p, overlapV, hexV_tetC_paras);
    std::vector<unsigned long> overlapC;
    CheckOverlapCells(hexReader.GetMesh(), overlapV, overlapC);

    //////////////////////////////////////////////////////
    GetOverlapHtp(hexMesh, overlapC, hexV_tetC_paras);
    InitGeodesicDistance("surface.tri.vtk");
    GetNeiboringHexV(hexMesh, hexV_tetC_paras);
    GetParas(tetMesh, polycubeTetMesh, hexMesh, p, hexV_tetC_paras);
    ModifyParas(hexMesh, hexV_tetC_paras, p, overlapC);
    std::vector<Vertex> hexV;
    hexMesh.V.resize(hexvlist.size() + hexMesh.V.size());
    MapbackToOrigTet(tetMesh, hexMesh, p, hexV);
    std::vector<Cell> hexC = hexMesh.C;
    for (int i = 0; i < newHexC.size(); i++)
        hexC.push_back(newHexC.at(i));
    MeshFileWriter origHexMeshFileWriter(hexV, hexC, "clean.hex.vtk", HEXAHEDRA);
    origHexMeshFileWriter.WriteFile();
    //////////////////////////////////////////////////////

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

void CheckOverlapCells(const Mesh& hexMesh, const std::vector<unsigned long>& overlapV, std::vector<unsigned long>& overlapC)
{
    for (int i = 0; i < hexMesh.C.size(); i++)
    {
        const Cell& cell = hexMesh.C.at(i);
        bool bInside = true;
        for (int j = 0; j < cell.size(); j++)
        {
            if (std::find(overlapV.begin(), overlapV.end(), cell.at(j)) == overlapV.end()) // not found
            {
                bInside = false;
                break;
            }
        }
        if (bInside)
            overlapC.push_back(i);
    }
//    overlapC.erase(overlapC.begin());
//    overlapC.erase(overlapC.begin());
//    overlapC.erase(overlapC.begin());
    std::cout << "\n\n###Overlap cell ( " ;
    for (int k = 0; k < overlapC.size(); k++)
        std::cout << overlapC.at(k) << " " ;
    std::cout << ")" << std::endl;
}

void GetParameter(const Mesh& tetMesh, const Mesh& hexMesh, TetHexParas& p, std::vector<unsigned long>& overlapV, std::vector<HexV_TetC_Paras>& hexV_tetC_paras)
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
        std::vector<unsigned long> t;
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

            if (rambda.x >= 0 && rambda.y >= 0 && rambda.z>=0 && (rambda.x + rambda.y + rambda.z) <= 1.00000){
//            if (rambda.x > -1e-3 && rambda.y > -1e-3 && rambda.z>-1e-3 && (rambda.x + rambda.y + rambda.z) < 1.001){
                p.paras.at(i) = rambda;
                p.tetIndex.at(i) = j;
                p.flag.at(i) = true;
                flag = true;
                l.push_back(rambda);
                t.push_back(j);
                //break;
//                if (1 - rambda.x < 1e-8 || 1 - rambda.y < 1e-8 || 1 - rambda.z < 1e-8 || (rambda.x + rambda.y + rambda.z) < 1e-8);
//                else
                    if (std::find(degenerated.begin(), degenerated.end(), j) == degenerated.end())
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
            for (int k = 0; k < t.size(); k++)
                std::cout << t.at(k) << " " ;
            std::cout << ")" << std::endl;
            overlapV.push_back(i);

            //////////////////////////////
            HexV_TetC_Paras htp;
            htp.hexVId = i;
            htp.tetCIds = t;
            htp.paras = l;
            hexV_tetC_paras.push_back(htp);
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

