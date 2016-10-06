/*
 * polycube_hex_extract.cpp
 *
 *  Created on: Sep 9, 2016
 *      Author: cotrik
 */

#include "include.h"

// ************************** //
//       Data Structure       //
// ************************** //
struct VertexHex
{
    VertexHex()
            : v(Vertex()), tetId(0), vid(0), uvw(0.25, 0.25, 0.25)
    {

    }
    VertexHex(const Vertex& v)
            : v(v), tetId(0), vid(0), uvw(0.25, 0.25, 0.25)
    {

    }
    VertexHex(const Vertex& v, unsigned long tetId, const unsigned long vid, const glm::vec3& uvw)
            : v(v), tetId(tetId), vid(vid), uvw(uvw)
    {

    }
    VertexHex(const VertexHex& v)
            : v(v.v), tetId(v.tetId), vid(v.vid), uvw(v.uvw)
    {

    }
    VertexHex& operator =(const VertexHex& right)
    {
        if (this != &right)
        {
            v = right.v;
            tetId = right.tetId;
            vid = right.vid;
            uvw = right.uvw;
        }

        return *this;
    }
    const bool operator ==(const VertexHex& right) const
    {
        return v == right.v;
    }
    const bool operator <(const VertexHex& right) const
    {
        return v < right.v;
    }
    Vertex v;
    unsigned long tetId;
    unsigned long vid;
    glm::vec3 uvw;
};

typedef std::vector<VertexHex> Tet;
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
std::vector<int> cellField;
std::vector<int> pointField;

// ************************** //
//       Functions            //
// ************************** //
void InitGeodesicDistance(const char* polydataFilename);
double ComputeGeodesicDistance(const Mesh& M, vtkPolyData *polyData, vtkIdType beginVertId, vtkIdType endVertId);
double GetGeoDis(const Mesh& origTet, unsigned long vId1, unsigned long vId2);
void GenerateIntegerVertices(const PolycubeMesh& polycubeTetMesh, std::vector<VertexHex>& V);
void GetOverlapTetIds(const std::vector<VertexHex>& V, std::vector<unsigned long>& overlapTetIds);
void RemoveSharedVertex(std::vector<VertexHex>& V, std::vector<unsigned long>& overlapTetIds, const PolycubeMesh& M, std::vector<VertexHex>& sortedV);
void ExtactOverlapTetIds(const std::multimap<VertexHex, unsigned long>& V_T, std::vector<unsigned long>& overlapTetIds,
        const PolycubeMesh& M, std::vector<VertexHex>& needEraseEntries);
void EraseSharedEntries(std::multimap<VertexHex, unsigned long>& vertex_tetId_multimap, const std::vector<VertexHex>& needEraseEntries);
void GetOrigSpaceLocation(const Mesh& tetMesh, const unsigned long tetIndex, const glm::vec3& lambda, glm::vec3& new_r);
unsigned long GetClosestSurfacePoint(const Mesh& tetMesh, const glm::vec3& p);
unsigned long GetTargetSurfaceV(const Mesh& origTet, const VertexHex& p);

void Untangle(const Mesh& M,
        std::multimap<VertexHex, unsigned long>& vertex_tetId_multimap,
        const std::vector<std::vector<std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
        std::multimap<VertexHex, unsigned long>::const_iterator> > >& cubes,
        std::vector<Cell >& hexes);
void GetCubes(PolycubeMesh& polycubeTetMesh,
        const std::multimap<VertexHex, unsigned long>& vertex_tetId_multimap,
        std::vector<std::vector<std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
        std::multimap<VertexHex, unsigned long>::const_iterator> > >& cubes);
int main(int argc, char* argv[])
{
    if (argc != 3 && argc != 4)
    {
        std::cout << "Usage: polycube_hex_extract <polycube_tet_file> <polycube_hex_file> <orig_tet_file>" << std::endl;
        return -1;
    }
    MeshFileReader polycubeTetMeshFileReader(argv[1]);
    PolycubeMesh polycubeTetMesh(polycubeTetMeshFileReader.GetMesh());
    polycubeTetMesh.GetVI_VI();
    polycubeTetMesh.ExtractSurface();
    polycubeTetMesh.NormalizeCoordinateValue();
    polycubeTetMesh.GetMaxMinCoordinates();
    polycubeTetMesh.GetCorners();

    polycubeTetMesh.GetVertexInfo();
    polycubeTetMesh.GetFaceType();
    polycubeTetMesh.GetNormalOfSurfaceFaces();
    polycubeTetMesh.GetNormalOfSurfaceVertices();
    polycubeTetMesh.GetVertexInfo();

    polycubeTetMesh.GetFaceAndNeighborFaces();
    polycubeTetMesh.GetFacePatches_N();

    polycubeTetMesh.GetFacePatches();
    polycubeTetMesh.GetEdgePatches();
    polycubeTetMesh.GetVertexPatches();
    polycubeTetMesh.GetPatchesLabel();

    polycubeTetMesh.GetPatches();
    polycubeTetMesh.SortPatches();
    //polycubeTetMesh.ModifyPatchesPosition();

    // step 1. Round Surface of polycubeTetMesh
    polycubeTetMesh.Magnify(2.5e-3);
    //polycubeTetMesh.RoundSurface(2.5e-2);
    MeshFileWriter polycubeTetMeshFileWriter(polycubeTetMesh, "round.tet.vtk");
    polycubeTetMeshFileWriter.WriteFile();

    // step 2. Generate integer vertices in polycubeTetMesh's space
    std::vector<VertexHex> V;
    GenerateIntegerVertices(polycubeTetMesh, V);

    // step 3. Sort to find the overlap Vertices and Tets
    std::vector<unsigned long> overlapTetIds;
    //GetOverlapTetIds(V, overlapTetIds);  // this method is faster
//    std::vector<VertexHex> sortedV;
//    RemoveSharedVertex(V, overlapTetIds, polycubeTetMesh, sortedV);
//    V = sortedV;
//    overlapTetIds.clear();
    std::multimap<VertexHex, unsigned long> vertex_tetId_multimap;
    for (size_t i = 0; i < V.size(); i++)
        vertex_tetId_multimap.insert(std::pair<VertexHex, unsigned long>(V.at(i), V.at(i).tetId));
    std::vector<VertexHex> needEraseEntries;
    pointField.resize(V.size(), 0);
    ExtactOverlapTetIds(vertex_tetId_multimap, overlapTetIds, polycubeTetMesh, needEraseEntries);
    EraseSharedEntries(vertex_tetId_multimap, needEraseEntries);

    // step 4. extract polycubeHexMesh from polycubeTetMesh
    std::vector<std::vector<std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
    std::multimap<VertexHex, unsigned long>::const_iterator> > > cubes;
    std::cout << "GetCubes(polycubeTetMesh, vertex_tetId_multimap, cubes) " << std::endl;
    GetCubes(polycubeTetMesh, vertex_tetId_multimap, cubes);
    std::vector<Cell> hexes;
    std::vector<std::vector<VertexHex> > UVW;
    std::map<unsigned long, VertexHex> V_paras;
    std::vector<Vertex> VV2(V.size());
    if (argc == 4)
    {
        MeshFileReader origTetMeshFileReader(argv[3]);
        Mesh origTetMesh(origTetMeshFileReader.GetMesh());
        origTetMesh.ExtractSurface();
        MeshFileWriter origSurfaceMeshFileWriter(origTetMesh.V, origTetMesh.surface, "orig.surface.ref.vtk", TRIANGLE);
        origSurfaceMeshFileWriter.WriteVtkPolyDataFile();
        std::cout << "Untangling" << std::endl;
        Untangle(origTetMesh, vertex_tetId_multimap, cubes, hexes/*, UVW, V_paras*/);
        std::cout << "Mapping back " << std::endl;
        for (size_t i = 0; i < V.size(); i++)
        {
            glm::vec3 p;
            GetOrigSpaceLocation(origTetMesh, V[i].tetId, V[i].uvw, p);
            VV2[i] = p;
        }
    }
    std::cout << "GetHexElements(cubes, C) " << std::endl;
    Cell cell(8, 0);
    std::vector<Cell> C(hexes.size(), cell);
    // GetHexElements(cubes, C)
    {
//        for (size_t i = 0; i < cubes.size(); i++)
//        {
//            if (i % 10000 == 0)
//            std::cout << "Generate hex " << i << "-->" << cubes.size() << std::endl;
//            Cell hex(8, 0);
//            for (size_t j = 0; j < 8; j++)
//            {
////                const std::multimap<VertexHex, unsigned long>::const_iterator iterBegin = vertex_tetId_multimap.begin();
////                const std::multimap<VertexHex, unsigned long>::const_iterator iter = cubes[i][j].first;
////                C[i][j] = std::distance(iterBegin, iter);
//                C[i][j] = cubes[i][j].first->first.vid;
//            }
//        }
        for (size_t i = 0; i < hexes.size(); i++)
        {
            if (i % 100000 == 0)
            std::cout << "Generate hex " << i << "-->" << hexes.size() << std::endl;
            Cell hex(8, 0);
            for (size_t j = 0; j < 8; j++)
            {
                C[i][j] = hexes[i][j];
            }
        }
    }

    std::vector<Vertex> VV(V.size());
    for (size_t i = 0; i < VV.size(); i++)
    {
        VV[i] = V[i].v;
    }

    MeshFileWriter polycubeHexMeshFileWriter(VV, C, argv[2], HEXAHEDRA);
    polycubeHexMeshFileWriter.WriteFile();
    polycubeHexMeshFileWriter.WritePointData(pointField);
    polycubeHexMeshFileWriter.WriteCellData(cellField);

    MeshFileWriter polycubeHexMeshFileWriter2(VV2, C, "orig.hex.vtk", HEXAHEDRA);
    polycubeHexMeshFileWriter2.WriteFile();
    polycubeHexMeshFileWriter2.WriteCellData(cellField);
    return 0;
}

static void GetMaxMinVertex(const Mesh& M, const Cell& tet, Vertex& max_v, Vertex& min_v)
{
    min_v = M.V[tet[0]];
    max_v = M.V[tet[0]];

    GeoUtil::GetMaxMinCoordinateValueOfVertex(M.V[tet[1]], max_v, min_v);
    GeoUtil::GetMaxMinCoordinateValueOfVertex(M.V[tet[2]], max_v, min_v);
    GeoUtil::GetMaxMinCoordinateValueOfVertex(M.V[tet[3]], max_v, min_v);
}
void GenerateIntegerVertices(const PolycubeMesh& polycubeTetMesh, std::vector<VertexHex>& V)
{
    //std::ofstream ofs("overlap.txt");
    std::cout << "------------------------------------ " << std::endl;
    const PolycubeMesh& M = polycubeTetMesh;
    for (size_t i = 0; i < M.C.size(); i++)
    {
        if (i % 1000 == 0)
        std::cout << "GenerateIntegerVertices for tet " << i << std::endl;
        const Cell& tet = M.C.at(i);
        Vertex max_v = M.V[tet[0]];
        Vertex min_v = M.V[tet[0]];
        GetMaxMinVertex(M, tet, max_v, min_v);
//        std::cout << "minv " << min_v.x << " " << min_v.y << " " << min_v.z << std::endl;
//        std::cout << "maxv " << max_v.x << " " << max_v.y << " " << max_v.z << std::endl << std::endl;
        for (int x = min_v.x - 0.5; x < max_v.x + 0.5; x++) {
            for (int y = min_v.y - 0.5; y < max_v.y + 0.5; y++) {
                for (int z = min_v.z - 0.5; z < max_v.z + 0.5; z++)
                {
                    const Vertex v(x, y, z);

                    const Vertex& v0 = M.V[tet[0]];
                    const Vertex& v1 = M.V[tet[1]];
                    const Vertex& v2 = M.V[tet[2]];
                    const Vertex& v3 = M.V[tet[3]];
                    glm::vec3 uvw;
                   // if (GeoUtil::IsVertexInsideTetrahedron_Robust(v, v0, v1, v2, v3, uvw))
                    if (GeoUtil::IsPointInTetrahedron(v, v0, v1, v2, v3, uvw))
                    {
                        //ofs << x << " " << y << " " << z << "\t" << i << std::endl;
                        VertexHex hexV(v, i, V.size(), uvw);
                        V.push_back(hexV);
                    }
                }
            }
        }
    }
    std::cout << "------------------------------------ " << std::endl;
}

void GetOverlapTetIds(const std::vector<VertexHex>& V, std::vector<unsigned long>& overlapTetIds)
{
    std::vector<VertexHex> sortedV(V);
    std::sort(sortedV.begin(), sortedV.end());
    std::vector<VertexHex> overlapV;
    //std::vector<unsigned long> overlapTetIds;
    for (size_t i = 0; i < sortedV.size() - 1; i++)
    {
        if (sortedV[i] == sortedV[i+1] && sortedV[i].tetId != sortedV[i+1].tetId)
        {
            overlapV.push_back(sortedV[i]);
            overlapTetIds.push_back(sortedV[i].tetId);
            overlapTetIds.push_back(sortedV[i + 1].tetId);
        }
    }
    std::sort(overlapTetIds.begin(), overlapTetIds.end());
    std::vector<unsigned long>::iterator iterTetIds = std::unique(overlapTetIds.begin(), overlapTetIds.end());
    overlapTetIds.resize(std::distance(overlapTetIds.begin(), iterTetIds));

    std::cout << "**** Overlap tetIds ****\n";
    for (size_t i = 0; i < overlapTetIds.size(); i++)
    {
        std::cout << overlapTetIds.at(i) << " ";
    }
    std::cout << "**** ******** ****\n";
}

void RemoveSharedVertex(std::vector<VertexHex>& V, std::vector<unsigned long>& overlapTetIds, const PolycubeMesh& M, std::vector<VertexHex>& sortedV)
{
    //std::vector<VertexHex> sortedV(V);
    sortedV = V;
    std::sort(sortedV.begin(), sortedV.end());
    std::vector<VertexHex> overlapV;
    //std::vector<unsigned long> overlapTetIds;
    for (std::vector<VertexHex>::iterator iter = sortedV.begin(); iter != sortedV.end() - 1 && iter != sortedV.end();)
    {
        static int i = 0;
        if (i++ % 1000 == 0)
            std::cout << "processing V " << i << std::endl;
        if ((*iter) == *(iter+1) && (*iter).tetId != (*(iter+1)).tetId)
        {
            bool share = false;
            std::vector<unsigned long> vids;
                const unsigned long tetid1 = (*iter).tetId;
                const Cell& tet1 = M.C.at(tetid1);
                for(int j = 0; j < 4; j++)
                    vids.push_back(tet1.at(j));

                const unsigned long tetid2 = (*(iter+1)).tetId;
                const Cell& tet2 = M.C.at(tetid2);
                for(int j = 0; j < 4; j++)
                    vids.push_back(tet2.at(j));

                std::sort(vids.begin(), vids.end());
                std::vector<unsigned long>::const_iterator iterV = std::unique(vids.begin(), vids.end());
                if (iterV != vids.end())
                    share = true;

            if (share)
            {
                iter = sortedV.erase(iter);
                continue;
            }
            else
            {
                overlapV.push_back((*iter));
                overlapTetIds.push_back((*iter).tetId);
                overlapTetIds.push_back((*(iter+1)).tetId);
            }
        }
        //else
        {
            ++iter;
        }
    }
    std::sort(overlapTetIds.begin(), overlapTetIds.end());
    std::vector<unsigned long>::iterator iterTetIds = std::unique(overlapTetIds.begin(), overlapTetIds.end());
    overlapTetIds.resize(std::distance(overlapTetIds.begin(), iterTetIds));

    std::cout << "**** Overlap tetIds ****\n";
    for (size_t i = 0; i < overlapTetIds.size(); i++)
    {
        std::cout << overlapTetIds.at(i) << std::endl;
    }
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

double ComputeGeodesicDistance(const Mesh& M, vtkPolyData *polyData, vtkIdType beginVertId, vtkIdType endVertId)
{
    std::vector<vtkIdType> ids;
    interpolator->InterpolateLine(reader->GetOutput(), beginVertId, endVertId, ids);

    double geodesicDistance = 0.0;
    int lastIndex = ids.at(0);
    for (int i = 1; i < ids.size(); i++)
    {
        int curIndex = ids.at(i);
        geodesicDistance += M.V.at(curIndex).Distance(M.V.at(lastIndex));
        lastIndex = curIndex;
    }
    return geodesicDistance;
}

double GetGeoDis(const Mesh& origTet, unsigned long vId1, unsigned long vId2)
{
    return ComputeGeodesicDistance(origTet, pd, vId1, vId2);
}

static void GetDuplicated(const std::vector<unsigned long>& vids, std::vector<unsigned long>& duplicated_vids)
{
    for (size_t i = 0; i < vids.size() - 1; i++)
        if (vids[i] == vids[i+1] && std::find(duplicated_vids.begin(), duplicated_vids.end(), vids[i]) == duplicated_vids.end())
            duplicated_vids.push_back(vids[i]);
}

bool FindConnectivtyInGroups(const PolycubeMesh& M, const unsigned long vid, std::vector<std::vector<unsigned long> >& gs)
{
    bool bRet = false;
    for (size_t i = 0; i < gs.size(); i++)
    {
        for (size_t j = 0; j < gs[i].size(); j++)
        {
            const std::pair<unsigned long, unsigned long> p(vid, gs[i][j]);
            std::pair<std::multimap<unsigned long, unsigned long>::const_iterator,
            std::multimap<unsigned long, unsigned long>::const_iterator> ret = M.VI_VI.equal_range(gs[i][j]);
            std::multimap<unsigned long, unsigned long>::const_iterator iter = ret.first;
            for (; iter != ret.second; ++iter)
            {
                if (iter->second == vid)
                {
                    gs[i].push_back(vid);
                    bRet = true;
                    break;
                }
            }
            if (bRet)
                break;
        }
        if (bRet)
            break;
    }

    return bRet;
}

static void Group(const PolycubeMesh& M, const std::vector<unsigned long>& vs, std::vector<std::vector<unsigned long> >& gs)
{
    std::vector<unsigned long> V = vs;
    std::vector<unsigned long> g;
//    g.push_back(*V.begin());
//    gs.push_back(g);
//    V.erase(V.begin());
    while (!V.empty())
    {
        bool hasConnection = false;
        std::vector<unsigned long>::iterator iter = V.begin();
        for (; iter != V.end();)
        {
            if (FindConnectivtyInGroups(M, *iter, gs))
            {
                hasConnection = true;
                //g.push_back(*iter);
                iter = V.erase(iter);
            }
            else
                ++iter;
        }
        if (!hasConnection){
//            gs.push_back(g);
//            g.clear();
            if (!V.empty())
            {
                g.clear();
                g.push_back(*V.begin());
                gs.push_back(g);
                V.erase(V.begin());
            }
        }
    }

//    for (size_t i = 0; i < vs.size(); i++)
//    {
//        const unsigned long vid = vs.at(i);
//        if (!FindConnectivtyInGroups(M, vid, gs))
//        {
//            std::vector<unsigned long> g;
//            g.push_back(vid);
//            gs.push_back(g);
//        }
//    }
}

std::ofstream ofs2("log2.txt");

void ProcessGroupVertexHex(const std::vector<std::vector<VertexHex> >& groupVertexHex,
        const std::multimap<VertexHex, unsigned long>& V_T, std::vector<unsigned long>& overlapTetIds,
        const PolycubeMesh& M, std::vector<VertexHex>& needEraseEntries)
{

    for (int i = 0; i < groupVertexHex.size(); ++i){
        bool share = false;
        std::vector<unsigned long> vids;
        for (int j = 0; j < groupVertexHex[i].size(); ++j){
            const unsigned long tetid = groupVertexHex[i][j].tetId;
            const Cell& tet = M.C.at(tetid);
            for(int i = 0; i < 4; i++)
                vids.push_back(tet.at(i));
        }
        std::sort(vids.begin(), vids.end());
        std::vector<unsigned long> duplicated_vids;
        GetDuplicated(vids, duplicated_vids);
        std::vector<unsigned long>::iterator iterV = std::unique(vids.begin(), vids.end());
        if (iterV != vids.end())
        {
            share = true;
        }
        ofs2 << "vids = ";
        for (int i = 0; i < vids.size(); ++i){
            ofs2 << vids[i] << " ";
        }
        ofs2 << endl;
        ofs2 << "duplicated_vids = ";
        for (int i = 0; i < duplicated_vids.size(); ++i){
            ofs2 << duplicated_vids[i] << " ";
        }
        ofs2 << endl << endl;

        if (!share)
        {
            //if (vids.size() > 4)
            for (int j = 0; j < groupVertexHex[i].size(); ++j){
                overlapTetIds.push_back(groupVertexHex[i][j].tetId);
                pointField[groupVertexHex[i][j].vid] = 2;
            }
        }
        else
        {
            pointField[groupVertexHex[i][0].vid] = 2;
            for (int j = 1; j < groupVertexHex[i].size(); ++j){
                needEraseEntries.push_back(groupVertexHex[i][j]);
                pointField[groupVertexHex[i][j].vid] = 1;
            }
        }
    }

//    for (int i = 0; i < groupVertexHex.size(); ++i)
//    {
//        pointField[groupVertexHex[i][0].vid] = 2;
//        for (int j = 1; j < groupVertexHex[i].size(); ++j)
//        {
//            needEraseEntries.push_back(groupVertexHex[i][j]);
//            pointField[groupVertexHex[i][j].vid] = 2;
//        }
//    }
}

static int GetGroupIndex(const PolycubeMesh& M, const std::vector<std::vector<unsigned long> >& groups,
        std::multimap<VertexHex, unsigned long>::const_iterator iter)
{
    int groupId = 0;
    for (int i = 0; i < groups.size(); ++i)
    {
        for (int j = 0; j < groups[i].size(); ++j)
        {
            const unsigned long tetid = iter->second;
            const Cell& tet = M.C.at(tetid);
            for (int k = 0; k < 4; k++)
            {
                if (groups[i][j] == tet[k])
                {
                    groupId = i;
                    return i;
                }
            }
        }
    }

    return groupId;
}

static void GroupTet(const PolycubeMesh& M, const std::vector<std::vector<unsigned long> >& groups,
        const std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
                        std::multimap<VertexHex, unsigned long>::const_iterator>& ret,
        std::vector<std::vector<VertexHex> >& groupVertexHex)
{
    groupVertexHex.resize(groups.size());
    for (std::multimap<VertexHex, unsigned long>::const_iterator iter = ret.first; iter != ret.second; ++iter)
    {
        int gId = GetGroupIndex(M, groups, iter);
        groupVertexHex[gId].push_back(iter->first);
    }
}

static void LogPrint(std::ofstream& ofs, const std::vector<std::vector<unsigned long> >& groups,
        const std::vector<unsigned long>& vids,
        const std::vector<unsigned long>& unique_vids,
        const std::vector<unsigned long>& duplicated_vids,
        const std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
                        std::multimap<VertexHex, unsigned long>::const_iterator>& ret)
{
    ofs << "groups = " << groups.size()<< endl;
    for (int i = 0; i < groups.size(); ++i){
        ofs << "group = " << i << " --> ";
        for (int j = 0; j < groups[i].size(); ++j)
            ofs << groups[i][j] << " ";
        ofs << endl;
    }

    ofs << "vids = ";
    for (int i = 0; i < vids.size(); ++i){
        ofs << vids[i] << " ";
    }
    ofs << endl;

    ofs << "unique_vids = ";
    for (int i = 0; i < unique_vids.size(); ++i){
        ofs << unique_vids[i] << " ";
    }
    ofs << endl;

    ofs << "duplicated_vids = ";
    for (int i = 0; i < duplicated_vids.size(); ++i){
        ofs << duplicated_vids[i] << " ";
    }
    ofs << endl;

    ofs << "tetIds = ";
    for (std::multimap<VertexHex, unsigned long>::const_iterator iter = ret.first; iter != ret.second; ++iter)
    {
        ofs << iter->second << " ";
    }
    ofs << endl << endl;
}
void ExtactOverlapTetIds(const std::multimap<VertexHex, unsigned long>& V_T, std::vector<unsigned long>& overlapTetIds,
        const PolycubeMesh& M, std::vector<VertexHex>& needEraseEntries)
{
    std::ofstream ofs("log.txt");
    const std::multimap<VertexHex, unsigned long>::const_iterator iterBegin = V_T.begin();
    const std::multimap<VertexHex, unsigned long>::const_iterator iterEnd = V_T.end();
    std::multimap<VertexHex, unsigned long>::const_iterator iter = iterBegin;
    for (; iter != iterEnd;)
    {
        const std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
                std::multimap<VertexHex, unsigned long>::const_iterator> ret =
                V_T.equal_range(iter->first);

        const int count = std::distance(ret.first, ret.second);
        if (count > 1){
//            bool share = false;
            std::vector<unsigned long> vids;
            for (iter = ret.first; iter != ret.second; ++iter)
            {
                const unsigned long tetid = iter->second;
                const Cell& tet = M.C.at(tetid);
                for(int i = 0; i < 4; i++)
                    vids.push_back(tet.at(i));
            }
            std::sort(vids.begin(), vids.end());
            std::vector<unsigned long> duplicated_vids;
            GetDuplicated(vids, duplicated_vids);

            std::vector<unsigned long>::iterator iterV = std::unique(vids.begin(), vids.end());
            if (iterV != vids.end())
            {
                std::vector<unsigned long> unique_vids;
                std::copy(vids.begin(), iterV, back_inserter(unique_vids));
                std::vector<std::vector<unsigned long> > groups;
                Group(M, unique_vids, groups);
                if (groups.size() > 1)
                {
                    LogPrint(ofs, groups, vids, unique_vids, duplicated_vids, ret);
                    std::vector<std::vector<VertexHex> > groupVertexHex(groups.size());
                    GroupTet(M, groups, ret, groupVertexHex);
                    ProcessGroupVertexHex(groupVertexHex, V_T, overlapTetIds, M, needEraseEntries);
                }
                else //if (groups.size() == 1)
                {
                    iter = ret.first;
                    pointField[iter->first.vid] = 1;
                    ++iter;
                    for (; iter != ret.second; ++iter){
                        needEraseEntries.push_back(iter->first);
                        pointField[iter->first.vid] = 1;
                    }
                }
            }
            else
            {
                for (iter = ret.first; iter != ret.second; ++iter){
                    overlapTetIds.push_back(iter->second);
                    pointField[iter->first.vid] = 2;
                }
            }
        }
        iter = ret.second;
    }

    std::sort(overlapTetIds.begin(), overlapTetIds.end());
    std::vector<unsigned long>::iterator iterTetIds = std::unique(overlapTetIds.begin(), overlapTetIds.end());
    overlapTetIds.resize(std::distance(overlapTetIds.begin(), iterTetIds));

    std::cout << "**** Overlap tetIds ****\n";
    for (size_t i = 0; i < overlapTetIds.size(); i++)
    {
        std::cout << overlapTetIds.at(i) << " ";
    }
    std::cout << "\n**** ******** ****\n";
}

void EraseSharedEntries(std::multimap<VertexHex, unsigned long>& vertex_tetId_multimap, const std::vector<VertexHex>& needEraseEntries)
{
    std::cout << "---- Erase before " << vertex_tetId_multimap.size() << std::endl;
    for (int i = 0; i< needEraseEntries.size(); i++)
    {
        const VertexHex& vh = needEraseEntries.at(i);
        std::pair<std::multimap<VertexHex, unsigned long>::iterator,
            std::multimap<VertexHex, unsigned long>::iterator> ret = vertex_tetId_multimap.equal_range(vh);
        std::multimap<VertexHex, unsigned long>::iterator iterBegin = ret.first;
        std::multimap<VertexHex, unsigned long>::iterator iterEnd = ret.second;
        std::multimap<VertexHex, unsigned long>::iterator iter = iterBegin;
        for (; iter != iterEnd;)
        {
            if (iter->first == vh && iter->first.vid == vh.vid && iter->first.tetId == vh.tetId/* && iter->second == vh.tetId*/)
                iter = vertex_tetId_multimap.erase(iter);
            else
                ++iter;
        }
    }
    std::cout << "---- Erase After " << vertex_tetId_multimap.size() << std::endl;
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

unsigned long GetTargetSurfaceV(const Mesh& origTet, const VertexHex& p)
{
    glm::vec3 P;
    GetOrigSpaceLocation(origTet, p.tetId, p.uvw, P);
    return GetClosestSurfacePoint(origTet, P);
}

//void Untangle(const Mesh& M,
//        std::multimap<VertexHex, unsigned long>& vertex_tetId_multimap,
//        const std::vector<std::vector<std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
//        std::multimap<VertexHex, unsigned long>::const_iterator> > >& cubes,
//        std::vector<Cell >& hexes, std::vector<std::vector<VertexHex> >& UVW,
//        std::map<unsigned long, VertexHex>& V_paras)
//{
//    InitGeodesicDistance("orig.surface.ref.vtk");
//
//    for (int i = 0; i < cubes.size(); i++)
//    {
//        int overlapCount = 0;
//        bool bOverlap[8] = {false};
//        for (int j = 0; j < 8; j++)
//            if (vertex_tetId_multimap.count(cubes[i][j].first->first) > 1)
//            {
//                bOverlap[j] = true;
//                overlapCount++;
//            }
////        if (overlapCount != 0)
////            std::cout << "Hex element " << i << " overlapCount = " << overlapCount << std::endl;
//        if (overlapCount == 0) // has no Overlap
//        {
//            Cell hex(8);
//            std::vector<VertexHex> uvw(8);
//            for (int j = 0; j < 8; j++){
//                hex[j] = cubes[i][j].first->first.vid;
//                uvw[j] = cubes[i][j].first->first;
//                V_paras[hex[j]] = uvw[j];
//            }
//            hexes.push_back(hex);
//            UVW.push_back(uvw);
//        }
//
//        if (overlapCount > 0 && overlapCount < 8) // Partially Overlap
//        {
//            unsigned long RefTetId = cubes[i][0].first->second;
//            int j = 0;
//            for (; j < 8; j++)
//                if (!bOverlap[j]){
//                    RefTetId = cubes[i][j].first->second;
//                    break;
//                }
//
//            unsigned long vid = GetTargetSurfaceV(M, cubes[i][j].first->first);
//            Cell hex(8);
//            std::vector<VertexHex> uvw(8);
//            for (j = 0; j < 8; j++)
//            {
//                if (!bOverlap[j])
//                {
//                    hex[j] = cubes[i][j].first->first.vid;
//                }
//                else
//                {
//                    double min_distance = 10000000000000;
//                    double max_distance = -10000000000000;
//                    std::multimap<VertexHex, unsigned long>::const_iterator iterPrefer = cubes[i][j].first;
//                    std::multimap<VertexHex, unsigned long>::const_iterator iter = cubes[i][j].first;
//                    for (; iter != cubes[i][j].second; ++iter)
//                    {
//                        unsigned long tetId = cubes[i][j].first->second;
//                        unsigned long v = GetTargetSurfaceV(M, cubes[i][j].first->first);
//                        double dis = GetGeoDis(M, v, vid);
//                        if (dis < min_distance)
//                        {
//                            min_distance = dis;
//                            iterPrefer = iter;
//                        }
//                    }
//                    hex[j] = iterPrefer->first.vid;
//                    uvw[j] = iterPrefer->first;
//                    V_paras[hex[j]] = uvw[j];
//                }
//            }
//            hexes.push_back(hex);
//            UVW.push_back(uvw);
//        }
//
//        if (overlapCount == 8) // Partially Overlap
//        {
//            unsigned long vid = GetTargetSurfaceV(M, cubes[i][0].first->first);
//            Cell hexMin(8);
//            Cell hexMax(8);
//            std::vector<VertexHex> uvwMin(8);
//            std::vector<VertexHex> uvwMax(8);
//            for (int j = 0; j < 8; j++)
//            {
//                double min_distance = 10000000000000;
//                double max_distance = -10000000000000;
//                std::multimap<VertexHex, unsigned long>::const_iterator iterPreferMin = cubes[i][j].first;
//                std::multimap<VertexHex, unsigned long>::const_iterator iterPreferMax = cubes[i][j].first;
//                std::multimap<VertexHex, unsigned long>::const_iterator iter = cubes[i][j].first;
//                for (; iter != cubes[i][j].second; ++iter)
//                {
//                    unsigned long tetId = cubes[i][j].first->second;
//                    unsigned long v = GetTargetSurfaceV(M, cubes[i][j].first->first);
//                    double dis = GetGeoDis(M, v, vid);
//                    if (dis < min_distance)
//                    {
//                        min_distance = dis;
//                        iterPreferMin = iter;
//                    }
//                    if (dis > max_distance)
//                    {
//                        max_distance = dis;
//                        iterPreferMax = iter;
//                    }
//                }
//                hexMin[j] = iterPreferMin->first.vid;
//                hexMax[j] = iterPreferMax->first.vid;
//                uvwMin[j] = iterPreferMin->first;
//                uvwMax[j] = iterPreferMax->first;
//                V_paras[hexMin[j]] = uvwMin[j];
//                V_paras[hexMax[j]] = uvwMax[j];
//            }
//            hexes.push_back(hexMin);
//            hexes.push_back(hexMax);
//            UVW.push_back(uvwMin);
//            UVW.push_back(uvwMax);
//        }
//    }
//}

void GetCubes(PolycubeMesh& polycubeTetMesh,
        const std::multimap<VertexHex, unsigned long>& vertex_tetId_multimap,
        std::vector<std::vector<std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
        std::multimap<VertexHex, unsigned long>::const_iterator> > >& cubes)
{
    const PolycubeMesh& M = polycubeTetMesh;
    polycubeTetMesh.GetMaxMinCoordinates();
    const Vertex& max_v = M.m_maxVertex;
    const Vertex& min_v = M.m_minVertex;
    for (double x = min_v.x + 0.5; x < max_v.x; x++) {
        for (double y = min_v.y + 0.5; y < max_v.y; y++) {
            for (double z = min_v.z + 0.5; z < max_v.z; z++)
            {
                const Vertex v0(x, y, z);
                const VertexHex vc(v0); // center of cube
                const VertexHex v[8] =
                {
                    Vertex(round(x - 0.5), round(y - 0.5), round(z + 0.5)),
                    Vertex(round(x + 0.5), round(y - 0.5), round(z + 0.5)),
                    Vertex(round(x + 0.5), round(y - 0.5), round(z - 0.5)),
                    Vertex(round(x - 0.5), round(y - 0.5), round(z - 0.5)),
                    Vertex(round(x - 0.5), round(y + 0.5), round(z + 0.5)),
                    Vertex(round(x + 0.5), round(y + 0.5), round(z + 0.5)),
                    Vertex(round(x + 0.5), round(y + 0.5), round(z - 0.5)),
                    Vertex(round(x - 0.5), round(y + 0.5), round(z - 0.5))
                };
                bool flag = true;
                std::vector<std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
                                std::multimap<VertexHex, unsigned long>::const_iterator> >vids(8);
                for (int i = 0; i < 8; i++)
                {
                    const std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
                          std::multimap<VertexHex, unsigned long>::const_iterator> ret = vertex_tetId_multimap.equal_range(v[i]);
                    if (ret.first == ret.second)
                    {
                        flag = false;
                        break;
                    }
                    vids[i] = ret;
                }

                if (flag)
                {
                    cubes.push_back(vids);
                }
            }
        }
    }
}


void Untangle(const Mesh& M,
        std::multimap<VertexHex, unsigned long>& vertex_tetId_multimap,
        const std::vector<std::vector<std::pair<std::multimap<VertexHex, unsigned long>::const_iterator,
        std::multimap<VertexHex, unsigned long>::const_iterator> > >& cubes,
        std::vector<Cell >& hexes)
{
    InitGeodesicDistance("orig.surface.ref.vtk");

    for (int i = 0; i < cubes.size(); i++)
    {
        int overlapCount = 0;
        bool bOverlap[8] = {false};
        for (int j = 0; j < 8; j++)
            if (vertex_tetId_multimap.count(cubes[i][j].first->first) > 1
                    && cubes[i][j].first->second != cubes[i][j].second->second)
            {
                bOverlap[j] = true;
                overlapCount++;
            }
//        if (overlapCount < 8 && overlapCount > 0)
//            std::cout << "Hex element " << i << " overlapCount = " << overlapCount << std::endl;
        if (overlapCount == 0) // has no Overlap
        {
            Cell hex(8);
            for (int j = 0; j < 8; j++){
                hex[j] = cubes[i][j].first->first.vid;
            }
            hexes.push_back(hex);
            cellField.push_back(0);
        }

        if (overlapCount > 0 && overlapCount < 8) // Partially Overlap
        {
            unsigned long RefTetId = cubes[i][0].first->second;
            int j = 0;
            for (; j < 8; j++)
                if (!bOverlap[j]){
                    RefTetId = cubes[i][j].first->second;
                    break;
                }

            unsigned long vid = GetTargetSurfaceV(M, cubes[i][j].first->first);
            Cell hex(8);
            //std::cout << "  origTet target vid = " << vid << std::endl;
            for (j = 0; j < 8; j++)
            {
                if (!bOverlap[j])
                {
                    hex[j] = cubes[i][j].first->first.vid;
                }
                else
                {
                    double min_distance = 10000000000000;
                    double max_distance = -10000000000000;
                    std::multimap<VertexHex, unsigned long>::const_iterator iterPrefer = cubes[i][j].first;
                    std::multimap<VertexHex, unsigned long>::const_iterator iter = cubes[i][j].first;
                    for (; iter != cubes[i][j].second; ++iter)
                    {
                        unsigned long tetId = iter->second;
                        unsigned long v = GetTargetSurfaceV(M, iter->first);
                        double dis = GetGeoDis(M, v, vid);
                        if (dis < min_distance)
                        {
                            min_distance = dis;
                            iterPrefer = iter;
                        }
                        //std::cout << "        tetId = " << tetId << " v = " << v << " dis = " << dis << " hex-vid = " << iter->first.vid << std::endl;
                    }
                    hex[j] = iterPrefer->first.vid;
                }
            }
            hexes.push_back(hex);
            cellField.push_back(1);
        }

        if (overlapCount == 8) // Fully Overlap
        {
            unsigned long vid = GetTargetSurfaceV(M, cubes[i][0].first->first);
            Cell hexMin(8);
            Cell hexMax(8);
            //std::cout << "  origTet target vid = " << vid << std::endl;
            for (int j = 0; j < 8; j++)
            {
                double min_distance = 10000000000000;
                double max_distance = -10000000000000;
                std::multimap<VertexHex, unsigned long>::const_iterator iterPreferMin = cubes[i][j].first;
                std::multimap<VertexHex, unsigned long>::const_iterator iterPreferMax = cubes[i][j].first;
                std::multimap<VertexHex, unsigned long>::const_iterator iter = cubes[i][j].first;
                //std::cout << "    j = " << j << std::endl;
                for (; iter != cubes[i][j].second; ++iter)
                {
                    unsigned long tetId = iter->second;
                    unsigned long v = GetTargetSurfaceV(M, iter->first);
                    double dis = GetGeoDis(M, v, vid);
                    if (dis < min_distance)
                    {
                        min_distance = dis;
                        iterPreferMin = iter;
                    }
                    if (dis > max_distance)
                    {
                        max_distance = dis;
                        iterPreferMax = iter;
                    }
//                    std::cout << "        tetId = " << tetId << " v = " << v << " dis = " << dis << " hex-vid = " << iter->first.vid << std::endl;
                }
                hexMin[j] = iterPreferMin->first.vid;
                hexMax[j] = iterPreferMax->first.vid;
            }

//            std::cout << "hexMin = ";
//            for (int j = 0; j < 8; j++)
//                std::cout << hexMin[j] << " ";
//            std::cout << "\n";
//
//            std::cout << "hexMax = ";
//            for (int j = 0; j < 8; j++)
//                std::cout << hexMax[j] << " ";
//            std::cout << "\n";

            hexes.push_back(hexMin);
            hexes.push_back(hexMax);
            cellField.push_back(2);
            cellField.push_back(2);
        }
    }
}
