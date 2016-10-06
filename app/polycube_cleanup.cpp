#include "include.h"

//int main(int argc, char* argv[])
//{
//	return PolycubeCleanup(argc, argv);
//}

int main(int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "Usage: polycube_cleanup <input.off> <output.off>" << std::endl;
		return -1;
	}

	MeshFileReader polycubeTriMeshFileReader(argv[1], MESH_TYPE_TRIANGLE_OFF);
	PolycubeMesh polycubeTriMesh(polycubeTriMeshFileReader.GetMesh());
	////////
	polycubeTriMesh.GetMaxMinCoordinates();
	polycubeTriMesh.ExtractSurface();
	polycubeTriMesh.GetFaceType();
	polycubeTriMesh.GetFaceAndNeighborFaces();
	polycubeTriMesh.GetFacePatches_N();
	polycubeTriMesh.GetNormalOfSurfaceFaces();
	polycubeTriMesh.GetNormalOfSurfaceVertices();
	polycubeTriMesh.GetFacePatches();
	polycubeTriMesh.GetEdgePatches();
	polycubeTriMesh.GetVertexPatches();
	//////////////////////////////
	polycubeTriMesh.GetPatches();
	polycubeTriMesh.SortPatches();
	polycubeTriMesh.ModifyPatchesPosition();
	polycubeTriMesh.SmoothPatches();
	MeshFileWriter polycubeTriMeshFileWriter(polycubeTriMesh, argv[2], MESH_TYPE_TRIANGLE_OFF);
	polycubeTriMeshFileWriter.WriteMeshFile(argv[2], MESH_TYPE_TRIANGLE_OFF);

	return 0;
}
