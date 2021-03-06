#ifndef __Mesh_File_Reader_H__
#define __Mesh_File_Reader_H__

#include <vector>
#include <string>
#include <fstream>
#include "Mesh.h"

class MeshFileReader
{
public:
//	MeshFileReader();
	MeshFileReader(const char* pFileName);
	~MeshFileReader();

public:
	const Mesh& GetMesh() const;
	const MESH_TYPE& GetMeshType() const;
	void GetScalarFields();

private:
	void ReadOffFile();
	void ReadMeshFile();
	void ReadVtkFile();
	void ReadObjFile();
	void ReadStlFile();

private:
	std::string m_strFileName;
	MESH_TYPE m_meshType;
	Mesh m_mesh;
};

#endif // __Mesh_File_Reader_H__
