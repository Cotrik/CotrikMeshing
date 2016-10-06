#ifndef __SIMPLE_GEOMETRIY_VIEWER_H__
#define __SIMPLE_GEOMETRIY_VIEWER_H__
#include "include.h"

void MyIdleFunc(void);
void RunIdleFunc(void);
void PauseIdleFunc(void);
void DrawSmoothShaded(void);
void DrawWireframe(void);
void DrawFlatShaded(void);
void DrawHiddenLine(void);
void ReshapeCallback(int width, int height);

void SetDisplayMenu(int value);
void SetMeshSelectionMenu(int value);
void SetMainMenu(int value);

void KeyboardCallback(unsigned char key, int x, int y);
void SpecialKeyboardCallback(int key, int x, int y);
void MouseCallback(int button, int state, int x, int y);
void MouseWheelCallBack(int button, int dir, int x, int y);
void MotionCallback(int x, int y);

void InitGL();
void InitMenu();

void InitGeometry(const Mesh& mesh);

int DisplayMesh(int argc, char **argv, const Mesh& mesh);

/////////////////////////////////////////////////////////
void CurrentFacePatchIndex_Subtract(void);
void CurrentFacePatchIndex_Plus(void);
void CurrentVertexPatchIndex_Subtract(void);
void CurrentVertexPatchIndex_Plus(void);

void DisplayVertexPatch(void);
void DisplayFacePatch(void);
void DisplayCorners(void);
void DisplayOverlappingVertices(void);
void DisplayBoundaryVertices(void);
void RelabelWedgePatch(void);
void RelabelIsolatedPatch(void);
void SavePolycubeTetMesh(void);
void MagnifyMesh(void);
void RelabelSurfaceFacePerFace(void);
void DesplaySelectedVertex(void);
void Exit();
/////////////////////////////////////////////////////////
extern bool bEnterPressed;
extern bool bKeyPressed[256];
extern std::vector<unsigned long> errorPointIndices;
#endif // __SIMPLE_GEOMETRIY_VIEWER_H__
