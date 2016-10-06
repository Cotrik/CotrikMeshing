/* Simple geometry viewer:  Left mouse: rotate;  Middle mouse:  zoom;  Right mouse:   menu;  ESC to quit
 The function InitGeometry() initializes  the geometry that will be displayed. */
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include "SimpleGeometryViewer.h"
#include "include.h"

#include "math3d.h"
#include "Selection.h"

//////////////////////////////////////
Selection selection;
M3DVector2f left_bottom, right_top;
bool bool_select_area = false;
int viewport[4];
//////////////////////////////////////
/* Viewer state */
float sphi = 90.0, stheta = 45.0;
float sdepth = 2;
float zNear = 0.01, zFar = 100.0;
float aspect = 5.0 / 4.0;
float xcam = 0, ycam = 0;
long xsize, ysize;
int downX, downY;
float moveX, moveY;
bool leftButton = false, middleButton = false, wheelUp = false, wheelDown = false;
GLfloat light0Position[] =
{ 0, 1, 0, 1.0 };
int displayMenu, meshMenu, mainMenu;

bool bKeyPressed[256] = {0};
enum
{
	WIREFRAME,
	HIDDENLINE,
	FLATSHADED,
	SMOOTHSHADED,
};

enum
{
	Origin_Tet_Mesh,
	Origin_Tet_Mesh_Surface,
	Polycube_Tet_Mesh,
	Polycube_Tet_Mesh_Surface,
	Polycube_Hex_Mesh,
	Magnified_Tet_Mesh,
	Magnified_Hex_Mesh,
	Deformed_Hex_Mesh,
	Result_Tet_Mesh
};

int displayMode = WIREFRAME;

#define GLUT_SCROLL_UP_BUTTON   0x00000003
#define GLUT_SCROLL_DOWN_BUTTON 0x00000004

const float htmlcolor[16][3] =
{
	{1.0,  1.0,  1.0 }, // White
	{0.75, 0.75, 0.75}, // Silver
	{0.5,  0.5,  0.5},  // Gray
	{0.0,  0.0,  0.0},  // Black
	{1.0,  0.0,  0.0},  // Red
	{0.5,  0.0,  0.0},  // Maroon
	{1.0,  1.0,  0.0},  // Yellow
	{0.5,  0.5,  0.0},  // Olive
	{0.0,  1.0,  0.0},  // Lime
	{0.0,  0.5,  0.0},  // Green
	{0.0,  1.0,  1.0},  // Aqua
	{0.0,  0.5,  0.5},  // Teal
	{0.0,  0.0,  1.0},  // Blue
	{0.0,  0.0,  0.5},  // Navy
	{1.0,  0.0,  1.0},  // Fuchsia
	{0.5,  0.0,  0.5}   // Purple
};
////////////////////////////////////////////////
void DoNothing(void)
{

}
typedef void(*KeyboardFunc)(void);
// Keyboard Function()
struct Key_Func
{
	unsigned char key;
	KeyboardFunc func;
};

const Key_Func keyFuncs[] =
{
	{'0', DoNothing},
	{'1', DoNothing},
	{'2', DoNothing},
	{'3', DoNothing},
	{'4', DoNothing},
	{'5', DoNothing},
	{'6', DoNothing},
	{'7', DoNothing},
	{'8', DoNothing},
	{'9', DoNothing},
	{'a', CurrentFacePatchIndex_Subtract},
	{'b', DisplayBoundaryVertices},
	{'c', DoNothing},
	{'d', CurrentFacePatchIndex_Plus},
	{'e', DoNothing},
	{'f', DisplayFacePatch},
	{'g', DoNothing},
	{'h', DoNothing},
	{'i', RelabelIsolatedPatch},
	{'j', CurrentVertexPatchIndex_Subtract},
	{'k', CurrentVertexPatchIndex_Plus},
	{'l', DoNothing},
	{'m', MagnifyMesh},
	{'n', DoNothing},
	{'o', DisplayOverlappingVertices},
	{'p', DoNothing},
	{'q', DisplayCorners},
	{'r', DoNothing},
	{'s', DesplaySelectedVertex},
	{'t', DoNothing},
	{'u', DoNothing},
	{'v', DisplayVertexPatch},
	{'w', RelabelWedgePatch},
	{'x', DoNothing},
	{'y', DoNothing},
	{'z', DoNothing},
	{'A', DoNothing},
	{'B', DoNothing},
	{'C', DoNothing},
	{'D', DoNothing},
	{'E', DoNothing},
	{'F', DoNothing},
	{'G', DoNothing},
	{'H', DoNothing},
	{'I', DoNothing},
	{'J', DoNothing},
	{'K', DoNothing},
	{'L', DoNothing},
	{'M', DoNothing},
	{'N', DoNothing},
	{'O', DoNothing},
	{'P', DoNothing},
	{'Q', DoNothing},
	{'R', DoNothing},
	{'S', DoNothing},
	{'T', DoNothing},
	{'U', DoNothing},
	{'V', DoNothing},
	{'W', DoNothing},
	{'X', DoNothing},
	{'Y', DoNothing},
	{'Z', DoNothing},
	{'\r', SavePolycubeTetMesh},
	{' ', RelabelSurfaceFacePerFace},
	{27, Exit}
};

void CallKeyFunc(unsigned char key)
{
	const int n = sizeof(keyFuncs) / sizeof(keyFuncs[0]);
	for (int i = 0; i < n; i++)
	{
		if (keyFuncs[i].key == key)
		{
			keyFuncs[i].func();
			break;
		}
	}
}
////////////////////////////////////////////////
void MyIdleFunc(void)
{
	glutPostRedisplay();
} /* things to do while idle */
void RunIdleFunc(void)
{
	glutIdleFunc(MyIdleFunc);
}
void PauseIdleFunc(void)
{
	glutIdleFunc(NULL);
}

void DrawSmoothShaded(void)
{
	const size_t cellSize = pMesh->C.size();
	for (unsigned int i = 0; i < cellSize; i++)
	{
		DrawCell(pMesh->C.at(i));
	}
}

void DrawWireframe(void)
{
	const size_t cellSize = pMesh->C.size();
	for (unsigned int i = 0; i < cellSize; i++)
	{
		DrawWireframe(pMesh->C.at(i));
	}
}

void DrawFlatShaded(void)
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	const size_t cellSize = pMesh->C.size();
	for (unsigned int i = 0; i < cellSize; i++)
	{
		DrawCell(pMesh->C.at(i));
	}
	if (pMesh == pHexMesh || pMesh == pHexPolycubeMesh)
	{
		const std::vector<Vertex>& v = pMesh->V;
		const Cell& c = pMesh->C.at(SELECTED_CUBE_INDEX);
		DrawHexahedron_S(v[c[0]], v[c[1]], v[c[2]], v[c[3]], v[c[4]], v[c[5]], v[c[6]], v[c[7]]);
	}
	glDisable(GL_POLYGON_OFFSET_FILL);
}

void DrawHiddenLine(void)
{
	glEnable(GL_POLYGON_OFFSET_FILL);
	glColor3f(0, 0, 0);
	glBegin( GL_TRIANGLES);

	glEnd();
	glDisable(GL_POLYGON_OFFSET_FILL);
	glColor3f(1.0, 1.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void ReshapeCallback(int width, int height)
{
	xsize = width;
	ysize = height;
	aspect = (float) xsize / (float) ysize;
	glViewport(0, 0, xsize, ysize);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glutPostRedisplay();
}

void SetDisplayMenu(int value)
{
	displayMode = value;
	switch (value)
	{
	case WIREFRAME:
		glShadeModel(GL_FLAT);
		glDisable(GL_LIGHTING);
		break;
	case HIDDENLINE:
		glShadeModel(GL_FLAT);
		glDisable(GL_LIGHTING);
		break;
	case FLATSHADED:
		glShadeModel(GL_FLAT);
		glDisable(GL_LIGHTING);
		break;
	case SMOOTHSHADED:
		glShadeModel(GL_SMOOTH);
		glEnable(GL_LIGHTING);
		break;
	}
	glutPostRedisplay();
}

void SetMeshSelectionMenu(int value)
{
	switch (value)
	{
	case Origin_Tet_Mesh:
		if (pTetMesh != NULL) pMesh = pTetMesh;
		break;
	case Origin_Tet_Mesh_Surface:
		if (pTetSurfaceMesh != NULL) pMesh = pTetSurfaceMesh;
		break;
	case Polycube_Tet_Mesh:
		if (pTetPolycubeMesh != NULL) pMesh = pTetPolycubeMesh;
		break;
	case Polycube_Tet_Mesh_Surface:
		if (pTetPolycubeSurfaceMesh != NULL) pMesh = pTetPolycubeSurfaceMesh;
		break;
	case Polycube_Hex_Mesh:
		if (pHexMesh != NULL) pMesh = pHexMesh;
		break;
	case Magnified_Tet_Mesh:
		if (pTetMagnifiedMesh != NULL) pMesh = pTetMagnifiedMesh;
		break;
	case Magnified_Hex_Mesh:
		if (pHexMagnifiedMesh != NULL) pMesh = pHexMagnifiedMesh;
		break;
	case Deformed_Hex_Mesh:
		if (pHexMesh != NULL) pMesh = pHexMesh;
		break;
	case Result_Tet_Mesh:
		if (pHexMesh != NULL) pMesh = pResultMesh;
		break;
	}
	glutPostRedisplay();
}

void SetMainMenu(int value)
{
	switch (value)
	{
	case 99:
		exit(0);
		break;
	}
}
void DrawAxes()
{
	glColor3f(0.0f, 0.0f, 1.0f); //指定线的颜色,蓝色

	glBegin( GL_LINES);
	{
		// x-axis
		glColor3f(1.0f, 0.0f, 0.0f); //指定线的颜色,蓝色
		glVertex3f(-1.0f, 0.0f, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);
		// x-axis arrow
		glVertex3f(1.0f, 0.0f, 0.0f);
		glVertex3f(0.930f, 0.030f, 0.0f);
		glVertex3f(1.0f, 0.0f, 0.0f);
		glVertex3f(0.930f, -0.03f, 0.0f);

		// y-axis
		glColor3f(0.0f, 1.0f, 0.0f);
		glVertex3f(0.0f, -1.0f, 0.0f);
		glVertex3f(0.0f, 1.0f, 0.0f);
		glVertex3f(0.0f, 1.0f, 0.0f);
		glVertex3f(0.03f, 0.90f, 0.0f);
		glVertex3f(0.0f, 1.0f, 0.0f);
		glVertex3f(-0.03f, 0.930f, 0.0f);

		// z-axis
		glColor3f(1.0f, 1.0f, 0.0f);
		glVertex3f(0.0f, 0.0f, -1.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);
		glVertex3f(0.0f, 0.0f, 1.0f);
		glVertex3f(0.030f, 0.0f, 0.930f);
		glVertex3f(0.0f, 0.0f, 1.0f);
		glVertex3f(-0.030f, 0.0f, 0.930f);
	}
	glEnd();
}
static void showedge()
{
	if (pMesh != NULL && pMesh == pTetPolycubeMesh && bKeyPressed['p'])
	{
		glm::vec3 colorx(1.0f, 0.0f, 0.0f);
		glm::vec3 colory(0.0f, 1.0f, 0.0f);
		glm::vec3 colorz(1.0f, 1.0f, 0.0f);
		const std::vector<Vertex>& V = pTetPolycubeMesh-> V;
		if (bKeyPressed['x'])
		for (int i = 0; i < pTetPolycubeMesh->m_vecXBoundaryLine.size(); i++)
		{
			const Line& l = pTetPolycubeMesh->m_vecXBoundaryLine.at(i);
			DrawLine(V[l.p0], V[l.p1], colorx);
		}
		if (bKeyPressed['y'])
		for (int i = 0; i < pTetPolycubeMesh->m_vecYBoundaryLine.size(); i++)
		{
			const Line& l = pTetPolycubeMesh->m_vecYBoundaryLine.at(i);
			DrawLine(V[l.p0], V[l.p1], colory);
		}
		if (bKeyPressed['z'])
		for (int i = 0; i < pTetPolycubeMesh->m_vecZBoundaryLine.size(); i++)
		{
			const Line& l = pTetPolycubeMesh->m_vecZBoundaryLine.at(i);
			DrawLine(V[l.p0], V[l.p1], colorz);
		}
	}
	if (pMesh != NULL && pMesh == pTetMesh && bKeyPressed['p'])
	{
		glm::vec3 colorx(1.0f, 0.0f, 0.0f);
		glm::vec3 colory(0.0f, 1.0f, 0.0f);
		glm::vec3 colorz(1.0f, 1.0f, 0.0f);
		const std::vector<Vertex>& V = pTetMesh-> V;
		if (bKeyPressed['x'])
		for (int i = 0; i < pTetPolycubeMesh->m_vecXBoundaryLine.size(); i++)
		{
			const Line& l = pTetPolycubeMesh->m_vecXBoundaryLine.at(i);
			DrawLine(V[l.p0], V[l.p1], colorx);
		}
		if (bKeyPressed['y'])
		for (int i = 0; i < pTetPolycubeMesh->m_vecYBoundaryLine.size(); i++)
		{
			const Line& l = pTetPolycubeMesh->m_vecYBoundaryLine.at(i);
			DrawLine(V[l.p0], V[l.p1], colory);
		}
		if (bKeyPressed['z'])
		for (int i = 0; i < pTetPolycubeMesh->m_vecZBoundaryLine.size(); i++)
		{
			const Line& l = pTetPolycubeMesh->m_vecZBoundaryLine.at(i);
			DrawLine(V[l.p0], V[l.p1], colorz);
		}
	}
}

bool hasReadHightLightLabel = false;
int current_face_patch_index = 0;
int current_vertex_patch_index = 0;

void DisplayXFacePatch(void)
{
	if (pMesh != NULL && (pMesh == pTriMesh || pMesh == pTriPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pTriPolycubeMesh->faceType.at(i) == FACE_X)
			{
				glm::vec3 color(
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][2]);
				DrawTriangle(V[f.at(0)], V[f.at(1)], V[f.at(2)], color);
			}
		}
	}

	if (pMesh != NULL && (pMesh == pTetMesh || pMesh == pTetPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pTetPolycubeMesh->faceType.at(i) == FACE_X)
			{
				glm::vec3 color(
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][2]);
				DrawTriangle(V[f.at(0)], V[f.at(1)], V[f.at(2)], color);
			}
		}
	}

	if (pMesh != NULL && (pMesh == pHexMesh || pMesh == pHexPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pHexPolycubeMesh->faceType.at(i) == FACE_X)
			{
				glm::vec3 color(
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][2]);
				DrawQuad(V[f.at(0)], V[f.at(1)], V[f.at(2)], V[f.at(3)], color);
			}
		}
	}
}

void DisplayYFacePatch(void)
{
	if (pMesh != NULL && (pMesh == pTriMesh || pMesh == pTriPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pTriPolycubeMesh->faceType.at(i) == FACE_Y)
			{
				glm::vec3 color(
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][2]);
				DrawTriangle(V[f.at(0)], V[f.at(1)], V[f.at(2)], color);
			}
		}
	}

	if (pMesh != NULL && (pMesh == pTetMesh || pMesh == pTetPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pTetPolycubeMesh->faceType.at(i) == FACE_Y)
			{
				glm::vec3 color(
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][2]);
				DrawTriangle(V[f.at(0)], V[f.at(1)], V[f.at(2)], color);
			}
		}
	}

	if (pMesh != NULL && (pMesh == pHexMesh || pMesh == pHexPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pHexPolycubeMesh->faceType.at(i) == FACE_Y)
			{
				glm::vec3 color(
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][2]);
				DrawQuad(V[f.at(0)], V[f.at(1)], V[f.at(2)], V[f.at(3)], color);
			}
		}
	}
}

void DisplayZFacePatch(void)
{
	if (pMesh != NULL && (pMesh == pTriMesh || pMesh == pTriPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pTriPolycubeMesh->faceType.at(i) == FACE_Z)
			{
				glm::vec3 color(
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pTriPolycubeMesh->faceLabel[i] % 16][2]);
				DrawTriangle(V[f.at(0)], V[f.at(1)], V[f.at(2)], color);
			}
		}
	}

	if (pMesh != NULL && (pMesh == pTetMesh || pMesh == pTetPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pTetPolycubeMesh->faceType.at(i) == FACE_Z)
			{
				glm::vec3 color(
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pTetPolycubeMesh->faceLabel[i] % 16][2]);
				DrawTriangle(V[f.at(0)], V[f.at(1)], V[f.at(2)], color);
			}
		}
	}

	if (pMesh != NULL && (pMesh == pHexMesh || pMesh == pHexPolycubeMesh))
	{
		const std::vector<Vertex>& V = pMesh->V;
		for (int i = 0; i < pMesh->surface.size(); i++)
		{
			const Face& f = pMesh->surface.at(i);
			if (pHexPolycubeMesh->faceType.at(i) == FACE_Y)
			{
				glm::vec3 color(
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][0],
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][1],
						htmlcolor[pHexPolycubeMesh->faceLabel[i] % 16][2]);
				DrawQuad(V[f.at(0)], V[f.at(1)], V[f.at(2)], V[f.at(3)], color);
			}
		}
	}
}

void DisplayFacePatch(void)
{
	if ((pMesh != NULL && pMesh == pTetMesh) || (pMesh != NULL && pMesh == pHexMesh) || (pMesh != NULL && pMesh == pTriMesh))
	{
		const std::vector<Vertex>& V = pMesh-> V;
		const std::vector<std::vector<unsigned long> >& facePatches = pMesh->facePatches;
		glm::vec3 color(1.0f, 0.0f, 0.0f);
		if (!pMesh->facePatches.empty())
		{
			glm::vec3 color(1.0f, 0.0f, 0.0f);
			const std::vector<unsigned long>& facePatch = facePatches.at(current_face_patch_index);
			for (int j = 0; j < facePatch.size(); j++)
			{
				unsigned long face_index = facePatches.at(current_face_patch_index).at(j);
				const Face& f = pMesh->surface.at(face_index);
				if (f.size() == 3)
				DrawTriangle(V[f.at(0)], V[f.at(1)], V[f.at(2)], color);
				else if (f.size() == 4)
				DrawQuad(V[f.at(0)], V[f.at(1)], V[f.at(2)], V[f.at(3)], color);
			}
		}
	}
}

static void DisplayHighLightFacePatch()
{
	if (pMesh != NULL && pMesh == pTetMesh && bKeyPressed['h'])
	{
		glm::vec3 color(1.0f, 0.0f, 0.0f);
		const std::vector<Vertex>& V = pTetMesh->V;
		const std::vector<std::vector<unsigned long> >& facePatches = pTetMesh->facePatches;
		for (int i = 0; i < pTetMesh->highlightPatchLabel.size(); i++)
		{
			const std::vector<unsigned long>& facePatch = facePatches.at(i);
			for (int j = 0; j < facePatch.size(); j++)
			{
				unsigned long face_index = facePatches.at(i).at(j);
				const Face& f = pTetMesh->surface.at(face_index);
				DrawTriangle(V[f.at(0)], V[f.at(1)], V[f.at(2)], color);
			}
		}
	}
}
static void showface()
{
	if (bKeyPressed['p'])
	{
		if (bKeyPressed['x']) DisplayXFacePatch();
		if (bKeyPressed['y']) DisplayYFacePatch();
		if (bKeyPressed['z']) DisplayZFacePatch();
	}
	if (bKeyPressed['f']) DisplayFacePatch();
	DisplayHighLightFacePatch();
}

std::vector<unsigned long> errorPointIndices;
void onInitialization( )
{
    glEnable( GL_POINT_SPRITE ); // GL_POINT_SPRITE_ARB if you're
                                 // using the functionality as an extension.

    glEnable( GL_POINT_SMOOTH );
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
    glPointSize( 10.0 );

    /* assuming you have setup a 32-bit RGBA texture with a legal name */
//    glActiveTexture(GL_TEXTURE0);
//    glEnable( GL_TEXTURE_2D );
//    glTexEnv(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
//    glTexEnv(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
//    glBindTexture(GL_TEXTURE_2D, texture_name);
}
bool hasReadSelectedPoints = false;

void DisplayCorners(void)
{
	if (pMesh == NULL)
	{
		return;
	}
	if (pMesh == pHexPolycubeMesh || pMesh == pTetPolycubeMesh)
	{
		glColor4f( 0.00f, 0.00f, 0.931f, 1.0f );
		const std::vector<Vertex>& v = pMesh->V;
		for (int i = 0; i < ((PolycubeMesh*)pMesh)->m_vec_xyPlaneCorner.size(); i++)
		{
			const unsigned long index =	((PolycubeMesh*)pMesh)->m_vec_xyPlaneCorner.at(i);
			glVertex3f(v[index].x, v[index].y, v[index].z);
		}
	}
}
void DisplayOverlappingVertices(void)
{
	if (pMesh == pHexPolycubeMesh || pMesh == pTetPolycubeMesh)
	{
		glColor4f( 1.00f, 1.00f, 0.50f, 1.0f );
		const std::vector<Vertex>& v = pMesh->V;
		for (int i = 0; i < ((PolycubeMesh*)pMesh)->overlappingVertices.size(); i++)
		{
			const unsigned long index =	((PolycubeMesh*)pMesh)->overlappingVertices.at(i);
			glVertex3f(v[index].x, v[index].y, v[index].z);
		}
	}
}
void DisplayBoundaryVertices(void)
{
	if (pMesh == NULL)
	{
		return;
	}
	const std::vector<Vertex>& v = pMesh->V;
	if (bKeyPressed['x'])
	for (int i = 0; i < ((PolycubeMesh*)pMesh)->m_vecXBoundaryLine.size(); i++)
	{
		const Line& line = ((PolycubeMesh*)pMesh)->m_vecXBoundaryLine[i];
		unsigned long index = line.p0;
		glVertex3f(v[index].x, v[index].y, v[index].z);
		index =	line.p1;
		glVertex3f(v[index].x, v[index].y, v[index].z);
	}
	if (bKeyPressed['y'])
	for (int i = 0; i < ((PolycubeMesh*)pMesh)->m_vecYBoundaryLine.size(); i++)
	{
		const Line& line = ((PolycubeMesh*)pMesh)->m_vecYBoundaryLine[i];
		unsigned long index =	line.p0;
		glVertex3f(v[index].x, v[index].y, v[index].z);
		index =	line.p1;
		glVertex3f(v[index].x, v[index].y, v[index].z);
	}
	if (bKeyPressed['z'])
	for (int i = 0; i < ((PolycubeMesh*)pMesh)->m_vecZBoundaryLine.size(); i++)
	{
		const Line& line = ((PolycubeMesh*)pMesh)->m_vecZBoundaryLine[i];
		unsigned long index =	line.p0;
		glVertex3f(v[index].x, v[index].y, v[index].z);
		index =	line.p1;
		glVertex3f(v[index].x, v[index].y, v[index].z);
	}
}
void DisplayVertexPatch(void)
{
	if (pMesh == NULL)
	{
		return;
	}
	const std::vector<Vertex>& v = pMesh->V;
	int count = 0;
	for (int i = 0; i < pMesh->m_vertex_patches.size(); i++)
	{
		for (int j = 0; j < pMesh->m_vertex_patches.at(i).size(); j++)
		{
			if (count++ == current_vertex_patch_index)
			{
				for (int k = 0;	k < pMesh->m_vertex_patches.at(i).at(j).size();	k++)
				{
					const unsigned long index =	pMesh->m_vertex_patches.at(i).at(j).at(k);
					glVertex3f(v[index].x, v[index].y, v[index].z);
				}
			}
		}
	}
}

void DesplaySelectedVertex(void)
{
    if (pMesh != NULL)
    {
    	if (!hasReadSelectedPoints)
    	{
			ifstream file("selected_points.txt");
			unsigned long index = 0;
			while (file >> index)
			{
				pMesh->outsidePoints.push_back(index);
			}
			hasReadSelectedPoints = true;
    	}
		for ( int i = 0; i < pMesh->outsidePoints.size(); ++i )
		{
			const std::vector<Vertex>& v = pMesh->V;
			const unsigned long index = pMesh->outsidePoints.at(i);
			glVertex3f(v[index].x, v[index].y, v[index].z);
		}
    }
}
void onDisplay()
{
//    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
//    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    glBegin( GL_POINTS );
    //glColor4f( 0.95f, 0.207, 0.031f, 1.0f );
    glColor4f( 0.00f, 0.93f, 0.031f, 1.0f );

    if (pHexMagnifiedMesh != NULL)
    {
		for ( int i = 0; i < errorPointIndices.size(); ++i )
		{
			const std::vector<Vertex>& v = pHexMagnifiedMesh->V;
			const unsigned long index = errorPointIndices.at(i);
			glVertex3f(v[index].x, v[index].y, v[index].z);
		}
    }

    if (bKeyPressed['s'])   DesplaySelectedVertex();
	if (bKeyPressed['v'])   DisplayVertexPatch();
	if (bKeyPressed['q'])	DisplayCorners();
	if (bKeyPressed['o'])	DisplayOverlappingVertices();
	if (bKeyPressed['b'])   DisplayBoundaryVertices();
    glEnd();
    glFinish();
//    glutSwapBuffers();
}

void ShowHexMesh_WIREFRAME(void)
{
	const Cell& c = pMesh->C.at(SELECTED_CUBE_INDEX);
	const std::vector<Vertex>& v = pMesh->V;
	const std::vector<VertexInfo>& vInfo = pMesh->vertexInfo;
	if (pMesh != NULL && pMesh == pHexMesh)
	{
		for (unsigned int i = 0; i < 12; i++)
		{
			DrawWireFrameLine_S(v[c[HexEdge[i][0]]], v[c[HexEdge[i][1]]]);
		}
		for (unsigned int i = 0; i < pMesh->lowQualityCellIndex.size(); i++)
		{
			const Cell& cell = pMesh->C.at(pMesh->lowQualityCellIndex.at(i));
			for (unsigned int i = 0; i < 12; i++)
			{
				//if (vInfo[cell[HexEdge[i][0]]].bSurface || vInfo[cell[HexEdge[i][1]]].bSurface)
				DrawWireFrameLine_Yellow(v[cell[HexEdge[i][0]]], v[cell[HexEdge[i][1]]]);
			}
		}
		for (unsigned int i = 0; i < pMesh->invertedCellIndex.size(); i++)
		{
			const Cell& cell = pMesh->C.at(pMesh->invertedCellIndex.at(i));
			for (unsigned int i = 0; i < 12; i++)
			{
				//if (vInfo[cell[HexEdge[i][0]]].bSurface || vInfo[cell[HexEdge[i][1]]].bSurface)
				DrawWireFrameLine_S(v[cell[HexEdge[i][0]]], v[cell[HexEdge[i][1]]]);
			}
		}
	}
}
void ShowHexMesh_SMOOTHSHADED(void)
{
	if (pHexMesh != NULL)
	{
		const Cell& c = pHexMesh->C.at(SELECTED_CUBE_INDEX);
		const std::vector<Vertex>& v = pHexMesh->V;
		for (unsigned int i = 0; i < 12; i++)
		{
			DrawWireFrameLine_S(v[c[HexEdge[i][0]]], v[c[HexEdge[i][1]]]);
		}
	}
}
void DisplayCallback(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(64.0, aspect, zNear, zFar);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(moveX, moveY, -sdepth);
	glRotatef(stheta, 1.0, 0.0, 0.0);
	glRotatef(-sphi, 0.0, 0.0, 1.0);
	DrawAxes();

	///////////////////////////////////////////
	// Selection
	M3DMatrix44f mat_proj, mat_modelview;
	int width = glutGet( GLUT_WINDOW_WIDTH);
	int height = glutGet(GLUT_WINDOW_HEIGHT);
	if (bKeyPressed['s'])
	{
		glViewport(0, 0, (GLsizei) width, (GLsizei) height);
		glGetIntegerv(GL_VIEWPORT, viewport);
		glClear(GL_COLOR_BUFFER_BIT);

		glPushAttrib(GL_POLYGON_BIT);
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(10, (GLfloat) width / (GLfloat) height, 1, 300);

		glGetFloatv(GL_PROJECTION_MATRIX, mat_proj);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(0, 0, 3, 0, 0, 0, 0, 1, 0);
		glColor4f(0.1, 0.4, 0.6, 0.7);
		glPushMatrix();
	}
	///////////////////////////////////////////
	switch (displayMode)
	{
	case WIREFRAME:
	{
		//glShadeModel(GL_FLAT);
		onInitialization();
		onDisplay();
		DrawWireframe();
		ShowHexMesh_WIREFRAME();
		showedge();
		showface();
	}
		break;
	case HIDDENLINE:
		DrawHiddenLine();
		break;
	case FLATSHADED:
		glShadeModel (GL_FLAT);
		DrawFlatShaded();
		//showedge();
		showface();
		break;
	case SMOOTHSHADED:
	{
		glShadeModel (GL_SMOOTH);
		DrawSmoothShaded();
		ShowHexMesh_SMOOTHSHADED();
		showedge();
	}
		break;
	}
	///////////////////////////////////////////
	// Selection
	if (bKeyPressed['s'])
	{
		glPopMatrix();
		glPopAttrib();
		selection.set_config(*pMesh, pMesh->V.size(), left_bottom, right_top, mat_modelview, mat_proj, viewport);
		if (bool_select_area)
		{
			selection.draw_area();
			selection.highlight_selected_pts();
		}
	}
	///////////////////////////////////////////
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}
bool bEnterPressed = false;

void MapbackFromDeformedCubes_S(Mesh& mesh, const std::vector<unsigned long>& vec_cubeIndex)
{
	//parameterization of new OrigTet to new Cube Mesh
	std::cout << "parameterization of new OrigTet to new Cube Mesh" << std::endl;
	for (int i = 0; i < mesh.V.size()/*vertices.size()*/; i++)
	{
		if (!paras_flag[i])
			continue;
		Vertex& v = mesh.V.at(i/*vertices.at(i)*/);
		// 1. judge v is in which cube
		const unsigned long cubeIndex = vec_cubeIndex.at(i);
		const Cell& cubeCell = mesh.m_vecPolycubeHexCell.at(cubeIndex);
		const Cell& tet = tets.at(i);

		const Vertex& origT0 = mesh.m_vecPolycubeHexVertex[tet.at(0)];
		const Vertex& origT1 = mesh.m_vecPolycubeHexVertex[tet.at(1)];
		const Vertex& origT2 = mesh.m_vecPolycubeHexVertex[tet.at(2)];
		const Vertex& origT3 = mesh.m_vecPolycubeHexVertex[tet.at(3)];

		glm::vec3 p03((origT0.x - origT3.x), (origT0.y - origT3.y), (origT0.z - origT3.z));
		glm::vec3 p13((origT1.x - origT3.x), (origT1.y - origT3.y), (origT1.z - origT3.z));
		glm::vec3 p23((origT2.x - origT3.x), (origT2.y - origT3.y), (origT2.z - origT3.z));
		glm::mat3x3 OrigT(p03, p13, p23);
		glm::vec3 new_r4(origT3.x, origT3.y, origT3.z);
		const glm::vec3& rambda = paras.at(i);
		glm::vec3 new_r = OrigT * rambda;
		new_r += new_r4;
		Vertex baryCenter(new_r.x, new_r.y, new_r.z);
		v.x = baryCenter.x; v.y = baryCenter.y; v.z = baryCenter.z;
		//std::cout << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << " " << rambda.x << " " << rambda.y << " " << rambda.z << std::endl;
	}
}
void MagnifyMesh()
{
	bEnterPressed = true;
	//exit(1);
	pTetMesh->NormalizeCoordinateValue();
	//	origTetMesh.GenerateHexahedralMesh_Magnifier(origTetMesh, strOutputPolyHexFileName.c_str(),
	//			strOutputResultHexFileName.c_str(), cubesize, dir);
	pTetMesh->m_vecPolycubeHexCell.clear();
	pTetMesh->m_vecPolycubeHexVertex.clear();
	pTetMesh->GenerateSmallCubes_Magifier(CUBESIZE);
	//std::vector<unsigned long> vec_cubeIndex;
	pTetMesh->GetParameterAndCubeIndicesOfVertexIndices();

	FocusContextMagnifier focusContextMagnifier(*pTetMesh);
	focusContextMagnifier.SolveLinearEquation();
	//pTetMesh->WriteHexahedralmesh(strMagnifiedHexCpMeshFileName.c_str());

	pTetMesh->MapbackFromDeformedCubes();
	MeshFileWriter magnifiedTetMeshFileWriter(*pTetMesh);
	//magnifiedTetMeshFileWriter.WriteMeshFile(strMagnifiedTetMeshFileName.c_str(), MESH_TYPE_TETRAHEDRON_VTK);
	static int count = 1;
	char filename[128] = { 0 };
	std::string str;
	sprintf(filename, "magnified.iter%d.tet.vtk", count++);
	magnifiedTetMeshFileWriter.WriteMeshFile(filename, MESH_TYPE_TETRAHEDRON_VTK);
}
void SavePolycubeTetMesh()
{
	//pTetPolycubeMesh->SmoothBoundaryLine(43, 1);
	pTetPolycubeMesh->SmoothBoundaryLine();
	pTetPolycubeMesh->SmoothSurface();
	pTetPolycubeMesh->SmoothVolume();
	MeshFileWriter meshFileWriter(*pTetPolycubeMesh);
	meshFileWriter.WriteMeshFile("polycube.cleanup.tet.vtk", MESH_TYPE_TETRAHEDRON_VTK);
	std::cout << "successfully write polycube.cleanup.tet.vtk" << std::endl;
}
void CurrentFacePatchIndex_Subtract(void)
{
	current_face_patch_index--;
	if (current_face_patch_index < 0)	current_face_patch_index = 0;
	std::cout << "current_face_patch_index = " << current_face_patch_index	<< std::endl;
}

void CurrentFacePatchIndex_Plus(void)
{
	current_face_patch_index++;
	if (current_face_patch_index >= pMesh->m_patchNumber) current_face_patch_index = pMesh->m_patchNumber - 1;
	std::cout << "current_face_patch_index = " << current_face_patch_index << std::endl;
}

void CurrentVertexPatchIndex_Subtract(void)
{
	current_vertex_patch_index--;
	if (current_vertex_patch_index < 0) current_vertex_patch_index = 0;
	std::cout << "current_vertex_patch_index = " << current_vertex_patch_index << std::endl;
}

void CurrentVertexPatchIndex_Plus(void)
{
	current_vertex_patch_index++;
	if (current_vertex_patch_index >= pMesh->m_vertex_patches_num) current_vertex_patch_index = pMesh->m_vertex_patches_num - 1;
	std::cout << "current_vertex_patch_index = " << current_vertex_patch_index << std::endl;
}

void RelabelWedgePatch()
{
	if (pMesh != NULL
		&& ((pMesh == pTetMesh && pTetPolycubeMesh != NULL) || pMesh == pTetPolycubeMesh))
	{
		pTetPolycubeMesh->RelabelWedgePatch();
		std::cout << "RelabelWedgePatch done!" << std::endl;
	}
}

void RelabelIsolatedPatch()
{
	if (pMesh != NULL
		&& ((pMesh == pTetMesh && pTetPolycubeMesh != NULL) || pMesh == pTetPolycubeMesh))
	{
		pTetPolycubeMesh->RelabelIsolatedPatch();
		std::cout << "RelabelIsolatedPatch done!" << std::endl;
	}
}

void RelabelSurfaceFacePerFace()
{
	// Polycube Mesh cleaning
	if (pMesh != NULL
		&& ((pMesh == pTetMesh && pTetPolycubeMesh != NULL) || pMesh == pTetPolycubeMesh))
	{
		const std::vector<Vertex>& PV = pTetPolycubeMesh->V;
		const std::vector<Face>& PF = pTetPolycubeMesh->surface;
		const std::vector<Vertex>& OV = pTetMesh->V;
		const std::vector<Face>& OF = pTetMesh->surface;
		std::vector<FACE_TYPE>& PFT = pTetPolycubeMesh->faceType;
		std::vector<FACE_TYPE>& OFT = pTetMesh->faceType;
		const size_t PF_SIZE = PF.size();
		pTetPolycubeMesh->RelabelSurfaceFacePerFace();
		OFT = PFT;
	}
}

void Exit(void)
{
	exit(1);
}
void KeyboardCallback(unsigned char key, int x, int y)
{
	//  Print what key the user is hitting
	printf("User is hitting the '%c' key. ASCII code is %d.\n", key, key);
	bKeyPressed[key] = !bKeyPressed[key];
	CallKeyFunc(key);
	glutPostRedisplay();
}

//-------------------------------------------------------------------------
//  This function is passed to the glutSpecialFunc and is called
//  whenever the user hits a special key.
//-------------------------------------------------------------------------
void SpecialKeyboardCallback(int key, int x, int y)
{
	switch (key)
	{
		case GLUT_KEY_F1 :
			printf ("F1 function key.\n");
			break;
		case GLUT_KEY_F2 :
			printf ("F2 function key. \n");
			break;
		case GLUT_KEY_F3 :
			printf ("F3 function key. \n");
			break;
		case GLUT_KEY_F4 :
			printf ("F4 function key. \n");
			break;
		case GLUT_KEY_F5 :
			printf ("F5 function key. \n");
			break;
		case GLUT_KEY_F6 :
			printf ("F6 function key. \n");
			break;
		case GLUT_KEY_F7 :
			printf ("F7 function key. \n");
			break;
		case GLUT_KEY_F8 :
			printf ("F8 function key. \n");
			break;
		case GLUT_KEY_F9 :
			printf ("F9 function key. \n");
			break;
		case GLUT_KEY_F10 :
			printf ("F10 function key. \n");
			break;
		case GLUT_KEY_F11 :
			printf ("F11 function key. \n");
			break;
		case GLUT_KEY_F12 :
			printf ("F12 function key. \n");
			break;
		case GLUT_KEY_LEFT :
			if (SELECTED_CUBE_INDEX > 13)
			SELECTED_CUBE_INDEX -= 13;
			SELECTED_CUBE_INDEX = SELECTED_CUBE_INDEX % 2197;
			printf ("Left directional key. \n");
			break;
		case GLUT_KEY_UP :
			SELECTED_CUBE_INDEX += 1;
			SELECTED_CUBE_INDEX = SELECTED_CUBE_INDEX % 2197;
			printf ("Up directional key. \nSELECTED_CUBE_INDEX = %d\n", SELECTED_CUBE_INDEX);
			break;
		case GLUT_KEY_RIGHT :
			SELECTED_CUBE_INDEX += 13;
			SELECTED_CUBE_INDEX = SELECTED_CUBE_INDEX % 2197;
			printf ("Right directional key. \nSELECTED_CUBE_INDEX = %d\n", SELECTED_CUBE_INDEX);
			break;
		case GLUT_KEY_DOWN :
			SELECTED_CUBE_INDEX -= 1;
			SELECTED_CUBE_INDEX = SELECTED_CUBE_INDEX % 2197;
			printf ("Down directional key. \nSELECTED_CUBE_INDEX = %d\n", SELECTED_CUBE_INDEX);
			break;
		case GLUT_KEY_PAGE_UP :
			SELECTED_CUBE_INDEX += 169;
			SELECTED_CUBE_INDEX = SELECTED_CUBE_INDEX % 2197;
			printf ("Page up directional key. \nSELECTED_CUBE_INDEX = %d\n", SELECTED_CUBE_INDEX);
			break;
		case GLUT_KEY_PAGE_DOWN :
			if (SELECTED_CUBE_INDEX > 169)
			SELECTED_CUBE_INDEX -= 169;
			SELECTED_CUBE_INDEX = SELECTED_CUBE_INDEX % 2197;
			printf ("Page down directional key. \nSELECTED_CUBE_INDEX = %d\n", SELECTED_CUBE_INDEX);
			break;
		case GLUT_KEY_HOME :
			printf ("Home directional key. \n");
			break;
		case GLUT_KEY_END :
			printf ("End directional key. \n");
			break;
		case GLUT_KEY_INSERT :
			printf ("Inset directional key. \n");
			break;
	}

	glutPostRedisplay ();
}

void MouseCallback(int button, int state, int x, int y)
{
	///////////////////////////
	// Selection
	if (bKeyPressed['s'])
	{
		int width = glutGet( GLUT_WINDOW_WIDTH);
		int height = glutGet(GLUT_WINDOW_HEIGHT);

		if (button == GLUT_LEFT_BUTTON)
		{
			if (state == GLUT_DOWN)
			{
				bool_select_area = true;
			}
			else if (state == GLUT_UP)
			{
				bool_select_area = false;
			}

			m3dLoadVector2(left_bottom, x, height - y);
			m3dLoadVector2(right_top, x, height - y);
		}

		glutPostRedisplay();
		return;
	}
	///////////////////////////
	switch (button)
	{
	case GLUT_LEFT_BUTTON:
		switch (state)
		{
		case GLUT_DOWN:
			//printf("Mouse Left Button Pressed (Down)...\n");
			break;
		case GLUT_UP:
			//printf("Mouse Left Button Released (Up)...\n");
			break;
		}
		break;
	case GLUT_MIDDLE_BUTTON:
		switch (state)
		{
		case GLUT_DOWN:
			//printf("Mouse Middle Button Pressed (Down)...\n");
			break;
		case GLUT_UP:
			//printf("Mouse Middle Button Released (Up)...\n");
			break;
		}
		break;
	case GLUT_SCROLL_UP_BUTTON:
		switch (state)
		{
		case GLUT_DOWN:
			//printf("Mouse Wheel Scoll Up...\n");
			wheelUp = true;
			break;
		case GLUT_UP:
			//printf("Mouse Wheel Released (Up)...\n");
			sdepth -= fabs(sdepth - zNear)/4;
			break;
		}
		break;
	case GLUT_SCROLL_DOWN_BUTTON:
		switch (state)
		{
		case GLUT_DOWN:
			//printf("Mouse Wheel Scoll Down...\n");
			wheelDown = true;
			break;
		case GLUT_UP:
			//printf("Mouse Wheel Released (Up)...\n");
			sdepth += fabs(sdepth - zNear)/4;
			break;
		}
		break;
	}
	downX = x;
	downY = y;
	leftButton = ((button == GLUT_LEFT_BUTTON) && (state == GLUT_DOWN));
	middleButton = ((button == GLUT_MIDDLE_BUTTON) && (state == GLUT_DOWN));
	glutPostRedisplay();
}

void MouseWheelCallBack(int button, int dir, int x, int y)
{
    if (dir > 0)
    {
        // Zoom in
    	printf("Zoom in...\n");
    	downX = 1.5*x;
    	downY = 1.5*y;
    }
    else
    {
        // Zoom out
    	downX = 0.75*x;
    	downY = 0.75*y;
    }
    middleButton = true;
    glutPostRedisplay();
}

void MotionCallback(int x, int y)
{
	if (bKeyPressed['s'])
	{
		int width = glutGet( GLUT_WINDOW_WIDTH ), height = glutGet( GLUT_WINDOW_HEIGHT );

		if( bool_select_area ){
			m3dLoadVector2(right_top, x, height - y);
		}

		glutPostRedisplay();
		return;
	}

	if (leftButton)
	{
		sphi += (float) (downX - x) / 4.0;
		stheta += (float) (downY - y) / 4.0;
	} // rotate
	if (middleButton)
	{
		moveX += ((float) (x - downX))/400.0;
		moveY += ((float) (downY - y))/400.0;
		//sdepth += (float) (downY - y) / 10.0;
	} // scale
	downX = x;
	downY = y;
	glutPostRedisplay();
}

void InitGL()
{
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500, 500);
	glutCreateWindow("Simple Geometry Viewer");
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glPolygonOffset(1.0, 1.0);
	glDisable(GL_CULL_FACE);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT, GL_DIFFUSE);
	glLightfv(GL_LIGHT0, GL_POSITION, light0Position);
	glEnable(GL_LIGHT0);
	glutReshapeFunc(ReshapeCallback);
	glutDisplayFunc(DisplayCallback);
	glutKeyboardFunc(KeyboardCallback);
	glutSpecialFunc(SpecialKeyboardCallback);
	glutMouseFunc(MouseCallback);
	//glutMouseWheelFunc(MouseWheelCallBack);
	glutMotionFunc(MotionCallback);
}

void InitMenu()
{
	displayMenu = glutCreateMenu(SetDisplayMenu);
	glutAddMenuEntry("Wireframe", WIREFRAME);
	glutAddMenuEntry("Hidden Line", HIDDENLINE);
	glutAddMenuEntry("Flat Shaded", FLATSHADED);
	glutAddMenuEntry("Smooth Shaded", SMOOTHSHADED);

	meshMenu = glutCreateMenu(SetMeshSelectionMenu);
	glutAddMenuEntry("Origin Tet Mesh",        Origin_Tet_Mesh);
	glutAddMenuEntry("Origin Tet Mesh_Surface",Origin_Tet_Mesh_Surface);
	glutAddMenuEntry("Polycube Tet Mesh",      Polycube_Tet_Mesh);
	glutAddMenuEntry("Polycube Tet Mesh_Surface",Polycube_Tet_Mesh_Surface);
	glutAddMenuEntry("Polycube Hex Mesh",      Polycube_Hex_Mesh);
	glutAddMenuEntry("Magnified Tet Mesh",      Magnified_Tet_Mesh);
	glutAddMenuEntry("Magnified Hex Mesh",      Magnified_Hex_Mesh);
	glutAddMenuEntry("Deformed Hex Mesh",      Deformed_Hex_Mesh);
	glutAddMenuEntry("Result Tet Mesh",        Result_Tet_Mesh);

	mainMenu = glutCreateMenu(SetMainMenu);
	glutAddSubMenu("Display", displayMenu);
	glutAddSubMenu("Select Mesh", meshMenu);
	glutAddMenuEntry("Exit", 99);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
}

void InitGeometry(const Mesh& mesh)
{
	pMesh = (Mesh*)&mesh;
	sdepth = mesh.m_maxVertex.z - mesh.m_minVertex.z;
	float x = mesh.m_maxVertex.x - mesh.m_minVertex.x;
	float y = mesh.m_maxVertex.y - mesh.m_minVertex.y;
	if (x > sdepth) sdepth = x;
	if (y > sdepth) sdepth = y;
	sdepth /= 1.2;
}

int DisplayMesh(int argc, char **argv, const Mesh& mesh)
{
	glutInit(&argc, argv);
	InitGL();
	InitMenu();
	InitGeometry(mesh);
	glutMainLoop();
	return 0;
}

