#include "Selection.h"
#include <GL/glut.h>

Selection::Selection(void)
{
	this->vec_selected_pts_index.clear();
	this->pts = NULL;
}


Selection::~Selection(void)
{

}

void Selection::draw_area()
{
	int width = glutGet( GLUT_WINDOW_WIDTH ), height = glutGet( GLUT_WINDOW_HEIGHT );
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, width, 0, height);

		glMatrixMode( GL_MODELVIEW );
		glPushMatrix();
			glLoadIdentity();

			/*** ������ ***/
			glColor4f(1.0, 1.0, 0, 0.2);
			glEnable( GL_BLEND );
			glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glRectf( m3dGetVectorX(left_bottom), m3dGetVectorY(left_bottom), m3dGetVectorX(right_top), m3dGetVectorY(right_top));


			/*** ���������� ***/
			glEnable( GL_LINE_STIPPLE );
			glColor4f(1, 0, 0, 0.5);
			glLineStipple(3, 0xAAAA);

			glBegin( GL_LINE_LOOP );
			glVertex2f(m3dGetVectorX(left_bottom), m3dGetVectorY(left_bottom));
			glVertex2f(m3dGetVectorX(right_top), m3dGetVectorY(left_bottom));
			glVertex2f(m3dGetVectorX(right_top), m3dGetVectorY(right_top));
			glVertex2f(m3dGetVectorX(left_bottom), m3dGetVectorY(right_top));
			glEnd();

		glPopMatrix();
	glPopMatrix();
	glDisable( GL_LINE_STIPPLE );
	glDisable( GL_BLEND );
}

void Selection::highlight_selected_pts()
{
	int i = 0;

	//������һ��opengl״̬
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
		glLoadMatrixf(proj);

		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
			glLoadMatrixf(model_view);

			glPushAttrib(GL_POINT_BIT);
			glPointSize(3);
			glColor4f(0, 1, 0, 1);

			glBegin(GL_POINTS);

			for(i = 0; i != vec_selected_pts_index.size(); i++){
				glVertex3fv(pts[ vec_selected_pts_index[i] ]);
			}

			glEnd();
			glPopAttrib();
		glPopMatrix();
	glPopMatrix();
}

void Selection::get_selected_pts_index( vector<int> &v )
{
	// assign ��ֵ���
	v.assign(vec_selected_pts_index.begin(), vec_selected_pts_index.end());
}

void Selection::set_config( M3DVector3f *_pts, int _nr, M3DVector2f _left_bottom, M3DVector2f _right_top, M3DMatrix44f _model_view, M3DMatrix44f _proj, int _viewport[] )
{
	pts = _pts;
	nr = _nr;

	m3dCopyVector2(left_bottom, _left_bottom);
	m3dCopyVector2(right_top, _right_top);

	m3dCopyMatrix44(model_view, _model_view);
	m3dCopyMatrix44(proj, _proj);
	memcpy(viewport, _viewport, sizeof(int) * 4);

	/*** ���㱻ѡ�е�index ***/

	cal_selected_index();
}

Mesh* pSelectionMesh = NULL;
M3DVector3f *_pts = NULL;

void Selection::set_config( const Mesh& mesh, int _nr, M3DVector2f _left_bottom, M3DVector2f _right_top, M3DMatrix44f _model_view, M3DMatrix44f _proj, int _viewport[] )
{
	if (pSelectionMesh != &mesh)
	{
		pSelectionMesh = (Mesh*)&mesh;
		if (_pts != NULL)
		{
			delete [] _pts;
		}
		_pts = new M3DVector3f[mesh.V.size()];
		for (int i = 0; i < mesh.V.size(); i++)
		{
			_pts[i][0] = mesh.V.at(i).x;
			_pts[i][1] = mesh.V.at(i).y;
			_pts[i][2] = mesh.V.at(i).z;
		}
		std::cout << "Enter point copy!" <<std::endl;
	}

	pts = _pts;
	nr = _nr;

	m3dCopyVector2(left_bottom, _left_bottom);
	m3dCopyVector2(right_top, _right_top);

	m3dCopyMatrix44(model_view, _model_view);
	m3dCopyMatrix44(proj, _proj);
	memcpy(viewport, _viewport, sizeof(int) * 4);

	cal_selected_index();
}

void Selection::cal_selected_index()
{
	vec_selected_pts_index.clear();

	int i = 0;
	for(i = 0; i != nr; i++){

		if( drop_in_area(pts[i]) ){
			vec_selected_pts_index.push_back(i);
		}

	}
}

bool Selection::drop_in_area( M3DVector3f x )
{
	M3DVector2f win_coord;

	m3dProjectXY(win_coord, model_view, proj, viewport, x);

	if( (win_coord[0] < left_bottom[0] && win_coord[0] <right_top[0]) || (win_coord[0] > left_bottom[0] && win_coord[0] > right_top[0] ))
		return false;

	if( (win_coord[1] < left_bottom[1] && win_coord[1] <right_top[1]) || (win_coord[1] > left_bottom[1] && win_coord[1] > right_top[1] ))
		return false;

	return true;
}
