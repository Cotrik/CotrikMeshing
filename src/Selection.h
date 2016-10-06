#ifndef SELECTION_H
#define SELECTION_H

#include <vector>
using namespace std;

#include "math3d.h"
#include "include.h"

class Selection
{
public:
	Selection(void);
	~Selection(void);

	void draw_area();

	void highlight_selected_pts();

	void get_selected_pts_index(vector<int> &);
	void set_config(M3DVector3f *pts, int _nr, M3DVector2f _left_bottom, M3DVector2f _right_top, M3DMatrix44f model_view, M3DMatrix44f proj, int viewport[]);

	void set_config(const Mesh& mesh, int _nr, M3DVector2f _left_bottom, M3DVector2f _right_top, M3DMatrix44f model_view, M3DMatrix44f proj, int viewport[]);

private:
	bool drop_in_area(M3DVector3f x);
	void cal_selected_index();

private:
	M3DVector3f *pts;
	int nr;

	M3DVector2f left_bottom, right_top;
	M3DMatrix44f model_view, proj;
    int viewport[4];

	vector<int> vec_selected_pts_index;
};

#endif // endif
