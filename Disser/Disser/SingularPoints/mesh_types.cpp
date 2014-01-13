#include "mesh_types.h"
#include <algorithm>
#include <math.h>

namespace molecule_descriptor
{

double Distance(const Vertice& vertice_1, const Vertice& vertice_2)
{
	return norm(vertice_1.coord - vertice_2.coord);
}

double Distance(const Triangle& triangle_1, const Triangle& triangle_2)
{
	cv::Point3d tr_cent_1 = triangle_1.a->coord + triangle_1.b->coord + triangle_1.c->coord;
	tr_cent_1.x /= 3.0;
	tr_cent_1.y /= 3.0;
	tr_cent_1.z /= 3.0;
	cv::Point3d tr_cent_2 = triangle_2.a->coord + triangle_2.b->coord + triangle_2.c->coord;
	tr_cent_2.x /= 3.0;
	tr_cent_2.y /= 3.0;
	tr_cent_2.z /= 3.0;
	return norm(tr_cent_1 - tr_cent_2);
}

int CountSameVertices(const Triangle& triangle_1, const Triangle& triangle_2)
{
	int same_vertices = 0;

	Vertice* vertices_1[] = {triangle_1.a, triangle_1.b, triangle_1.c};
	Vertice* vertices_2[] = {triangle_2.a, triangle_2.b, triangle_2.c};

	for (int vert_1_ind = 0; vert_1_ind < 3; vert_1_ind++)
	{
		for (int vert_2_ind = 0; vert_2_ind < 3; vert_2_ind++)
		{
			same_vertices += (vertices_1[vert_1_ind] == vertices_2[vert_2_ind]);
		}
	}

	return same_vertices;
}

}