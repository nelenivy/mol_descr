#pragma once
#include "opencv2\core\core.hpp"
#include "CommonUtilities/attributes_container.h"
#include "GraphLib\graph_structures.h"

namespace molecule_descriptor
{
struct Vertice
{
	Vertice() : coord(cv::Point3d(0, 0, 0)), normal(cv::Point3d(0, 0, 0)) { }
	
	Vertice(cv::Point3d coord_in,	cv::Point3d normal_in)
		: coord(coord_in), normal(normal_in) { }
	cv::Point3d Center() const { return coord; }
	cv::Point3d Normal() const { return normal; }

	cv::Point3d coord;
	cv::Point3d normal;
	mutable AttributesContainer attr;
};

struct Triangle
{
	Triangle() : a(nullptr), b(nullptr), c(nullptr) { }
	Triangle(Vertice* vert_1, Vertice* vert_2, Vertice* vert_3)
		: a(vert_1), b(vert_2), c(vert_3) { }

	Triangle(Vertice* vertices_in[3])
		: a(vertices_in[0]), b(vertices_in[1]), c(vertices_in[2]) { }
	cv::Point3d Center() const { 
		cv::Point3d tr_cent = a->coord +b->coord + c->coord;
		tr_cent.x /= 3.0;
		tr_cent.y /= 3.0;
		tr_cent.z /= 3.0; 
		return tr_cent;
	}

	cv::Point3d Normal() const { 
		cv::Point3d tr_normal = a->normal +b->normal + c->normal;
		tr_normal.x /= 3.0;
		tr_normal.y /= 3.0;
		tr_normal.z /= 3.0; 
		return tr_normal;
	}

	Vertice* a;
	Vertice* b;
	Vertice* c;
	mutable AttributesContainer attr;
};

double Distance(const Vertice& vertice_1, const Vertice& vertice_2);
double Distance(const Triangle& triangle_1, const Triangle& triangle_2);
int CountSameVertices(const Triangle& triangle_1, const Triangle& triangle_2);

typedef GraphNode<Vertice> MeshVertice;
typedef GraphNode<Triangle> MeshTriangle;

}