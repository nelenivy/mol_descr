#pragma once
#include "opencv2\core\core.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/subgraph.hpp"
#include "boost/graph/properties.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_iterator.hpp"

#define BOOST_DEF_PROPERTY(KIND, NAME) \
	enum KIND##_##NAME##_t { KIND##_##NAME }; \
	BOOST_INSTALL_PROPERTY(KIND, NAME)

namespace boost
{
	BOOST_DEF_PROPERTY(vertex, info_3d);
	BOOST_DEF_PROPERTY(vertex, segment_num);
	BOOST_DEF_PROPERTY(vertex, parent);
}
#undef BOOST_DEF_PROPERTY

namespace molecule_descriptor
{

//for proxy container map
template <typename VertexType>
struct GetCoord;
template <typename VertexType>
struct GetNormal;

struct Vertice
{
	Vertice() : coord(cv::Point3d(0, 0, 0)), normal(cv::Point3d(0, 0, 0)) { }
	
	Vertice(cv::Point3d coord_in,	cv::Point3d normal_in)
		: coord(coord_in), normal(normal_in) { }
	cv::Point3d Center() const { return coord; }
	cv::Point3d Normal() const { return normal; }

	cv::Point3d coord;
	cv::Point3d normal;
};

template <>
struct GetCoord<Vertice>
{
	cv::Point3d& operator()(Vertice& vert) const 
	{
		return vert.coord;
	} 
	const cv::Point3d& operator()(const Vertice& vert) const {return vert.coord;} 
};
template <>
struct GetNormal<Vertice>
{
	cv::Point3d& operator()(Vertice& vert) const 
	{
		return vert.normal;
	} 
	const cv::Point3d& operator()(const Vertice& vert) const {return vert.normal;} 
};

inline bool operator==(const Vertice& vert_1, const Vertice& vert_2)
{
	return (vert_1.coord == vert_2.coord) && (vert_1.normal == vert_2.normal);
}
inline bool operator!=(const Vertice& vert_1, const Vertice& vert_2)
{
	return !(vert_1 == vert_2);
}

template <typename VertexDescriptor, class Graph>
class Triangle
{
public:
	Triangle() : a(VertexDescriptor()), b(VertexDescriptor()), 
		c(VertexDescriptor()), m_graph(nullptr) { }
	Triangle(const VertexDescriptor vert_1, 
		const VertexDescriptor vert_2, const VertexDescriptor vert_3, Graph& graph)
		: a(vert_1), b(vert_2), c(vert_3), m_graph(&graph) 
	{
		SetVerticesAndGraph(vert_1, vert_2, vert_3, graph);
	}

	void SetVerticesAndGraph(const VertexDescriptor vert_1, 
		const VertexDescriptor vert_2, const VertexDescriptor vert_3, Graph& graph)
	{
		a = vert_1;
		b = vert_2;
		c = vert_3;
		m_graph = &graph;
		RecalculateProperties();
	}
	void RecalculateProperties()
	{
		CalcCenter();
		CalcNormal();
		CalcArea();
	}
	double Area() const {
		return area;
	}
	const cv::Point3d& Normal() const { 
		return normal;
	}
	const cv::Point3d& Center() const {
		return center;
	}
	const VertexDescriptor& GetA() const
	{
		return a;
	}
	const VertexDescriptor& GetB() const
	{
		return b;
	}
	const VertexDescriptor& GetC() const
	{
		return c;
	}
private:
	void CalcCenter() { 
		const auto pmap = get(boost::vertex_info_3d, *m_graph);
		cv::Point3d tr_cent = pmap[a].Center() + pmap[b].Center() + pmap[c].Center();
		tr_cent *= 1.0 / 3.0;
		center = tr_cent;
	}

	void CalcNormal() { 
		const auto pmap =  get(boost::vertex_info_3d, *m_graph);
		cv::Point3d tr_normal = pmap[a].Normal() + pmap[b].Normal() + pmap[c].Normal();
		tr_normal.x *= 1.0 / 3.0;
		normal = tr_normal;
	}
	void CalcArea() 
	{
		const auto pmap = get(boost::vertex_info_3d, *m_graph);
		const double side_a = cv::norm(pmap[a].Center() - pmap[b].Center());
		const double side_b = cv::norm(pmap[a].Center() - pmap[c].Center());
		const double side_c = cv::norm(pmap[b].Center() - pmap[c].Center());
		const double half_perim = (side_a + side_b + side_c) / 2.0;
		area = sqrt(half_perim * (half_perim - side_a) * (half_perim - side_b) * (half_perim - side_c));
	}
	VertexDescriptor a;
	VertexDescriptor b;
	VertexDescriptor c;
	cv::Point3d normal;
	cv::Point3d center;
	double area;

	Graph* m_graph;
};

inline double Distance(const Vertice& vertice_1, const Vertice& vertice_2)
{
	return norm(vertice_1.coord - vertice_2.coord);
}

template <typename VertexDescriptor, class Graph>
double Distance(const Triangle<VertexDescriptor, Graph>& triangle_1, 
				const Triangle<VertexDescriptor, Graph>& triangle_2)
{
	return norm(triangle_1.Center() - triangle_2.Center());
}

template <typename VertexDescriptor, class Graph>
int CountSameVertices(const Triangle<VertexDescriptor, Graph>& triangle_1, 
					  const Triangle<VertexDescriptor, Graph>& triangle_2)
{
	int same_vertices = 0;

	VertexDescriptor vertices_1[] = {triangle_1.GetA(), triangle_1.GetB(), triangle_1.GetC()};
	VertexDescriptor vertices_2[] = {triangle_2.GetA(), triangle_2.GetB(), triangle_2.GetC()};

	for (int vert_1_ind = 0; vert_1_ind < 3; vert_1_ind++)
	{
		for (int vert_2_ind = 0; vert_2_ind < 3; vert_2_ind++)
		{
			same_vertices += (vertices_1[vert_1_ind] == vertices_2[vert_2_ind]);
		}
	}

	return same_vertices;
}

typedef boost::subgraph<
	boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property<boost::vertex_info_3d_t, Vertice>, 
	boost::property<boost::edge_index_t, size_t>, 
	boost::no_property, boost::vecS>> VerticesGraph;
typedef boost::graph_traits<VerticesGraph>::vertex_descriptor VertexDescriptor;
typedef boost::graph_traits<VerticesGraph>::vertex_iterator VertexIter;
typedef boost::graph_traits<VerticesGraph>::adjacency_iterator AdjacencyIter;

typedef boost::subgraph<
	boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property<boost::vertex_info_3d_t, Vertice,
	boost::property<boost::vertex_parent_t, VertexDescriptor>>, 
	boost::property<boost::edge_index_t, size_t>, 
	boost::no_property, boost::vecS>> SingularPointsGraph;

typedef Triangle<VertexDescriptor, VerticesGraph> MeshTriangle;
typedef boost::subgraph<
	boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property<boost::vertex_info_3d_t, MeshTriangle>, 
	boost::property<boost::edge_index_t, size_t>, 
	boost::no_property, boost::vecS>> TrianglesGraph;
typedef boost::graph_traits<TrianglesGraph>::vertex_descriptor 
	TriangleDescriptor;
typedef boost::graph_traits<TrianglesGraph>::vertex_iterator TriangleIter;

struct Mesh
{
	VerticesGraph vertices;
	TrianglesGraph triangles;
};

template <>
struct GetCoord<MeshTriangle>
{
	const cv::Point3d& operator()(const MeshTriangle& vert) const 
	{
		return vert.Center();
	} 
};

}