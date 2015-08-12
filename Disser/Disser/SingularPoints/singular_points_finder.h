#pragma once 

#include "i_singular_points_finder.h"
#include <utility>
#include <vector>
#include <stdint.h>
#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "GraphLib/array_property_map.h"
#include "GraphLib\proxy_property_map.h"
#include "mesh_constructor.h"
#include "mesh_types.h"

namespace molecule_descriptor
{

class SngPtsFinderSegmentation : public ISingularPointsFinder
{
public:
	SngPtsFinderSegmentation() { }
	virtual void Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
		const std::vector<cv::Point3i>& triangles, 
		const std::vector<std::pair<cv::Point3d, double>>& charges, const std::vector<std::pair<cv::Point3d, double>>& wdv_radii,
		const bool calc_prop_as_average);
	virtual void GetNonMarkedSingularPoints(std::pair<std::vector<NonMarkedSingularPoint>, std::vector<size_t>>& non_marked_singular_points);
	virtual void GetNonMarkedSingularPointsLevels(std::pair<std::vector<std::vector<NonMarkedSingularPoint>>, std::vector<std::vector<size_t>>>& non_marked_singular_points);
	virtual void GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers);
	virtual void GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types);
	virtual void InitParams(int argc, char** argv) {}
private:
	Mesh& GetMesh(){ return const_cast<Mesh&>(m_mesh_keeper.GetMesh());/*m_filtered_mesh;*/}
	const Mesh& GetMesh() const { return const_cast<const SngPtsFinderSegmentation*>(this)->GetMesh();}
	void Clear();//call before processing new data
	void CalculateVerticesSurfaceType(const Mesh& mesh);
	void CalculateCurvature(const Mesh& mesh);
	void SegmentMolecularSurface(const size_t max_segm_size, const Mesh& mesh);
	void FindSegmentsGraphAndCenters(const Mesh& mesh);
	void CalculatePropsInSingPts(const bool calc_prop_as_average);
	void CalcSegmentsArea(/*const Mesh& mesh*/);
	void CalcShiftsMaximums();
	void CalculateDistanceMaps(const Mesh& mesh);
	void CalcTangentBasis(const Mesh& mesh);

	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> VetrticesChargeMap;
	typedef ContPropMap<VerticesGraph, std::vector<size_t>, VERTEX> VetrticesTypeMap;
	typedef ContPropMap<VerticesGraph, std::vector<uint8_t>, VERTEX> VetrticesCurvMap;
	typedef ContPropMap<VerticesGraph, std::vector<size_t>, VERTEX> VetrticesSegmMap;
	typedef ContPropMap<VerticesGraph, std::vector<std::map<size_t, uint8_t>>, VERTEX> VetrticesSegmScoreMap;
	typedef ContPropMap<VerticesGraph, std::vector<Vertice>, VERTEX> FilteredCoordMap;

	typedef ContPropMap<SingularPointsGraph, std::vector<double>, VERTEX> SingPtsDoublePropMap;

	typedef ContPropMap<TrianglesGraph, std::vector<uint8_t>, VERTEX> TrianglesCurvMap;
	typedef ContPropMap<TrianglesGraph, std::vector<size_t>, VERTEX> TrianglesSegmMap;
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> DoubleVertGraphProp;

	template<typename GraphPropMap>
	void CalcPropInSingPts(const GraphPropMap& graph_prop_map, const bool calc_prop_as_average, SingPtsDoublePropMap& sing_pts_prop_map);

	MeshKeeper m_mesh_keeper;
	DoubleVertGraphProp m_mean_curvature;
	DoubleVertGraphProp m_gaussian_curvature;

	std::vector<VertexDescriptor> m_maximums;
	VetrticesChargeMap m_vertex_charge_map;
	VetrticesChargeMap m_vertex_lennard_jones_map;
	VetrticesCurvMap m_vertex_curv_type;
	VetrticesSegmMap m_vertex_segm;
	VetrticesSegmScoreMap m_vertex_score_map;
	TrianglesCurvMap m_triangle_curv_type;
	TrianglesSegmMap m_triangles_segm;
	std::array<TrianglesSegmMap, 3>  m_type_triangle_segm;
	SingularPointsGraph m_singular_points_graph;
	SingPtsDoublePropMap m_sing_pts_potential;
	SingPtsDoublePropMap m_sing_pts_lennard_jones;
	SingPtsDoublePropMap m_sing_pts_segm_area;

	//distance maps
	typedef ProxyPropMap<boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type, GetCoord<Vertice>> CoordMap;
	typedef ProxyPropMap<boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type, GetNormal<Vertice>> NormalMap;

	cv::Mat_<double> m_vert_vert_dist;
	cv::Mat_<double> m_vert_tr_dist;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> TangentBasisMap;
	TangentBasisMap m_tangent_basis_map;

};

}