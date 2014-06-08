#pragma once 

#include "i_singular_points_finder.h"
#include <utility>
#include <vector>
#include <stdint.h>
#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "GraphLib/array_property_map.h"
#include "mesh_constructor.h"
#include "mesh_types.h"

namespace molecule_descriptor
{

class SingularPointsFinder : public ISingularPointsFinder
{
public:
	SingularPointsFinder() { }
	virtual void Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
		const std::vector<cv::Point3i>& triangles, 
		const std::vector<std::pair<cv::Point3d, double>>& charges, const std::vector<std::pair<cv::Point3d, double>>& wdv_radii,
		const bool calc_prop_as_average);
	virtual void GetMarkedSingularPoints(std::vector<MarkedSingularPoint>& marked_singular_points);
	virtual void GetNonMarkedSingularPoints(std::pair<std::vector<NonMarkedSingularPoint>, std::vector<size_t>>& non_marked_singular_points);
	virtual void GetSingularPointsHisto(std::vector<HistogramSingularPoint<kHistSize>>& singular_points_hist);
	//TODO: fix compile-time properties num. Maybe put pairs into this class
	virtual void GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers);
	virtual void GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types);
	virtual void Release() { delete this; }
private:
	Mesh& GetMesh(){ return const_cast<Mesh&>(m_mesh_keeper.GetMesh());/*m_filtered_mesh;*/}
	const Mesh& GetMesh() const { return const_cast<const SingularPointsFinder*>(this)->GetMesh();}
	void Clear();//call before processing new data
	void CalculateAllPotentials(const std::vector<std::pair<cv::Point3d, double>>& charges, const Mesh& mesh);
	void CalculateLennardJonesPotentials(const std::vector<std::pair<cv::Point3d, double>>& wdv_radii, const Mesh& mesh);
	void CalculateVerticesSurfaceType(const Mesh& mesh);
	void SegmentMolecularSurface(const size_t max_segm_size, const Mesh& mesh);
	void FindSegmentsGraphAndCenters(const Mesh& mesh);
	void CalculateSingularPointsHistograms();
	void CalculatePropsInSingPts(const bool calc_prop_as_average);

	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> VetrticesChargeMap;
	typedef ContPropMap<VerticesGraph, std::vector<size_t>, VERTEX> VetrticesTypeMap;
	typedef ContPropMap<VerticesGraph, std::vector<uint8_t>, VERTEX> VetrticesCurvMap;
	typedef ContPropMap<VerticesGraph, std::vector<size_t>, VERTEX> VetrticesSegmMap;
	typedef ContPropMap<VerticesGraph, std::vector<std::map<size_t, uint8_t>>, VERTEX> VetrticesSegmScoreMap;
	typedef ContPropMap<VerticesGraph, std::vector<Vertice>, VERTEX> FilteredCoordMap;
	typedef ContPropMap<SingularPointsGraph, 
	std::vector<std::array<uint8_t, kHistSize>>, VERTEX> SingPtsHistMap;

	typedef ContPropMap<SingularPointsGraph, std::vector<double>, VERTEX> SingPtsDoublePropMap;

	typedef ContPropMap<TrianglesGraph, std::vector<uint8_t>, VERTEX> TrianglesCurvMap;
	typedef ContPropMap<TrianglesGraph, std::vector<size_t>, VERTEX> TrianglesSegmMap;

	template<typename GraphPropMap>
	void CalcPropInSingPts(const GraphPropMap& graph_prop_map, const bool calc_prop_as_average, SingPtsDoublePropMap& sing_pts_prop_map);

	MeshKeeper m_mesh_keeper;
	Mesh m_filtered_mesh;
	ContPropMap<VerticesGraph, std::vector<double>, VERTEX> m_curvature_1;
	ContPropMap<VerticesGraph, std::vector<double>, VERTEX> m_curvature_0;

	//FilteredCoordMap m_filtered_coord;
	//VetrticesTypeMap m_vertices_type_map;
	VetrticesChargeMap m_vertex_charge_map;
	VetrticesChargeMap m_vertex_lennard_jones_map;
	VetrticesCurvMap m_vertex_curv_type;
	VetrticesSegmMap m_vertex_segm;
	VetrticesSegmScoreMap m_vertex_score_map;
	TrianglesCurvMap m_triangle_curv_type;
	TrianglesSegmMap m_triangles_segm;
	std::array<TrianglesSegmMap, 3>  m_type_triangle_segm;
	std::vector<VertexDescriptor> m_singular_points;
	SingularPointsGraph m_singular_points_graph;
	SingPtsHistMap m_sing_pts_histo;
	SingPtsDoublePropMap m_sing_pts_potential;
	SingPtsDoublePropMap m_sing_pts_lennard_jones;
};

}