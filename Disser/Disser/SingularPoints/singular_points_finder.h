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
		const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges);
	virtual void GetMarkedSingularPoints(std::vector<MarkedSingularPoint>& marked_singular_points);
	virtual void GetNonMarkedSingularPoints(std::vector<NonMarkedSingularPoint>& non_marked_singular_points);
	virtual void GetSingularPointsHisto(std::vector<HistogramSingularPoint<kHistSize>>& singular_points_hist);
	//TODO: fix compile-time properties num. Maybe put pairs into this class
	virtual size_t GetTypesNum();
	virtual void GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers);
	virtual void GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types);
	virtual void Release() { delete this; }
private:
	void Clear();//call before processing new data
	void CalculateAllPotentials(const std::vector<std::pair<cv::Point3d, double>>& charges);
	void CalculateVerticesSurfaceType();
	void SegmentMolecularSurface(const size_t max_segm_size);
	void FindSegmentsGraphAndCenters();
	void CalculateSingularPointsTypes();
	void CalculateSingularPointsHistograms();
	
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> VetrticesChargeMap;
	typedef ContPropMap<VerticesGraph, std::vector<size_t>, VERTEX> VetrticesTypeMap;
	typedef ContPropMap<VerticesGraph, std::vector<uint8_t>, VERTEX> VetrticesCurvMap;
	typedef ContPropMap<VerticesGraph, std::vector<size_t>, VERTEX> VetrticesSegmMap;
	typedef ContPropMap<VerticesGraph, std::vector<std::map<size_t, uint8_t>>, VERTEX> VetrticesSegmScoreMap;

	typedef ContPropMap<SingularPointsGraph, 
	std::vector<std::array<uint8_t, kHistSize>>, VERTEX> SingPtsHistMap;

	typedef ContPropMap<TrianglesGraph, std::vector<uint8_t>, VERTEX> TrianglesCurvMap;
	typedef ContPropMap<TrianglesGraph, std::vector<size_t>, VERTEX> TrianglesSegmMap;

	MeshKeeper m_mesh_keeper;
	VetrticesTypeMap m_vertices_type_map;
	VetrticesChargeMap m_vertex_charge_map;
	VetrticesCurvMap m_vertex_curv_type;
	VetrticesSegmMap m_vertex_segm;
	VetrticesSegmScoreMap m_vertex_score_map;
	TrianglesCurvMap m_triangle_curv_type;
	TrianglesSegmMap m_triangles_segm;
	std::array<TrianglesSegmMap, 3>  m_type_triangle_segm;
	std::vector<VertexDescriptor> m_singular_points;
	SingularPointsGraph m_singular_points_graph;
	SingPtsHistMap m_sing_pts_histo;
};

}