#pragma once 

#include "i_singular_points_finder.h"
#include <utility>
#include <vector>
#include <stdint.h>
#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
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
	void CalculateVerticesSurfaceType();
	void SegmentMolecularSurface(const size_t max_segm_size);
	void FindSegmentsGraphAndCenters();
	void CalculateSingularPointsTypes();
	void CalculateSingularPointsHistograms();

	MeshKeeper m_mesh_keeper;
	std::vector<Vertice> m_singular_points;
	std::vector<GraphNode<Vertice>> m_singular_points_graph;
	std::vector<std::array<uint8_t, kHistSize>> m_sing_pts_histo;
};

}