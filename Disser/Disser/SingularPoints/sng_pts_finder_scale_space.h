#pragma once 

#include "i_singular_points_finder.h"
#include <utility>
#include <vector>
#include <stdint.h>
#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "GraphLib/array_property_map.h"
#include "GraphLib\proxy_property_map.h"
#include "scale_space_blur.h"
#include "mesh_constructor.h"
#include "mesh_types.h"

namespace molecule_descriptor
{

class SngPtsFinderScaleSpace : public ISingularPointsFinder
{
public:
	SngPtsFinderScaleSpace() { }
	virtual void Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
		const std::vector<cv::Point3i>& triangles, 
		const std::vector<std::pair<cv::Point3d, double>>& charges, const std::vector<std::pair<cv::Point3d, double>>& wdv_radii,
		const bool calc_prop_as_average);
	virtual void CalcOnlyProps(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
		const std::vector<cv::Point3i>& triangles, 
		const std::vector<std::pair<cv::Point3d, double>>& charges, 
		const std::vector<std::pair<cv::Point3d, double>>& wdv_radii);
	virtual void GetNonMarkedSingularPoints(std::pair<std::vector<NonMarkedSingularPoint>, std::vector<size_t>>& non_marked_singular_points);
	virtual void GetNonMarkedSingularPointsLevels(std::pair<std::vector<std::vector<NonMarkedSingularPoint>>, 
		std::vector<std::vector<size_t>>>& non_marked_singular_points);
	virtual void GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types);
	virtual void GetVerticesWithTypesLevels(std::vector<std::vector<std::pair<cv::Point3d, size_t>>>& vertices_with_types);
	virtual void AppendProp(std::vector<double>& prop, SurfProperty prop_type)
	{
		if (prop_type == ISingularPointsFinder::GAUSS_CURV)
		{
			prop.insert(prop.end(), m_gaussian_curvature.begin(),m_gaussian_curvature.end());
		}
		else if (prop_type == ISingularPointsFinder::MEAN_CURV)
		{	
			prop.insert(prop.end(), m_mean_curvature.begin(),m_mean_curvature.end());
		}
		else if (prop_type == ISingularPointsFinder::ELECTR_POTENT)
		{	
			prop.insert(prop.end(), m_vertex_charge_map.begin(),m_vertex_charge_map.end());
		}
		else if (prop_type == ISingularPointsFinder::STERIC_POTENT)
		{	
			prop.insert(prop.end(), m_vertex_lennard_jones_map.begin(),m_vertex_lennard_jones_map.end());
		}	
	}
	virtual void SetMeanAndSigma(std::vector<std::vector<double>> mean_and_sigma)
	{
		m_mean_and_sigma = mean_and_sigma;
	}
	virtual void GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers)
	{
	};
	void InitParams(int argc, char** argv);

private:
	Mesh& GetMesh(){ return const_cast<Mesh&>(m_mesh_keeper.GetMesh());/*m_filtered_mesh;*/}
	const Mesh& GetMesh() const { return const_cast<const SngPtsFinderScaleSpace*>(this)->GetMesh();}
	void Clear();//call before processing new data
	void CalculateVerticesSurfaceType(const Mesh& mesh);
	void CalcSingPtsFromCurvatureScales(const Mesh& mesh);
	void CalculatePropsInSingPts(const bool calc_prop_as_average);
	void CalcShiftsMaximums();
	void CalcCurvature(const Mesh& mesh);
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> VetrticesChargeMap;
	typedef ContPropMap<VerticesGraph, std::vector<size_t>, VERTEX> VetrticesTypeMap;
	typedef ContPropMap<VerticesGraph, std::vector<uint8_t>, VERTEX> VetrticesCurvMap;
	typedef ContPropMap<VerticesGraph, std::vector<Vertice>, VERTEX> FilteredCoordMap;

	typedef ContPropMap<SingularPointsGraph, std::vector<double>, VERTEX> SingPtsDoublePropMap;
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> DoubleVertGraphProp;

	MeshKeeper m_mesh_keeper;
	Mesh m_filtered_mesh;
	std::vector<Mesh> m_filtered_mesh_levels;
	std::vector<Mesh> m_filtered_mesh_levels_2;
	DoubleVertGraphProp m_mean_curvature;
	DoubleVertGraphProp m_gaussian_curvature;

	std::vector<DoubleVertGraphProp> m_mean_curvature_scales;
	std::vector<DoubleVertGraphProp> m_gaussian_curvature_scales;

	std::vector<VertexDescriptor> m_maximums;
	std::vector<std::vector<VertexDescriptor>> m_maximums_with_levels;
	//FilteredCoordMap m_filtered_coord;
	//VetrticesTypeMap m_vertices_type_map;
	VetrticesChargeMap m_vertex_charge_map;
	VetrticesChargeMap m_vertex_lennard_jones_map;
	VetrticesCurvMap m_vertex_curv_type;
	std::vector<VetrticesCurvMap> m_vertex_curv_type_mesh_levels;

	std::vector<std::vector<double>> m_mean_and_sigma;

	int m_sing_pts_levels_num;
	int m_scale_space_levels_num;
	int m_diff_between_sing_pts_levels_and_scale_space_levels;
	bool m_detect_blobs;
	bool m_combine_channels;

	double m_curv_sigma;
	double m_init_curv_sigma;
	double m_sigma_max;

	typedef ProxyPropMap<
		boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type, GetCoord<Vertice>> CoordMap;
	ScaleSpaceBlurrer<VerticesGraph, CoordMap, GaussianKernel<cv::Point3d, double>> m_scale_space_blurrer;
	std::vector<std::vector<DoubleVertGraphProp>> m_output_scale_space;


};

}