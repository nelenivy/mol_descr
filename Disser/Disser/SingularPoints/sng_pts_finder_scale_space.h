#pragma once 

#include "i_singular_points_finder.h"
#include <utility>
#include <vector>
#include <stdint.h>
#include <cctype>
#include "opencv2/core/core.hpp"
#include "../Common/singular_point.h"
#include "GraphLib/array_property_map.h"
#include "GraphLib\proxy_property_map.h"
#include "GraphLib\gradient_calculator.h"
#include "scale_space_blur.h"
#include "mesh_constructor.h"
#include "mesh_types.h"

namespace molecule_descriptor
{

enum DetectorFunctionType
{
	LOG, DOG, HESS_DET
};

DetectorFunctionType StringToDetectorFunctionType(const std::string& enum_name);

class SngPtsFinderScaleSpace : public ISingularPointsFinder
{
public:
	SngPtsFinderScaleSpace() { }
	virtual void CalcSingPtsFromCalculatedProperties(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
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
	virtual void GetVerticesWithDblProp(std::vector<std::pair<cv::Point3d, double>>& vertices_with_prop, const ISingularPointsFinder::SurfProperty prop_type);
	virtual void GetVerticesWithDblPropLevels(std::vector<std::vector<std::pair<cv::Point3d, double>>>& vertices_with_props_lev, 
		const ISingularPointsFinder::SurfProperty prop_type);
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
	virtual void SetMeanAndSigma(const std::vector<std::vector<double>>& mean_and_sigma)
	{
		m_mean_and_sigma = mean_and_sigma;
	}
	virtual void GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers)
	{
	};
	virtual size_t GetScaleSpaceLevelsNum() const 
	{
		return m_scale_space_levels_num;
	}
	virtual size_t GetDetectorFuncLevelsNum() const 
	{
		return m_blob_response_levels_num;
	}

	virtual void InitParams(int argc, char** argv);

	virtual void SetScaleSpace(const std::vector<std::vector<std::vector<double>>>& blurred_functions);
	virtual void SetDetectorFunction(const std::vector<std::vector<std::vector<double>>>& detector_functions);
	virtual void SetEigRatio(const std::vector<std::vector<double>>& eig_functions);
	virtual std::vector<double> GetSigmaValues();
private:
	Mesh& GetMesh(){ return const_cast<Mesh&>(m_mesh_keeper.GetMesh());/*m_filtered_mesh;*/}
	const Mesh& GetMesh() const { return const_cast<SngPtsFinderScaleSpace*>(this)->GetMesh();}
	const VerticesGraph& Vertices() const {  return GetMesh().vertices;}
	const TrianglesGraph& Triangles() const {  return GetMesh().triangles;}
	void Clear();//call before processing new data
	void CalculateVerticesSurfaceType();
	void CalcSingPtsFromCurvatureScales();
	void CalculatePropsInSingPts(const bool calc_prop_as_average);
	void CalcCurvature();
	void GetCoordinateMaps();
	void CalcTangentBasisConsistent();
	void CalcTangentBasisFast();
	//void CalcScaleSpacePropsHessianRatio();
	void CalcHessianOfProjectedLog();
	void CalculateDistanceMaps();
	void CalcScaleSpaceHessian();
	void CalcScalarBlobResponse();
	void CalcManifoldDeterminantBlobResponse();
	void ResccaleInputFunctions();
	void CalculateDOG();
	void CalculateLOG();

	void CalculateProjectedVectors();
	void FindSingPtsAsMaximumsOfScaleSpace();
	void FindSingPtsAsSeparateMaximumsOfLOG();
	void FindSingPtsAsCombinedMaximumsOfLOG();
	void FilterByEigenvaluesRatio();
	void FilterByNormedLogValue();

	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> VetrticesChargeMap;
	typedef ContPropMap<VerticesGraph, std::vector<size_t>, VERTEX> VetrticesTypeMap;
	typedef ContPropMap<VerticesGraph, std::vector<uint8_t>, VERTEX> VetrticesCurvMap;
	typedef ContPropMap<VerticesGraph, std::vector<Vertice>, VERTEX> FilteredCoordMap;

	typedef ContPropMap<SingularPointsGraph, std::vector<double>, VERTEX> SingPtsDoublePropMap;
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> DoubleVertGraphProp;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Point3d>, VERTEX> Point3dVertGraphProp;
	typedef ContPropMap<VerticesGraph, std::vector<uint8_t>, VERTEX> UInt8VertGraphProp;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> MatrixVertGraphProp;

	MeshKeeper m_mesh_keeper;
	DoubleVertGraphProp m_mean_curvature;
	DoubleVertGraphProp m_gaussian_curvature;

	std::vector<VertexDescriptor> m_maximums;
	std::vector<std::vector<VertexDescriptor>> m_maximums_with_levels;
	//FilteredCoordMap m_filtered_coord;
	//VetrticesTypeMap m_vertices_type_map;
	VetrticesChargeMap m_vertex_charge_map;
	VetrticesChargeMap m_vertex_lennard_jones_map;
	Point3dVertGraphProp m_vertex_charge_dir_map;
	Point3dVertGraphProp m_vertex_lennard_jones_dir_map;
	VetrticesCurvMap m_vertex_curv_type;
	std::vector<VetrticesCurvMap> m_vertex_curv_type_mesh_levels;

	std::vector<std::vector<double>> m_mean_and_sigma;

	int m_sing_pts_levels_num;
	int m_scale_space_levels_num;
	int m_blob_response_levels_num;
	int m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls;
	bool m_detect_blobs;
	DetectorFunctionType m_detector_type;
	bool m_use_euclid_distance;
	bool m_one_ring_neighb;
	bool m_combine_channels;

	double m_init_curv_sigma;
	double m_sigma_max;

	typedef ProxyPropMapVal<boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type, GetCoord<Vertice>> CoordMap;
	typedef ProxyPropMapVal<boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type, GetNormal<Vertice>> NormalMap;
	typedef ProxyPropMapVal<boost::property_map<const TrianglesGraph, boost::vertex_info_3d_t>::const_type, GetCoord<MeshTriangle>> CoordMapTriangle;
	ScaleSpaceBlurrer<VerticesGraph, CoordMap, GaussianKernel<cv::Point3d, double>> m_scale_space_blurrer;
	std::vector<std::vector<DoubleVertGraphProp>> m_components_blob_response;
	bool m_components_blob_response_calculated;
	std::vector<std::vector<DoubleVertGraphProp>> m_output_scale_space;
	std::vector<std::vector<Point3dVertGraphProp>> m_output_spher_scale_space;
	bool m_use_spherical;
	bool m_scale_space_calculated;
	std::vector<DoubleVertGraphProp> m_input_prop_map;
	std::vector<Point3dVertGraphProp> m_input_spherical_prop_map; 
	bool m_hessian_ratio_calculated;

	struct PCAProjecter
	{
		cv::PCA pca;
		cv::Mat_<double> vect_to_project_on;
	};
	typedef ContPropMap<VerticesGraph, std::vector<PCAProjecter>, VERTEX> PCAProjecterMap;
	std::vector<PCAProjecterMap> m_scale_space_projecter;
	bool m_use_central_projector;
	bool m_filter_by_eigenvalues;
	bool m_scale_extr;

	std::vector<DoubleVertGraphProp> m_blob_response;
	std::vector<DoubleVertGraphProp> m_scale_space_projected;
	std::vector<std::vector<DoubleVertGraphProp>> m_projecters_coords;

	typedef ContPropMap<VerticesGraph, std::vector<cv::Mat_<double>>, VERTEX> TangentBasisMap;
	TangentBasisMap m_tangent_basis_map;

	std::vector<HessianMatrixCalculator> m_hessian_map_calculators;
	std::vector<HessianMatrixCalculatorSpherical> m_spherical_hessian_map_calculators;
	std::vector<std::vector<DoubleVertGraphProp>> m_props_hessian_ratio;
	std::vector<std::vector<DoubleVertGraphProp>> m_grad_dx;
	std::vector<std::vector<DoubleVertGraphProp>> m_grad_dy;
	std::vector<DoubleVertGraphProp> m_projected_grad_x;
	std::vector<DoubleVertGraphProp> m_projected_grad_y;
	std::vector<DoubleVertGraphProp> m_projected_grad_norm;
	std::vector<DoubleVertGraphProp> m_props_hessian_ratio_of_proj;
	std::vector<DoubleVertGraphProp> m_scalar_blob_response;
	std::vector<DoubleVertGraphProp> m_manifold_det_blob_response;
	std::vector<std::vector<DoubleVertGraphProp>> m_hessian_det;
	std::vector<std::vector<MatrixVertGraphProp>> m_hessian;
	std::vector<std::vector<MatrixVertGraphProp>> m_hessian_spherical;
	double m_ratio_thresh;

	std::vector<DoubleVertGraphProp> m_props_hessian_ratio_of_proj_LOG;
	//distance maps
	cv::Mat_<double> m_vert_vert_dist;
	cv::Mat_<double> m_vert_tr_dist;

	CoordMap m_coord_map;
	NormalMap m_norm_map;
	CoordMapTriangle coord_map_tr;
	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map;
	boost::property_map<const TrianglesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map_tr;

	enum ChannelsCombination
	{
		PCA, NORM, DETECTOR_NORM, SCALAR_BLOB_RESPONSE, MANIFOLD_DET_BLOB_RESPONSE
	};

	ChannelsCombination m_channel_combination;

	/*UInt8VertGraphProp*/
	DoubleVertGraphProp m_uncoincided_vertices;
};

}