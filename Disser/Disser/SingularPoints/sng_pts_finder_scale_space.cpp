#include "sng_pts_finder_scale_space.h"

#include <algorithm>
#include <utility>
#include "opencv2/core/core.hpp"
#include "opencv2/core/core_c.h"
#include "SingularPoints/mesh_types.h"
#include "SingularPoints\mesh_operations.h"
#include "SingularPoints\pca_projections.h"
#include "GraphLib\connected_components_segmentator.h"
#include "GraphLib\graph_operations.h"
#include "GraphLib\graph_filter.h"
#include "GraphLib\curvature_calculator.h"
#include "GraphLib\graph_functions.h"
#include "GraphLib\coordinates_transform.h"
#include "GraphLib\graph_dist_calculate.h"
#include "CommonUtilities\common_functions.h"
#include "InputOutput/params_reader.h"

#include "points_keeper.h"
#include "properties_calculator.h"

namespace molecule_descriptor
{

using namespace std;
using namespace cv;

DetectorFunctionType StringToDetectorFunctionType(const std::string& enum_name)
{
	std::string lowercase_name = enum_name;
	std::transform(lowercase_name.begin(), lowercase_name.end(), lowercase_name.begin(), std::tolower);
	if (enum_name == std::string("log"))
	{
		return DetectorFunctionType::LOG;
	}
	else if (enum_name == std::string("dog"))
	{
		return DetectorFunctionType::DOG;
	}
	else if (enum_name == std::string("hess_det"))
	{
		return DetectorFunctionType::HESS_DET;
	}
	else
	{
		CV_Assert(0);
	}
}

void SngPtsFinderScaleSpace::InitParams(int argc, char** argv)
{
	ReadParamFromCommandLineWithDefault(argc, argv, "-mesh_levels_num", m_sing_pts_levels_num, 10);
	ReadParamFromCommandLineWithDefault(argc, argv, "-detect_blobs", m_detect_blobs, true);
	ReadParamFromCommandLineWithDefault(argc, argv, "-use_euclid_distance", m_use_euclid_distance, false);
	ReadParamFromCommandLineWithDefault(argc, argv, "-combine_channels", m_combine_channels, true);
	ReadParamFromCommandLineWithDefault(argc, argv, "-one_ring_neighb", m_one_ring_neighb, true);

	ReadParamFromCommandLineWithDefault(argc, argv, "-init_curv_sigma", m_init_curv_sigma, 0.33);
	ReadParamFromCommandLineWithDefault(argc, argv, "-sigma_max", m_sigma_max, 2.89);
	ReadParamFromCommandLineWithDefault(argc, argv, "-use_central_projector", m_use_central_projector, false);
	ReadParamFromCommandLineWithDefault(argc, argv, "-filter_by_eigenvalues", m_filter_by_eigenvalues, false);
	ReadParamFromCommandLineWithDefault(argc, argv, "-use_spherical", m_use_spherical, false);
	std::string detector_type_name;
	ReadParamFromCommandLineWithDefault(argc, argv, "-detector_type", detector_type_name, std::string("hess_det"));
	m_detector_type = StringToDetectorFunctionType(detector_type_name);
	ReadParamFromCommandLineWithDefault(argc, argv, "-scale_extr", m_scale_extr, true);
	ReadParamFromCommandLineWithDefault(argc, argv, "-ratio_thresh", m_ratio_thresh, 8.0);
	
	std::string cnannels_combining_name;
	ReadParamFromCommandLineWithDefault(argc, argv, "-cnannels_combining", cnannels_combining_name, std::string("PCA"));
	m_channel_combination = SCALAR_BLOB_RESPONSE;
	
	if (cnannels_combining_name == "PCA")
	{
		m_channel_combination = PCA;
	}
	else if (cnannels_combining_name == "NORM")
	{
		m_channel_combination = NORM;
	}
	else if (cnannels_combining_name == "DETECTOR_NORM")
	{
		m_channel_combination = DETECTOR_NORM;
	}
	else if (cnannels_combining_name == "GENERALIZED_HESSIAN")
	{
		m_channel_combination = SCALAR_BLOB_RESPONSE;
	}
	else if (cnannels_combining_name == "MANIFOLD_DET_BLOB_RESPONSE")
	{
		m_channel_combination = MANIFOLD_DET_BLOB_RESPONSE;
	}
	

	if (m_detector_type == HESS_DET)
	{
		CV_Assert(m_channel_combination == SCALAR_BLOB_RESPONSE || m_channel_combination == MANIFOLD_DET_BLOB_RESPONSE);
	}

	m_scale_space_levels_num = m_sing_pts_levels_num + (m_detect_blobs ? 3 : 2);
	m_blob_response_levels_num = m_scale_space_levels_num - 1;
	m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls = 1;

	//INIT SCALE-SPACE
	const bool kAdditive = false;
	const double d_mult_sigma = pow(m_sigma_max / m_init_curv_sigma, 1.0 / (m_sing_pts_levels_num));
	std::cout << "d_mult_sigma " << d_mult_sigma << std::endl;
	m_scale_space_blurrer.Init(m_init_curv_sigma, d_mult_sigma, 0.0, kAdditive, false, 1.0);
}

template <typename Graph, typename PropMap>
void CopyVectorToMap(const std::vector<double>& input, const Graph& graph, PropMap& prop_map)
{
	//std::cout << input.size() << " " << num_vertices(graph) << "\n";
	CV_Assert(input.size() == num_vertices(graph));
	prop_map.SetGraph(graph);
	for (size_t ind = 0; ind < input.size(); ++ind)
	{
		prop_map[vertex(ind, graph)] = input[ind];
	}
}

void SngPtsFinderScaleSpace::SetScaleSpace(const std::vector<std::vector<std::vector<double>>>& blurred_functions)
{
	CV_Assert(blurred_functions.size() == m_scale_space_levels_num);
	m_output_scale_space.resize(m_scale_space_levels_num);

	for (size_t curr_lev = 0; curr_lev < m_scale_space_levels_num; ++curr_lev)
	{
		//std::cout << curr_lev << "\n";
		CV_Assert(blurred_functions[curr_lev].size() == ISingularPointsFinder::SURF_PROPS_NUM);
		m_output_scale_space[curr_lev].resize(ISingularPointsFinder::SURF_PROPS_NUM);

		for (int curr_prop = 0; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
		{
			//std::cout << curr_prop << "\n";

			const auto& input_cont = blurred_functions[curr_lev][curr_prop];
			auto& output_cont = m_output_scale_space[curr_lev][curr_prop];
			CopyVectorToMap(input_cont, Vertices(), output_cont);
		}
	}

	m_scale_space_calculated = true;
}

void SngPtsFinderScaleSpace::SetDetectorFunction(const std::vector<std::vector<std::vector<double>>>& detector_functions)
{
	CV_Assert(detector_functions.size() == m_blob_response_levels_num);
	m_components_blob_response.resize(m_blob_response_levels_num);

	for (size_t curr_lev = 0; curr_lev < m_blob_response_levels_num; ++curr_lev)
	{
		CV_Assert(detector_functions[curr_lev].size() == ISingularPointsFinder::SURF_PROPS_NUM);
		m_components_blob_response[curr_lev].resize(ISingularPointsFinder::SURF_PROPS_NUM);

		for (int curr_prop = 0; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
		{
			const auto& input_cont = detector_functions[curr_lev][curr_prop];
			auto& output_cont = m_components_blob_response[curr_lev][curr_prop];
			CopyVectorToMap(input_cont, Vertices(), output_cont);
		}
	}

	m_components_blob_response_calculated = true;
}

void SngPtsFinderScaleSpace::SetEigRatio(const std::vector<std::vector<double>>& eig_ratio_levels)
{
	CV_Assert(eig_ratio_levels.size() == m_blob_response_levels_num);
	m_props_hessian_ratio_of_proj.resize(m_blob_response_levels_num);

	for (size_t curr_lev = 0; curr_lev < m_blob_response_levels_num; ++curr_lev)
	{		
		const auto& input_cont = eig_ratio_levels[curr_lev];
		auto& output_cont = m_props_hessian_ratio_of_proj[curr_lev];
		CopyVectorToMap(input_cont, Vertices(), output_cont);
	}

	m_hessian_ratio_calculated = true;
}

void SngPtsFinderScaleSpace::GetCoordinateMaps()
{
	coord_3d_map = get(boost::vertex_info_3d, Vertices());
	m_coord_map = GetProxyPropMapVal(coord_3d_map, GetCoord<Vertice>());
	m_norm_map = GetProxyPropMapVal(coord_3d_map, GetNormal<Vertice>());
	coord_3d_map_tr = get(boost::vertex_info_3d, Triangles());
	coord_map_tr = GetProxyPropMapVal(coord_3d_map_tr, GetCoord<MeshTriangle>());
}
void SngPtsFinderScaleSpace::CalcOnlyProps(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
						   const std::vector<cv::Point3i>& triangles, 
						   const std::vector<std::pair<cv::Point3d, double>>& charges, 
						   const std::vector<std::pair<cv::Point3d, double>>& wdv_radii)
{
	m_mesh_keeper.ConstructMesh(vertices, normals, triangles);
	GetCoordinateMaps();
	CalcTangentBasisFast();
	CalcCurvature();
	CalculateAllPotentials(charges, m_mesh_keeper.GetMesh(), m_vertex_charge_map, m_vertex_charge_dir_map);
	CalculateLennardJonesPotentials(wdv_radii, m_mesh_keeper.GetMesh(), m_vertex_lennard_jones_map, m_vertex_lennard_jones_dir_map);
	m_scale_space_calculated = false;
	m_components_blob_response_calculated = false;
	m_hessian_ratio_calculated = false;
}
void SngPtsFinderScaleSpace::CalcSingPtsFromCalculatedProperties(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
	const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges, 
	const std::vector<std::pair<cv::Point3d, double>>& wdv_radii,
	const bool calc_prop_as_average)
{
	CalculateDistanceMaps();
	//CalcTangentBasisConsistent();
	CalcTangentBasisFast();
	m_uncoincided_vertices.SetGraph(Vertices());
	MarkUncoincidedBasises(Vertices(), m_tangent_basis_map, m_uncoincided_vertices);
	CalcSingPtsFromCurvatureScales();
}

void SngPtsFinderScaleSpace::CalcTangentBasisConsistent()
{
	//calculate tangent basis	
	m_tangent_basis_map.SetGraph(Vertices());
	//CalcTangentCoordSystemMap_SystemSolving(Vertices(), m_coord_map, m_norm_map, m_tangent_basis_map);

	CalcTangentCoordSystemMap(Vertices(), m_coord_map, m_norm_map, 0.1, m_tangent_basis_map);
}

void SngPtsFinderScaleSpace::CalcTangentBasisFast()
{
	//calculate tangent basis	
	m_tangent_basis_map.SetGraph(Vertices());
	VectorOrientedCoordFinder coord_finder;
	for (auto curr_vert = boost::vertices(Vertices()).first, end_vert = boost::vertices(Vertices()).second; curr_vert != end_vert; ++curr_vert)
	{
		cv::Mat_<double> curr_normal, new_basis;
		Point_ToMat_Transposed(m_norm_map[*curr_vert], curr_normal);

		coord_finder.Process(curr_normal, new_basis);
		m_tangent_basis_map[*curr_vert] = new_basis.t();//inv(DECOMP_SVD);
	}
}

void SngPtsFinderScaleSpace::GetNonMarkedSingularPointsLevels(std::pair<std::vector<std::vector<NonMarkedSingularPoint>>, 
	std::vector<std::vector<size_t>>>& non_marked_singular_points)
{
	const size_t levels_num = m_maximums_with_levels.size();
	CV_Assert(levels_num == m_sing_pts_levels_num);
	non_marked_singular_points.first.resize(levels_num);
	non_marked_singular_points.second.resize(levels_num);

	for (size_t level = 0; level <m_sing_pts_levels_num; ++level)
	{
		const size_t points_num = m_maximums_with_levels[level].size();
		non_marked_singular_points.first[level].resize(points_num);
		non_marked_singular_points.second[level].resize(points_num);
		for (size_t ind = 0; ind < points_num; ++ind)
		{
			auto& curr_point = non_marked_singular_points.first[level][ind];
			const VertexDescriptor curr_descr = m_maximums_with_levels[level][ind];
			curr_point.Property().SurfaceType() = m_vertex_curv_type_mesh_levels[level
 + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls][curr_descr];
			curr_point.Property().Charge() = Sign(m_vertex_charge_map[curr_descr]);
			curr_point.Property().ElectricPotential() = m_output_scale_space[
				level+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls][ELECTR_POTENT][curr_descr]
			* m_mean_and_sigma[ELECTR_POTENT][1]+ m_mean_and_sigma[ELECTR_POTENT][0] - 3 * m_mean_and_sigma[ELECTR_POTENT][1];
			
			curr_point.Property().LennardJones() = m_output_scale_space[
				level+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls][STERIC_POTENT][curr_descr]
			* m_mean_and_sigma[STERIC_POTENT][1]+ m_mean_and_sigma[STERIC_POTENT][0] - 3 * m_mean_and_sigma[STERIC_POTENT][1];;//m_vertex_lennard_jones_map[curr_descr];
			curr_point.Property().Area() = 0;
			curr_point.Coord() = m_coord_map[curr_descr];
		}
	}
}

void SngPtsFinderScaleSpace::GetNonMarkedSingularPoints(std::pair<std::vector<NonMarkedSingularPoint>, std::vector<size_t>>& non_marked_singular_points)
{
	const size_t sing_pts_num = m_maximums.size();//num_vertices(m_singular_points_graph);
	non_marked_singular_points.first.resize(sing_pts_num);
	non_marked_singular_points.second.resize(sing_pts_num);

	for (size_t ind = 0; ind < m_maximums.size(); ++ind)
	{
		auto& curr_point = non_marked_singular_points.first[ind];
		const VertexDescriptor curr_descr = m_maximums[ind];
		curr_point.Property().SurfaceType() = m_vertex_curv_type[curr_descr];
		curr_point.Property().Charge() = Sign(m_vertex_charge_map[curr_descr]);
		curr_point.Property().ElectricPotential() = m_vertex_charge_map[curr_descr];
		curr_point.Property().LennardJones() = m_vertex_lennard_jones_map[curr_descr];
		curr_point.Property().Area() = 0;
		curr_point.Coord() = m_coord_map[curr_descr];
	}
}
void SngPtsFinderScaleSpace::GetVerticesWithDblProp(std::vector<std::pair<cv::Point3d, double>>& vertices_with_prop, 
	const ISingularPointsFinder::SurfProperty prop_type)
{
	vertices_with_prop.clear();
	vertices_with_prop.reserve(num_vertices(Vertices()));

	const DoubleVertGraphProp& prop_map = (prop_type <= ISingularPointsFinder::SURF_PROPS_NUM) ? 
		m_input_prop_map[prop_type] : m_uncoincided_vertices; 
	for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		const double curr_prop = prop_map[*curr_vertice];
		const cv::Point3d curr_coord = m_coord_map[*curr_vertice];
		vertices_with_prop.push_back(std::make_pair(curr_coord, curr_prop));
	}
}

void SngPtsFinderScaleSpace::GetVerticesWithDblPropLevels(std::vector<std::vector<std::pair<cv::Point3d, double>>>& vertices_with_props_lev, 
	const ISingularPointsFinder::SurfProperty prop_type)
{
	const size_t levels_num = 
		(prop_type <= ISingularPointsFinder::SURF_PROPS_NUM) ?
	m_scale_space_levels_num : m_blob_response_levels_num;

	vertices_with_props_lev.clear();
	vertices_with_props_lev.resize(levels_num);

	for (size_t lev = 0; lev < levels_num; ++lev)
	{
		vertices_with_props_lev[lev].reserve(num_vertices(Vertices()));
		const DoubleVertGraphProp& prop_map = (prop_type <= ISingularPointsFinder::SURF_PROPS_NUM) ? 
			m_output_scale_space[lev][prop_type] :

		(prop_type >= ISingularPointsFinder::FIRST_LOG_PROP && prop_type <= ISingularPointsFinder::LAST_LOG_PROP ? 
			m_components_blob_response[lev]
		[prop_type - ISingularPointsFinder::SURF_PROPS_NUM - 1] : 

		(prop_type >= ISingularPointsFinder::FIRST_PCA_PROP && prop_type <= ISingularPointsFinder::LAST_PCA_PROP ?
			m_projecters_coords[lev]
		[prop_type - ISingularPointsFinder::FIRST_PCA_PROP] : 
		
		(prop_type >= ISingularPointsFinder::FIRST_EIG_PROP && prop_type <= ISingularPointsFinder::LAST_EIG_PROP ?
			m_props_hessian_ratio[lev]
		[prop_type - ISingularPointsFinder::FIRST_EIG_PROP] : 

		(prop_type == ISingularPointsFinder::PCA_LOG ? 
			m_blob_response[lev] :

		(prop_type == ISingularPointsFinder::PCA_GRAD ?
			m_projected_grad_norm[lev] :

		(prop_type == ISingularPointsFinder::PCA_EIG ?
			m_props_hessian_ratio_of_proj[lev] :

		(prop_type == ISingularPointsFinder::PCA_EIG_LOG ?
			m_props_hessian_ratio_of_proj_LOG[lev] :

		m_scale_space_projected[lev])))))));

		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			vertices_with_props_lev[lev].push_back(std::make_pair(m_coord_map[*curr_vertice], prop_map[*curr_vertice]));
		}
	}
}
void SngPtsFinderScaleSpace::GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types)
{
	vertices_with_types.clear();
	vertices_with_types.reserve(num_vertices(Vertices()));

	for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		vertices_with_types.push_back(std::make_pair(m_coord_map[*curr_vertice], m_vertex_curv_type[*curr_vertice]));
	}
}

void SngPtsFinderScaleSpace::GetVerticesWithTypesLevels(std::vector<std::vector<std::pair<cv::Point3d, size_t>>>& vertices_with_types_lev)
{
	vertices_with_types_lev.clear();
	vertices_with_types_lev.resize(m_sing_pts_levels_num);

	for (size_t lev = 0; lev < m_sing_pts_levels_num; ++lev)
	{
		vertices_with_types_lev[lev].reserve(num_vertices(Vertices()));

		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			const size_t curr_curv_type = m_vertex_curv_type_mesh_levels[lev + 
				m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls][*curr_vertice];
			vertices_with_types_lev[lev].push_back(std::make_pair(m_coord_map[*curr_vertice], curr_curv_type));
		}
	}
}

void SngPtsFinderScaleSpace::Clear()
{
	m_vertex_charge_map.Clear();
	m_vertex_curv_type.Clear();
}
void SngPtsFinderScaleSpace::CalcCurvature()
{
	//calculate curvature
	m_mean_curvature.SetGraph(Vertices());
	m_gaussian_curvature.SetGraph(Vertices());
	const int kMaxNeighbNum = 10;
	CurvatureCalculator<VerticesGraph> curvature_calculator(kMaxNeighbNum);

	for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		Curvature curvatures;
		curvature_calculator.CalculateCurvatureCubic(Vertices(), get(boost::vertex_info_3d, Vertices()), 
			m_tangent_basis_map, *curr_vertice, curvatures);
		m_gaussian_curvature[*curr_vertice] = curvatures.gaussian_curv;
		m_mean_curvature[*curr_vertice] = curvatures.mean_curv;
	}

	m_vertex_curv_type.SetGraph(Vertices());
	for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		m_vertex_curv_type[*curr_vertice] = GetCurvatureType(m_gaussian_curvature[*curr_vertice], m_mean_curvature[*curr_vertice]);
	}
}

void SngPtsFinderScaleSpace::CalcScaleSpaceHessian()
{
	m_props_hessian_ratio.resize(m_blob_response_levels_num);
	m_grad_dx.resize(m_blob_response_levels_num);
	m_grad_dy.resize(m_blob_response_levels_num);	
	m_hessian_det.resize(m_blob_response_levels_num);
	m_hessian.resize(m_blob_response_levels_num);
	m_scalar_blob_response.resize(m_blob_response_levels_num);
	m_props_hessian_ratio_of_proj.resize(m_blob_response_levels_num);
	m_hessian_map_calculators.resize(m_blob_response_levels_num);
	if (m_use_spherical)
	{
		m_spherical_hessian_map_calculators.resize(m_blob_response_levels_num);
		m_hessian_spherical.resize(m_blob_response_levels_num);
	}

	m_projected_grad_norm.resize(m_blob_response_levels_num);

	//calculate gradients and hessians fro each level of each property
	const int prev_threads_num = omp_get_num_threads();
	omp_set_num_threads(8);
#pragma omp parallel for
	for (int curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
	{
		for (size_t curr_prop = 0; curr_prop < m_output_scale_space[curr_level].size(); ++curr_prop)
		{
			for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
				curr_vertice != end_vertices; ++curr_vertice)
			{
				m_output_scale_space[curr_level][curr_prop][*curr_vertice] = 
					m_output_scale_space[curr_level][curr_prop][*curr_vertice] * m_scale_space_blurrer.GetSigma(curr_level);
			}	
		}

		m_props_hessian_ratio[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_grad_dx[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_grad_dy[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_hessian_det[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_hessian[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		if (m_use_spherical)
		{
			m_hessian_spherical[curr_level].resize(2 *(ISingularPointsFinder::LAST_SPHER_PROP - ISingularPointsFinder::FIRST_SPHER_PROP + 1));
		}
		for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
		{
			m_props_hessian_ratio[curr_level][curr_prop].SetGraph(Vertices());
			m_hessian_det[curr_level][curr_prop].SetGraph(Vertices());
			m_hessian[curr_level][curr_prop].SetGraph(Vertices());

			m_hessian_map_calculators[curr_level].Process(Vertices(), m_coord_map, m_output_scale_space[curr_level][curr_prop],
				m_vert_vert_dist, 1.0, m_tangent_basis_map, true, m_props_hessian_ratio[curr_level][curr_prop]);
			m_grad_dx[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_dx;
			m_grad_dy[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_dy;
			m_hessian_det[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_hess_det;
			m_hessian[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_hessian_map;
		}	
		//Calculate spherical hessians if we use them
		if (m_use_spherical)
		{
			for (int curr_prop = 0; curr_prop < ISingularPointsFinder::LAST_SPHER_PROP - ISingularPointsFinder::FIRST_SPHER_PROP + 1; ++curr_prop)
			{
				m_hessian_spherical[curr_level][2 * curr_prop].SetGraph(Vertices());
				m_hessian_spherical[curr_level][2 * curr_prop + 1].SetGraph(Vertices());
				m_spherical_hessian_map_calculators[curr_level].Process(Vertices(), m_coord_map, m_output_spher_scale_space[curr_level][curr_prop],
					m_vert_vert_dist, 1.0, m_tangent_basis_map);
				for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
					curr_vertice != end_vertices; ++curr_vertice)
				{
					m_spherical_hessian_map_calculators[curr_level].m_hessian_map_x[*curr_vertice] *= m_scale_space_blurrer.GetSigma(curr_level);// rescale the hessian
					m_spherical_hessian_map_calculators[curr_level].m_hessian_map_y[*curr_vertice] *= m_scale_space_blurrer.GetSigma(curr_level);// rescale the hessian
				}
				m_hessian_spherical[curr_level][2 * curr_prop] = m_spherical_hessian_map_calculators[curr_level].m_hessian_map_x;
				m_hessian_spherical[curr_level][2 * curr_prop + 1] = m_spherical_hessian_map_calculators[curr_level].m_hessian_map_y;
			}		
		}
		m_projected_grad_norm[curr_level].SetGraph(Vertices());
		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			m_projected_grad_norm[curr_level][*curr_vertice] = 
				sqrt(Sqr(m_grad_dx[curr_level][ISingularPointsFinder::GAUSS_CURV][*curr_vertice])
				+ Sqr(m_grad_dy[curr_level][ISingularPointsFinder::GAUSS_CURV][*curr_vertice]));
		}
	}	

	std::vector<DoubleVertGraphProp> filtered_det(ISingularPointsFinder::SURF_PROPS_NUM, DoubleVertGraphProp(Vertices()));
	for (int curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
	{
	//post-filter		
		MedianKernelDist<double> med_kernel(0.5);
		FilterMeshWeightedFuncMultiThread(Vertices(), Triangles(), m_vert_vert_dist, m_vert_tr_dist, 
			m_hessian_det[curr_level], true, 8,	med_kernel, filtered_det);
		m_hessian_det[curr_level] = filtered_det;
	}
	
	m_components_blob_response = m_hessian_det;
	omp_set_num_threads(prev_threads_num);	
}

/*
void SngPtsFinderScaleSpace::CalcScaleSpacePropsHessianRatio()
{
	
	m_props_hessian_ratio.resize(m_blob_response_levels_num);
	m_grad_dx.resize(m_blob_response_levels_num);
	m_grad_dy.resize(m_blob_response_levels_num);

	m_projected_grad_x.resize(m_blob_response_levels_num);
	m_projected_grad_y.resize(m_blob_response_levels_num);;
	m_projected_grad_norm.resize(m_blob_response_levels_num);
	m_props_hessian_ratio_of_proj.resize(m_blob_response_levels_num);
	m_hessian_map_calculators.resize(m_blob_response_levels_num);
	m_hessian_det.resize(m_blob_response_levels_num);
	//calculate gradients and hessians fro each level of each property
	omp_set_num_threads(8);
#pragma omp parallel for
	for (int curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
	{
		m_props_hessian_ratio[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_grad_dx[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_grad_dy[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_hessian_det[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);

		for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
		{
			m_props_hessian_ratio[curr_level][curr_prop].SetGraph(Vertices());
			m_hessian_det[curr_level][curr_prop].SetGraph(Vertices());

			m_hessian_map_calculators[curr_level].Process(Vertices(), m_coord_map, m_output_scale_space[curr_level][curr_prop],
				m_vert_vert_dist, 1.0, m_tangent_basis_map, true, m_props_hessian_ratio[curr_level][curr_prop]);
			m_grad_dx[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_dx;
			m_grad_dy[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_dy;
			m_hessian_det[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_hess_det;
		}

		m_projected_grad_x[curr_level].SetGraph(Vertices());
		m_projected_grad_y[curr_level].SetGraph(Vertices());
		m_projected_grad_norm[curr_level].SetGraph(Vertices());
		m_props_hessian_ratio_of_proj[curr_level].SetGraph(Vertices());
		//find projections
		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			m_projected_grad_x[curr_level][*curr_vertice] = 
				ProjectPCADiff(m_grad_dx[curr_level], m_coord_map, m_scale_space_projecter[curr_level], *curr_vertice);
			m_projected_grad_y[curr_level][*curr_vertice] = 
				ProjectPCADiff(m_grad_dy[curr_level], m_coord_map, m_scale_space_projecter[curr_level], *curr_vertice);
			m_projected_grad_norm[curr_level][*curr_vertice] = 
				sqrt(Sqr(m_grad_dx[curr_level][ISingularPointsFinder::GAUSS_CURV][*curr_vertice])
				+ Sqr(m_grad_dy[curr_level][ISingularPointsFinder::GAUSS_CURV][*curr_vertice]));
				//sqrt(Sqr(m_projected_grad_x[curr_level][*curr_vertice]) + Sqr(m_projected_grad_y[curr_level][*curr_vertice]));
		}
		//find hessian of projections
		m_hessian_map_calculators[curr_level].ProcessFromGradient(Vertices(), m_coord_map, 
			m_projected_grad_x[curr_level],m_projected_grad_y[curr_level],
			m_vert_vert_dist, 1.0, m_tangent_basis_map, true, m_props_hessian_ratio_of_proj[curr_level]);

		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			m_projected_grad_norm[curr_level][*curr_vertice] = 0;
			for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
			{
				m_projected_grad_norm[curr_level][*curr_vertice] += m_hessian_det[curr_level][curr_prop][*curr_vertice];	
			}
		}
	}
}
*/
void SngPtsFinderScaleSpace::CalcHessianOfProjectedLog()
{	//find hessian of LOG projections
	m_props_hessian_ratio_of_proj_LOG.resize(m_blob_response_levels_num);

	for (size_t curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
	{
		m_props_hessian_ratio_of_proj_LOG[curr_level].SetGraph(Vertices());
		m_hessian_map_calculators[0].Process(Vertices(), m_coord_map, m_blob_response[curr_level],
			m_vert_vert_dist, 1.0, m_tangent_basis_map, true, m_props_hessian_ratio_of_proj_LOG[curr_level]);
	}
}

void SngPtsFinderScaleSpace::CalculateDistanceMaps()
{
	//mean edge length
	double mean_dist = 0.0;
	double edges_num = 0.0;
	std::vector<double> edges_lengths;
	for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		for (auto curr_neighb = adjacent_vertices(*curr_vertice, Vertices()).first, 
			end_neighb = adjacent_vertices(*curr_vertice, Vertices()).second; 
			curr_neighb != end_neighb; ++curr_neighb)
		{
			const double curr_dist = cv::norm(m_coord_map[*curr_vertice] - m_coord_map[*curr_neighb]);
			edges_lengths.push_back(curr_dist);
			mean_dist += curr_dist;
			++edges_num;
		}
	}	
	std::nth_element(edges_lengths.begin(), edges_lengths.begin() + Round(0.9 * edges_lengths.size()), edges_lengths.end());
	std::cout << mean_dist / edges_num << " " << edges_lengths[Round(0.9 * edges_lengths.size())] << "\n";

	//distance 
	//calculate distances
	m_vert_vert_dist.create(num_vertices(Vertices()), num_vertices(Vertices()));
	m_vert_tr_dist.create(num_vertices(Vertices()), num_vertices(Triangles()));

	if (m_use_euclid_distance)
	{
		CalcDistBetweenGraphs(Vertices(), m_coord_map, Vertices(), m_coord_map, m_vert_vert_dist);
		CalcDistBetweenGraphs(Vertices(), m_coord_map, Triangles(), coord_map_tr, m_vert_tr_dist);
	}
	else
	{//USE GEODESIC DISTANCE
		DijkstraDistMapCalculator<VerticesGraph, double> djikstra_dist;
		djikstra_dist.Calc(Vertices(), m_vert_vert_dist);

		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			for (auto curr_tr = vertices(Triangles()).first, end_tr = vertices(Triangles()).second; 
				curr_tr != end_tr; ++curr_tr)
			{
				m_vert_tr_dist(*curr_vertice, *curr_tr) = (m_vert_vert_dist(*curr_vertice, coord_3d_map_tr[*curr_tr].GetA()) + 
					m_vert_vert_dist(*curr_vertice, coord_3d_map_tr[*curr_tr].GetB()) + 
					m_vert_vert_dist(*curr_vertice, coord_3d_map_tr[*curr_tr].GetC())) / 3.0;
			}
		}
	}

	double min_val = 0, max_val = 0;
	cv::minMaxLoc(m_vert_vert_dist, &min_val, &max_val);
	std::cout << min_val << " " << max_val << "\n";
	cv::minMaxLoc(m_vert_tr_dist, &min_val, &max_val);
	std::cout << min_val << " " << max_val << "\n";

}

void SngPtsFinderScaleSpace::CalcScalarBlobResponse()
{
	for (int curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
	{
		m_scalar_blob_response[curr_level].SetGraph(Vertices());

		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			m_scalar_blob_response[curr_level][*curr_vertice] = 0;

			for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
			{
				m_scalar_blob_response[curr_level][*curr_vertice] += m_hessian_det[curr_level][curr_prop][*curr_vertice];	
			}

			if (m_use_spherical)
			{
				for (int curr_prop = 0; curr_prop < ISingularPointsFinder::LAST_SPHER_PROP - ISingularPointsFinder::FIRST_SPHER_PROP + 1; ++curr_prop)
				{
					for (int i = 0; i < 2; ++i)
					{
						const Mat_<double> curr_hess = m_hessian_spherical[curr_level][2 * curr_prop + i][*curr_vertice];
						const double det = curr_hess(0,0) * curr_hess(1,1) - curr_hess(1,0) * curr_hess(0,1);
						m_scalar_blob_response[curr_level][*curr_vertice] += det;	
					}
				}
			}
		}
	}

	m_blob_response = m_scalar_blob_response;
}

void SngPtsFinderScaleSpace::CalcManifoldDeterminantBlobResponse()
{
	m_manifold_det_blob_response.resize(m_blob_response_levels_num);

	for (int curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
	{
		m_manifold_det_blob_response[curr_level].SetGraph(Vertices());

		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			
			cv::Mat_<double> combined = cv::Mat_<double>::zeros(2, 2);

			for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
			{
				combined += m_hessian[curr_level][curr_prop][*curr_vertice] * m_hessian[curr_level][curr_prop][*curr_vertice];	
			}

			if (m_use_spherical)
			{
				for (int curr_prop = 0; curr_prop < ISingularPointsFinder::LAST_SPHER_PROP - ISingularPointsFinder::FIRST_SPHER_PROP + 1; ++curr_prop)
				{
					for (int i = 0; i < 2; ++i)
					{
						combined += m_hessian_spherical[curr_level][2 * curr_prop + i][*curr_vertice] * m_hessian_spherical[curr_level][2 * curr_prop + i][*curr_vertice];
					}
				}
			}

			m_manifold_det_blob_response[curr_level][*curr_vertice] = sqrt(combined(0,0) * combined(1,1) - combined(1,0) * combined(0,1));
		}
	}

	m_blob_response = m_manifold_det_blob_response;
}
void SngPtsFinderScaleSpace::CalcSingPtsFromCurvatureScales()
{
	//DATA PREPARATION
	m_input_prop_map.resize(ISingularPointsFinder::SURF_PROPS_NUM);
	m_input_prop_map[ISingularPointsFinder::GAUSS_CURV] = m_gaussian_curvature;
	m_input_prop_map[ISingularPointsFinder::MEAN_CURV] = m_mean_curvature;
	m_input_prop_map[ISingularPointsFinder::ELECTR_POTENT] = m_vertex_charge_map;
	m_input_prop_map[ISingularPointsFinder::STERIC_POTENT] = m_vertex_lennard_jones_map;
	if (m_use_spherical)
	{
		m_input_spherical_prop_map.resize(ISingularPointsFinder::LAST_SPHER_PROP - ISingularPointsFinder::FIRST_SPHER_PROP + 1);
		m_input_spherical_prop_map[ISingularPointsFinder::ELECTR_FORCE_DIR - ISingularPointsFinder::FIRST_SPHER_PROP] = m_vertex_charge_dir_map;
		m_input_spherical_prop_map[ISingularPointsFinder::STERIC_FORCE_DIR - ISingularPointsFinder::FIRST_SPHER_PROP] = m_vertex_lennard_jones_dir_map;
	}
	ResccaleInputFunctions();
	//MAKE SCALE-SPACE
	const bool use_post_filter = !m_detect_blobs;
	std::cout << "scale-space" << std::endl;

	if (!m_scale_space_calculated)
	{
		m_scale_space_blurrer.MakeScaleSpace(Vertices(), m_vert_vert_dist, Triangles(), m_vert_tr_dist, m_input_prop_map, 
			m_scale_space_levels_num, use_post_filter, m_output_scale_space);
		if (m_use_spherical)
		{
			m_scale_space_blurrer.MakeScaleSpaceSpherical(Vertices(), m_vert_vert_dist, Triangles(), m_vert_tr_dist, m_input_spherical_prop_map, 
				m_scale_space_levels_num, m_output_spher_scale_space);
		}
		m_scale_space_calculated = true;
	}
	//DETECT TYPES OF VERTICES FOR EACH LEVEL, ACCORDING TO THE CURVATURE	
	m_vertex_curv_type_mesh_levels.assign(m_scale_space_levels_num, VetrticesCurvMap());
	for (size_t curr_level = 0; curr_level < m_scale_space_levels_num; ++curr_level)
	{
		m_vertex_curv_type_mesh_levels[curr_level].SetGraph(Vertices());
		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			m_vertex_curv_type_mesh_levels[curr_level][*curr_vertice] = GetCurvatureType(
				m_output_scale_space[curr_level][GAUSS_CURV][*curr_vertice] * m_mean_and_sigma[GAUSS_CURV][1]+ 
				m_mean_and_sigma[GAUSS_CURV][0] - 3 * m_mean_and_sigma[GAUSS_CURV][1], 
				m_output_scale_space[curr_level][MEAN_CURV][*curr_vertice]* m_mean_and_sigma[MEAN_CURV][1]+ 
				m_mean_and_sigma[MEAN_CURV][0] - 3 * m_mean_and_sigma[MEAN_CURV][1]);
		}	
	}	
	//DETECT FEATURES
	m_scale_space_projecter.resize(m_blob_response_levels_num);

	if (m_detect_blobs)
	{		
		std::cout << "LOG calculation" << std::endl;

		if (!m_components_blob_response_calculated)
		{
			m_components_blob_response.clear();
			m_components_blob_response.resize(m_blob_response_levels_num);

			if (m_detector_type == DOG)
			{
				CalculateDOG();
			}
			else if (m_detector_type == LOG)
			{
				CalculateLOG();
			}
			else if (m_detector_type == HESS_DET)
			{
				CalcScaleSpaceHessian();
				m_hessian_ratio_calculated = true;
			}
			else
			{
				CV_Assert(0);
			}

			m_components_blob_response_calculated = true;
		}

		//END LOG CALCULATION

		if (m_channel_combination == PCA || m_channel_combination == NORM || m_channel_combination == DETECTOR_NORM)
		{
			//FIND VECTORS ON WHICH WE WILL PROJECT
			std::cout << "find PCA" << std::endl;
			m_scale_space_projecter.resize(m_blob_response_levels_num);
			const int prev_threads = omp_get_num_threads();
			omp_set_num_threads(8);
	#pragma omp parallel for
			for (int lev = 0; lev < m_blob_response_levels_num; ++lev)
			{
				m_scale_space_projecter[lev].SetGraph(Vertices());

				if (m_channel_combination == PCA)
				{
					FindVectorsForProjection(Vertices(), m_coord_map, m_output_scale_space[lev], 
						m_vert_vert_dist, 2.0 * m_scale_space_blurrer.GetSigma(lev), m_scale_space_projecter[lev]);
				}
				else if (m_channel_combination == NORM)
				{
					SetPCAAsMeanOfFunction(Vertices(), m_output_scale_space[lev], m_vert_vert_dist, 
						2.0 * m_scale_space_blurrer.GetSigma(lev), m_scale_space_projecter[lev]);
				}
				else if (m_channel_combination == DETECTOR_NORM)
				{
					SetPCAAsMeanOfFunction(Vertices(), m_components_blob_response[lev], m_vert_vert_dist, 
						2.0 * m_scale_space_blurrer.GetSigma(lev), m_scale_space_projecter[lev]);
				}
			}

			omp_set_num_threads(prev_threads);
			CalculateProjectedVectors();

			if (!m_hessian_ratio_calculated)
			{
				std::cout << "hessian" << std::endl;
				//CalcScaleSpacePropsHessianRatio();
				CalcHessianOfProjectedLog();
				m_hessian_ratio_calculated = true;
			}
		}
		else if (m_channel_combination == SCALAR_BLOB_RESPONSE)
		{
			CalcScalarBlobResponse();
		}
		else if (m_channel_combination == MANIFOLD_DET_BLOB_RESPONSE)
		{
			CalcManifoldDeterminantBlobResponse();
		}
		else
		{
			CV_Assert(0);
		}
		//////////////////////////////////////////////////////////////////////////
		//DETECT POINTS
		m_maximums_with_levels.clear();

		if (m_combine_channels)
		{
			std::vector<DoubleVertGraphProp> input(1, DoubleVertGraphProp(Vertices()));
			std::vector<DoubleVertGraphProp> filtered_detector(1, DoubleVertGraphProp(Vertices()));

			for (size_t curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
			{
				input[0] = m_blob_response[curr_level];
				/*MedianKernelWeightedDist*/MedianKernelDist<double> med_kernel(0.5);
				//GaussianKernelWeightedDistTable<double> av_kernel(curr_sigma, m_exp_approx);
				FilterMeshWeightedFuncMultiThread(Vertices(), Triangles(), m_vert_vert_dist, m_vert_tr_dist, 
					input, true, 8,	med_kernel, filtered_detector);
				//m_detector_function[curr_level] = filtered_detector[0];
				const DoubleVertGraphProp& prop_map_to_process = m_blob_response[curr_level];//filtered_detector[0];
				double mean = 0.0;
				for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
					curr_vertice != end_vertices; ++curr_vertice)
				{
					mean += prop_map_to_process[*curr_vertice];
				}
				mean /= num_vertices(Vertices());
				double dev = 0.0;
				for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
					curr_vertice != end_vertices; ++curr_vertice)
				{
					dev += Sqr(mean - prop_map_to_process[*curr_vertice]);
				}
				dev /= num_vertices(Vertices());
				std::cout << dev << " ";
				for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
					curr_vertice != end_vertices; ++curr_vertice)
				{
					m_blob_response[curr_level][*curr_vertice] = prop_map_to_process[*curr_vertice] / sqrt(dev);
				}
			}

			FindSingPtsAsCombinedMaximumsOfLOG();
		}
		else
		{
			FindSingPtsAsSeparateMaximumsOfLOG();
		}
	}
	else
	{
		FindSingPtsAsMaximumsOfScaleSpace();
	}
	///////
	PointsKeeper<VerticesGraph> points_keeper(1.0, Vertices());

	for (size_t curr_level = 0; curr_level < m_sing_pts_levels_num; ++curr_level)
	{
		std::cout << m_maximums_with_levels[curr_level].size() << " ";	

		points_keeper.AddPoints(m_maximums_with_levels[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
			+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
	}
	m_maximums = points_keeper.Points();

	std::cout << "maximums " << m_maximums.size() << " ";	
}

void SngPtsFinderScaleSpace::ResccaleInputFunctions()
{
	if (!m_mean_and_sigma.empty())
	{//rescale functions
		for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++ curr_prop)
		{
			const double min_val = m_mean_and_sigma[curr_prop][0] - 3.0 * m_mean_and_sigma[curr_prop][1];
			const double max_val = m_mean_and_sigma[curr_prop][0] + 3.0 * m_mean_and_sigma[curr_prop][1];

			for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
				curr_vertice != end_vertices; ++curr_vertice)
			{
				m_input_prop_map[curr_prop][*curr_vertice] = std::min(m_input_prop_map[curr_prop][*curr_vertice], max_val);
				m_input_prop_map[curr_prop][*curr_vertice] = std::max(m_input_prop_map[curr_prop][*curr_vertice], min_val);
				m_input_prop_map[curr_prop][*curr_vertice] -= min_val;
				m_input_prop_map[curr_prop][*curr_vertice] /= m_mean_and_sigma[curr_prop][1];
			}
		}
	}
}

void SngPtsFinderScaleSpace::CalculateLOG()
{
	for (size_t curr_level = 0; curr_level <m_blob_response_levels_num; ++curr_level)
	{
		const double coeff = m_scale_space_blurrer.GetSigma(curr_level);
		m_components_blob_response[curr_level].resize(m_output_scale_space[curr_level].size());

		for (size_t curr_prop = 0; curr_prop < m_components_blob_response[curr_level].size(); ++curr_prop)
		{
			m_components_blob_response[curr_level][curr_prop].SetGraph(Vertices());
			/*SimpleLaplacian(vertices_graph, m_output_scale_space[curr_level][curr_prop], 
				output_scale_space_diff[curr_level][curr_prop]);*/
		}

		//LaplaceBeltramiKernelWeightedDist<double, SignedDistFunc<double>> kernel(0.5/*(m_scale_space_blurrer.GetSigma(curr_level) + 1.7) / 4.0*//*0.3*/, SignedDistFunc<double>());
		LoGKernelWeightedDist<double> kernel(m_scale_space_blurrer.GetSigma(curr_level));
		FilterMeshWeightedFunc(Vertices(), Triangles(), 
			m_vert_vert_dist, m_vert_tr_dist, m_input_prop_map, true,
			kernel, m_components_blob_response[curr_level]);
		/*FilterMeshWeightedFunc(vertices_graph, triangles_graph, 
		vert_vert_dist, vert_tr_dist, m_output_scale_space[curr_level], false,
		kernel, output_scale_space_diff[curr_level]);*/
					
		for (size_t curr_prop = 0; curr_prop < m_components_blob_response[curr_level].size(); ++curr_prop)
		{
			for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
				curr_vertice != end_vertices; ++curr_vertice)
			{
				m_components_blob_response[curr_level][curr_prop][*curr_vertice] = m_components_blob_response[curr_level][curr_prop][*curr_vertice] * coeff;
			}	
		}
	}		
}

void SngPtsFinderScaleSpace::CalculateDOG()
{
	ScaleSpaceBlurrer<VerticesGraph, CoordMap, GaussianKernel<cv::Point3d, double>> scale_space_blurrer_shifted;//for calculation of derivative in t
	const double d_mult_sigma = pow(m_sigma_max / m_init_curv_sigma, 1.0 / (m_sing_pts_levels_num));
	scale_space_blurrer_shifted.Init(m_init_curv_sigma, d_mult_sigma, 0.1, false, false, 0.1);
	std::vector<std::vector<DoubleVertGraphProp>> output_scale_space_shifted;
	const bool use_post_filter = !m_detect_blobs;
	scale_space_blurrer_shifted.MakeScaleSpace(Vertices(), m_vert_vert_dist, Triangles(), m_vert_tr_dist, m_input_prop_map, 
		m_scale_space_levels_num, use_post_filter, output_scale_space_shifted);
	//Calculate scales difference which approximates laplace operator
	for (size_t curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
	{
		const double coeff = (scale_space_blurrer_shifted.GetSigma(curr_level) + m_scale_space_blurrer.GetSigma(curr_level)) / 2.0 
			/ (scale_space_blurrer_shifted.GetSigma(curr_level) - m_scale_space_blurrer.GetSigma(curr_level));

		/*const double coeff = (m_scale_space_blurrer.GetSigma(curr_level + 1) + m_scale_space_blurrer.GetSigma(curr_level)) / 2.0 
			/ (m_scale_space_blurrer.GetSigma(curr_level + 1) - m_scale_space_blurrer.GetSigma(curr_level));*/

		m_components_blob_response[curr_level].resize(m_output_scale_space[curr_level].size());
		const auto& shifted_scale_space_lev = /*m_output_scale_space[curr_level + 1]*/output_scale_space_shifted[curr_level];
		for (size_t curr_prop = 0; curr_prop < m_components_blob_response[curr_level].size(); ++curr_prop)
		{
			m_components_blob_response[curr_level][curr_prop].SetGraph(Vertices());
			for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
				curr_vertice != end_vertices; ++curr_vertice)
			{
				m_components_blob_response[curr_level][curr_prop][*curr_vertice] = coeff * 				 
					(shifted_scale_space_lev[curr_prop][*curr_vertice] - m_output_scale_space[curr_level][curr_prop][*curr_vertice]);

			}	
		}
	}
}

void SngPtsFinderScaleSpace::CalculateProjectedVectors()
{
	m_blob_response.clear();
	m_blob_response.resize(m_blob_response_levels_num);
	m_scale_space_projected.clear();
	m_scale_space_projected.resize(m_blob_response_levels_num);
	m_projecters_coords.clear();
	m_projecters_coords.resize(m_blob_response_levels_num);

	for (int lev = 0; lev < m_blob_response_levels_num; ++lev)
	{
		m_blob_response[lev].SetGraph(Vertices());
		m_scale_space_projected[lev].SetGraph(Vertices());
		const size_t kPcaPropsNum = ISingularPointsFinder::LAST_PCA_PROP - ISingularPointsFinder::FIRST_PCA_PROP + 1;
		m_projecters_coords[lev].resize(kPcaPropsNum);
		for (size_t prop_num = 0; prop_num < kPcaPropsNum; ++prop_num)
		{
			m_projecters_coords[lev][prop_num].SetGraph(Vertices());
		}
		for (auto curr_vertice = vertices(Vertices()).first, end_vertices = vertices(Vertices()).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			m_blob_response[lev][*curr_vertice] = 
				ProjectPCADiff(m_components_blob_response[lev], m_coord_map, m_scale_space_projecter[lev], *curr_vertice);
			m_scale_space_projected[lev][*curr_vertice] = 
				ProjectPCA(m_output_scale_space[lev], m_coord_map, m_scale_space_projecter[lev], *curr_vertice);

			for (size_t prop_num = 0; prop_num < kPcaPropsNum; ++prop_num)
			{
				m_projecters_coords[lev][prop_num][*curr_vertice]
				= std::abs(m_scale_space_projecter[lev][*curr_vertice].pca.eigenvectors.at<double>(3, 3 + prop_num));
			}
		}			
	}
}

void SngPtsFinderScaleSpace::FindSingPtsAsCombinedMaximumsOfLOG()
{
	m_maximums_with_levels.clear();
	const double kMaxRadius = 1.0;
	if (!m_use_central_projector)
	{
		FindLocalMaximumsOnLevels(Vertices(),m_blob_response,  
			 std::greater<double>(), std::less<double>(), m_vert_vert_dist,kMaxRadius, 
			 m_scale_extr, m_maximums_with_levels);
	}
	else
	{
		FindLocalMaximumsOnLevelsVect(Vertices(),m_coord_map, m_components_blob_response, m_scale_space_projecter,  
			m_vert_vert_dist, kMaxRadius, std::greater_equal<double>(), std::less_equal<double>(), m_scale_extr, 
			m_use_central_projector, m_maximums_with_levels);
	}

	for (size_t curr_level = 0; curr_level < m_maximums_with_levels.size(); ++curr_level)
	{
		std::cout << m_maximums_with_levels[curr_level].size() << " ";
	}
	std::cout << "\n";
	if (m_filter_by_eigenvalues)
	{
		FilterByEigenvaluesRatio();
	}
	/*const auto maximums = m_maximums_with_levels;
	FilterByEigenvaluesRatio();
	m_maximums_with_levels = maximums;*/
	//FilterByNormedLogValue();
	//m_maximums_with_levels = maximums;
	m_maximums_with_levels.resize(m_sing_pts_levels_num);
	for (size_t curr_level = 1; curr_level < m_maximums_with_levels.size(); ++curr_level)
	{
		/*m_maximums_with_levels[0].insert(m_maximums_with_levels[0].end(), 
			m_maximums_with_levels[curr_level].begin(), m_maximums_with_levels[curr_level].end());*/
	}		
}

void SngPtsFinderScaleSpace::FilterByNormedLogValue()
{
	//find detector function max value
	//max value hass usually small absolute value, so it can be used as a threshold for detector function valueus.
	//maximums within this threshold usually concide to edges in detector function, so must be discarded
	std::vector<double> vect_for_thresh_search;

	for (size_t curr_level = 0; curr_level < m_maximums_with_levels.size(); ++curr_level)
	{
		//find detector function threshold
		const auto& curr_detect_func_values = m_blob_response[curr_level + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls];
		//vect_for_thresh_search.clear();

		//for (auto curr_val = curr_detect_func_values.begin(), end_val = curr_detect_func_values.end(); curr_val != end_val; ++curr_val)
		//{
		//	if (*curr_val > 0.0)
		//	{
		//		vect_for_thresh_search.push_back(*curr_val);
		//	}
		//}

		//double detect_func_thresh = 0.0;

		//if (!vect_for_thresh_search.empty())
		//{
		//	const size_t thresh_ind = vect_for_thresh_search.size() * 0.90;
		//	std::nth_element(vect_for_thresh_search.begin(), vect_for_thresh_search.begin() + thresh_ind, vect_for_thresh_search.end());
		//	detect_func_thresh = vect_for_thresh_search[thresh_ind];
		//}	
		//detect_func_thresh = 0.1;//std::min(0.4, detect_func_thresh);
		////std::cout << detect_func_thresh << "\n";
		const double thresh_low = 0.15;
		const double thresh_high = 0.35;
		std::vector<VertexDescriptor> not_small_detect_func_vals;

		for (int ind1 = 0; ind1 < m_maximums_with_levels[curr_level].size();++ind1)
		{
			const double curr_detect_val = curr_detect_func_values[m_maximums_with_levels[curr_level][ind1]];
			bool good_point = true;
			if (abs(curr_detect_val) < thresh_low)
			{
				good_point = false;
			}
			else if (abs(curr_detect_val) < thresh_high)
			{
				std::vector<VertexDescriptor> nearest_vertices;
				GetVerticesWithinDistPlusAdjacent(m_maximums_with_levels[curr_level][ind1], Vertices(), m_vert_vert_dist, 0.4, nearest_vertices);

				for (auto curr_neighb = nearest_vertices.begin(), end_neighb = nearest_vertices.end(); curr_neighb != end_neighb; ++curr_neighb)
				{
					if (abs(curr_detect_func_values[*curr_neighb]) < thresh_low)
					{
						good_point = false;
						break;
					}
				}
			}

			if (good_point)
			{
				not_small_detect_func_vals.push_back(m_maximums_with_levels[curr_level][ind1]);
			}
		}
		m_maximums_with_levels[curr_level] = not_small_detect_func_vals;
	}
	for (size_t curr_level = 0; curr_level < m_maximums_with_levels.size(); ++curr_level)
	{
		std::cout << m_maximums_with_levels[curr_level].size() << " ";
	}
	std::cout << "\n";
}

void SngPtsFinderScaleSpace::FilterByEigenvaluesRatio()
{	
	for (size_t curr_level = 0; curr_level < m_maximums_with_levels.size(); ++curr_level)
	{
		std::vector<VertexDescriptor> not_near_ridge;

		for (int ind1 = 0; ind1 < m_maximums_with_levels[curr_level].size();++ind1)
		{
			bool is_near = false;

			std::vector<VertexDescriptor> nearest_vertices;
			GetVerticesWithinDistPlusAdjacent(m_maximums_with_levels[curr_level][ind1], Vertices(), m_vert_vert_dist, 0.4, nearest_vertices);

			for (auto neighb_it = nearest_vertices.begin(),	end_neighb = nearest_vertices.end(); neighb_it != end_neighb; ++neighb_it)					
			{
				if (m_props_hessian_ratio_of_proj[curr_level + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls][*neighb_it] >= m_ratio_thresh)
				{
					is_near = true;
					break;
				}	
			}
			if (!is_near)
			{
				not_near_ridge.push_back(m_maximums_with_levels[curr_level][ind1]);
			}
		}
		m_maximums_with_levels[curr_level] = not_near_ridge;
	}
	for (size_t curr_level = 0; curr_level < m_maximums_with_levels.size(); ++curr_level)
	{
		std::cout << m_maximums_with_levels[curr_level].size() << " ";
	}
	std::cout << "\n";
}

void SngPtsFinderScaleSpace::FindSingPtsAsSeparateMaximumsOfLOG()
{
	m_maximums_with_levels.clear();
	//copy in another container, because we need another order in detection function
	std::vector<std::vector<DoubleVertGraphProp>> output_scale_space_diff_another_order(ISingularPointsFinder::SURF_PROPS_NUM);
	for (int curr_prop = 0; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
	{
		output_scale_space_diff_another_order[curr_prop].resize(m_blob_response_levels_num);
		for (size_t curr_level = 0; curr_level < m_blob_response_levels_num; ++curr_level)
		{
			output_scale_space_diff_another_order[curr_prop][curr_level] = m_components_blob_response[curr_level][curr_prop];
		}
	}
	m_vertex_curv_type_mesh_levels.resize(m_blob_response_levels_num);
	//find points
	std::vector<std::vector<std::vector<VertexDescriptor>>> sing_points(ISingularPointsFinder::SURF_PROPS_NUM);//Prop -> level -> points set

	for (int curr_prop = 0; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
	{
		FindLocalMaximumsOnLevels(Vertices(),output_scale_space_diff_another_order[curr_prop],  
			std::greater_equal<double>(), std::less_equal<double>(), m_vert_vert_dist, 0.5, 
			true, sing_points[curr_prop]);
	}

	for (size_t curr_level = 0; curr_level < m_sing_pts_levels_num; ++curr_level)
	{
		PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, Vertices());

		for (int curr_prop = 0; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
		{
			points_keeper_curr_scale.AddPoints(sing_points[curr_prop][curr_level], m_vertex_curv_type_mesh_levels[curr_level 
				+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
		}
		m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
		std::cout <<points_keeper_curr_scale.Points().size() << "\n";
	}
}

void SngPtsFinderScaleSpace::FindSingPtsAsMaximumsOfScaleSpace()
{
	m_maximums_with_levels.clear();
	//copy in another container, because we need another order in detection function
	std::vector<std::vector<DoubleVertGraphProp>> output_scale_space_another_order(ISingularPointsFinder::SURF_PROPS_NUM);
	for (int curr_prop = 0; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
	{
		output_scale_space_another_order[curr_prop].resize(m_scale_space_levels_num);
		for (size_t curr_level = 0; curr_level < m_scale_space_levels_num; ++curr_level)
		{
			output_scale_space_another_order[curr_prop][curr_level] = m_output_scale_space[curr_level][curr_prop];
		}
	}
	//find points
	std::vector<std::vector<std::vector<VertexDescriptor>>> sing_points(ISingularPointsFinder::SURF_PROPS_NUM);//Prop -> level -> points set

	for (int curr_prop = 0; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
	{
		FindScaleSingularPointsOnFunc(Vertices(), m_coord_map, output_scale_space_another_order[curr_prop],
			m_vertex_curv_type_mesh_levels, 1.0, sing_points[curr_prop], true);
	}

	for (size_t curr_level = 0; curr_level < m_sing_pts_levels_num; ++curr_level)
	{
		PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, Vertices());
		for (int curr_prop = 0; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
		{
			points_keeper_curr_scale.AddPoints(sing_points[curr_prop][curr_level], m_vertex_curv_type_mesh_levels[curr_level 
				+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
		}
		m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
		std::cout <<points_keeper_curr_scale.Points().size() << "\n";
	}
}

std::vector<double> SngPtsFinderScaleSpace::GetSigmaValues()
{
	std::vector<double> sing_pts_sigmas(m_sing_pts_levels_num);
	for (int i = 0; i < m_sing_pts_levels_num; ++i)
	{
		sing_pts_sigmas[i] = m_scale_space_blurrer.GetSigma(i + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls);
	}

	return sing_pts_sigmas;
}

template <typename Iter>
void CalcMeanAndDev(const Iter it_beg, const Iter it_end, double& mean, double& dev)
{
	mean = 0;
	int length = 0;
	for (Iter iter = it_beg; iter != it_end; ++iter)
	{
		mean += *iter;
		++length;
	}

	if (length == 0)
	{
		return;
	}
	mean /= static_cast<double>(length);

	dev = 0;
	for (Iter iter = it_beg; iter != it_end; ++iter)
	{
		dev += abs(*iter - mean);
	}
	dev /= static_cast<double>(length);
}

}