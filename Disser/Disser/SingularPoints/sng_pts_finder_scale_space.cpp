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

void SngPtsFinderScaleSpace::InitParams(int argc, char** argv)
{
	ReadParamFromCommandLineWithDefault(argc, argv, "-mesh_levels_num", m_sing_pts_levels_num, 10);
	ReadParamFromCommandLineWithDefault(argc, argv, "-detect_blobs", m_detect_blobs, true);
	ReadParamFromCommandLineWithDefault(argc, argv, "-use_DOG_as_LOG_approximation", m_use_DOG_as_LOG_approximation, true);
	ReadParamFromCommandLineWithDefault(argc, argv, "-use_euclid_distance", m_use_euclid_distance, false);
	ReadParamFromCommandLineWithDefault(argc, argv, "-combine_channels", m_combine_channels, true);
	ReadParamFromCommandLineWithDefault(argc, argv, "-one_ring_neighb", m_one_ring_neighb, true);

	ReadParamFromCommandLineWithDefault(argc, argv, "-init_curv_sigma", m_init_curv_sigma, 0.33);
	ReadParamFromCommandLineWithDefault(argc, argv, "-sigma_max", m_sigma_max, 2.89);

	m_scale_space_levels_num = m_sing_pts_levels_num + (m_detect_blobs ? 3 : 2);
	m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls = 1;

	//INIT SCALE-SPACE
	const bool kAdditive = false;
	const double d_mult_sigma = pow(m_sigma_max / m_init_curv_sigma, 1.0 / (m_sing_pts_levels_num));
	std::cout << "d_mult_sigma " << d_mult_sigma << std::endl;
	m_scale_space_blurrer.Init(m_init_curv_sigma, d_mult_sigma, 0.0, kAdditive);
}

void SngPtsFinderScaleSpace::CalcOnlyProps(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
						   const std::vector<cv::Point3i>& triangles, 
						   const std::vector<std::pair<cv::Point3d, double>>& charges, 
						   const std::vector<std::pair<cv::Point3d, double>>& wdv_radii)
{
	m_mesh_keeper.ConstructMesh(vertices, normals, triangles);
	const Mesh& mesh_to_use = /*m_mesh_keeper.GetMesh();*/GetMesh();
	CalcTangentBasis(mesh_to_use);
	CalcCurvature(mesh_to_use);
	CalculateAllPotentials(charges, m_mesh_keeper.GetMesh()/*mesh_to_use*/, m_vertex_charge_map);
	CalculateLennardJonesPotentials(wdv_radii, m_mesh_keeper.GetMesh()/*mesh_to_use*/, m_vertex_lennard_jones_map);
}
void SngPtsFinderScaleSpace::Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
	const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges, 
	const std::vector<std::pair<cv::Point3d, double>>& wdv_radii,
	const bool calc_prop_as_average)
{
	m_mesh_keeper.ConstructMesh(vertices, normals, triangles);
	//FilterMesh(m_mesh_keeper.GetMesh(), GaussianKernel<cv::Point3d, cv::Point3d>(1.0), m_filtered_mesh);
	/*m_filtered_mesh_levels.resize(m);
	m_filtered_mesh_levels_2.resize(mesh_levels_num);*/
	//m_filtered_mesh_levels[0] = m_mesh_keeper.GetMesh();
	/*const double kInitSigma = 0.3;
	const double kDSigma = 0.1;*/
	/*for (size_t curr_level = 0; curr_level < mesh_levels_num; ++curr_level)
	{
		FilterMesh(m_mesh_keeper.GetMesh(), GaussianKernel<cv::Point3d, cv::Point3d>(kInitSigma + (kDSigma) *(curr_level)), 
			m_filtered_mesh_levels[curr_level]);
		FilterMesh(m_mesh_keeper.GetMesh(), GaussianKernel<cv::Point3d, cv::Point3d>(2.0 * (kInitSigma + (kDSigma) *(curr_level))), 
			m_filtered_mesh_levels_2[curr_level]);
	}*/
	const Mesh& mesh_to_use = /*m_mesh_keeper.GetMesh();*/GetMesh();
	CalculateDistanceMaps(mesh_to_use);
	//CalculateVerticesSurfaceType(mesh_to_use);
	CalcTangentBasis(mesh_to_use);
	CalcCurvature(mesh_to_use);
	CalculateAllPotentials(charges, m_mesh_keeper.GetMesh()/*mesh_to_use*/, m_vertex_charge_map);
	CalculateLennardJonesPotentials(wdv_radii, m_mesh_keeper.GetMesh()/*mesh_to_use*/, m_vertex_lennard_jones_map);
	CalcSingPtsFromCurvatureScales(m_mesh_keeper.GetMesh());
	//CalcShiftsMaximums();
	const int kMaxSegmSize = 500;
}

void SngPtsFinderScaleSpace::CalcTangentBasis(const Mesh& mesh)
{
	//calculate tangent basis
	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
		get(boost::vertex_info_3d, mesh.vertices);
	CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());
	NormalMap norm_map = GetProxyPropMap(coord_3d_map, GetNormal<Vertice>());
	m_tangent_basis_map.SetGraph(mesh.vertices);
	CalcTangentCoordSystemMap(mesh.vertices, coord_map, norm_map, m_tangent_basis_map);
}

void SngPtsFinderScaleSpace::GetNonMarkedSingularPointsLevels(std::pair<std::vector<std::vector<NonMarkedSingularPoint>>, 
	std::vector<std::vector<size_t>>>& non_marked_singular_points)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;
	const size_t levels_num = m_maximums_with_levels.size();
	CV_Assert(levels_num == m_sing_pts_levels_num);
	non_marked_singular_points.first.resize(levels_num);
	non_marked_singular_points.second.resize(levels_num);
	auto coord_map_all = get(boost::vertex_info_3d, vertices_graph);
	const auto coord_map = GetProxyPropMap(coord_map_all, GetCoord<Vertice>());
	AverageFinder<VerticesGraph, decltype(coord_map)> av_finder;

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
			curr_point.Coord() = get(boost::vertex_info_3d, vertices_graph, curr_descr).Center();
		}
	}
}

void SngPtsFinderScaleSpace::GetNonMarkedSingularPoints(std::pair<std::vector<NonMarkedSingularPoint>, std::vector<size_t>>& non_marked_singular_points)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;
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
		curr_point.Coord() = get(boost::vertex_info_3d, vertices_graph, curr_descr).Center();
	}
}
void SngPtsFinderScaleSpace::GetVerticesWithDblProp(std::vector<std::pair<cv::Point3d, double>>& vertices_with_prop, 
	const ISingularPointsFinder::SurfProperty prop_type)
{
	const VerticesGraph& vertices_graph = GetMesh().vertices;
	vertices_with_prop.clear();
	vertices_with_prop.reserve(num_vertices(vertices_graph));

	const DoubleVertGraphProp& prop_map = /*(prop_type <= ISingularPointsFinder::SURF_PROPS_NUM) ? */
		m_input_prop_map[prop_type]; 
	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		const double curr_prop = prop_map[*curr_vertice];
		const cv::Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *curr_vertice).Center();
		vertices_with_prop.push_back(std::make_pair(curr_coord, curr_prop));
	}
}

void SngPtsFinderScaleSpace::GetVerticesWithDblPropLevels(std::vector<std::vector<std::pair<cv::Point3d, double>>>& vertices_with_props_lev, 
	const ISingularPointsFinder::SurfProperty prop_type)
{
	const VerticesGraph& vertices_graph = GetMesh().vertices;
	vertices_with_props_lev.clear();
	vertices_with_props_lev.resize(m_sing_pts_levels_num);

	for (size_t lev = 0; lev < m_sing_pts_levels_num; ++lev)
	{
		vertices_with_props_lev[lev].reserve(num_vertices(vertices_graph));
		const DoubleVertGraphProp& prop_map = (prop_type <= ISingularPointsFinder::SURF_PROPS_NUM) ? 
			m_output_scale_space[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls][prop_type] :

		(prop_type >= ISingularPointsFinder::FIRST_LOG_PROP && prop_type <= ISingularPointsFinder::LAST_LOG_PROP ? 
			m_output_scale_space_diff[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]
		[prop_type - ISingularPointsFinder::SURF_PROPS_NUM - 1] : 

		(prop_type >= ISingularPointsFinder::FIRST_PCA_PROP && prop_type <= ISingularPointsFinder::LAST_PCA_PROP ?
			m_projecters_coords[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]
		[prop_type - ISingularPointsFinder::FIRST_PCA_PROP] : 
		
		(prop_type >= ISingularPointsFinder::FIRST_EIG_PROP && prop_type <= ISingularPointsFinder::LAST_EIG_PROP ?
			m_props_hessian_ratio[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]
		[prop_type - ISingularPointsFinder::FIRST_EIG_PROP] : 

		(prop_type == ISingularPointsFinder::PCA_LOG ? 
			m_detector_function_projected[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls] :

		(prop_type == ISingularPointsFinder::PCA_GRAD ?
			m_projected_grad_norm[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls] :

		(prop_type == ISingularPointsFinder::PCA_EIG ?
			m_props_hessian_ratio_of_proj[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls] :

		(prop_type == ISingularPointsFinder::PCA_EIG_LOG ?
			m_props_hessian_ratio_of_proj_LOG[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls] :

		m_scale_space_projected[lev + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls])))))));

		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			const double curr_prop = prop_map[*curr_vertice];
			const cv::Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *curr_vertice).Center();
			vertices_with_props_lev[lev].push_back(std::make_pair(curr_coord, curr_prop));
		}
	}
}
void SngPtsFinderScaleSpace::GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types)
{
	const VerticesGraph& vertices_graph = GetMesh().vertices;
	vertices_with_types.clear();
	vertices_with_types.reserve(num_vertices(vertices_graph));

	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		const size_t curr_curv_type = m_vertex_curv_type[*curr_vertice];
		const cv::Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *curr_vertice).Center();
		vertices_with_types.push_back(std::make_pair(curr_coord, curr_curv_type));
	}
}

void SngPtsFinderScaleSpace::GetVerticesWithTypesLevels(std::vector<std::vector<std::pair<cv::Point3d, size_t>>>& vertices_with_types_lev)
{
	const VerticesGraph& vertices_graph = GetMesh().vertices;
	vertices_with_types_lev.clear();
	vertices_with_types_lev.resize(m_sing_pts_levels_num);

	for (size_t lev = 0; lev < m_sing_pts_levels_num; ++lev)
	{
		vertices_with_types_lev[lev].reserve(num_vertices(vertices_graph));

		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			const size_t curr_curv_type = m_vertex_curv_type_mesh_levels[lev + 
				m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls][*curr_vertice];
			const cv::Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *curr_vertice).Center();
			vertices_with_types_lev[lev].push_back(std::make_pair(curr_coord, curr_curv_type));
		}
	}
}

void SngPtsFinderScaleSpace::Clear()
{
	m_vertex_charge_map.Clear();
	m_vertex_curv_type.Clear();
}
void SngPtsFinderScaleSpace::CalcCurvature(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	//calculate curvature
	m_mean_curvature.SetGraph(vertices_graph);
	m_gaussian_curvature.SetGraph(vertices_graph);
	const int kMaxNeighbNum = 10;
	CurvatureCalculator<VerticesGraph> curvature_calculator(kMaxNeighbNum);

	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		Curvature curvatures;
		curvature_calculator.CalculateCurvatureCubic(vertices_graph, get(boost::vertex_info_3d, vertices_graph), 
			m_tangent_basis_map, *curr_vertice, curvatures);
		m_gaussian_curvature[*curr_vertice] = curvatures.gaussian_curv;
		m_mean_curvature[*curr_vertice] = curvatures.mean_curv;
	}

	m_vertex_curv_type.SetGraph(vertices_graph);
	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		m_vertex_curv_type[*curr_vertice] = GetCurvatureType(m_gaussian_curvature[*curr_vertice], m_mean_curvature[*curr_vertice]);
	}
}

void SngPtsFinderScaleSpace::CalcScaleSpacePropsHessianRatio(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
		get(boost::vertex_info_3d, vertices_graph);
	CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());
	m_props_hessian_ratio.resize(m_scale_space_levels_num);
	m_grad_dx.resize(m_scale_space_levels_num);
	m_grad_dy.resize(m_scale_space_levels_num);

	m_projected_grad_x.resize(m_scale_space_levels_num);
	m_projected_grad_y.resize(m_scale_space_levels_num);;
	m_projected_grad_norm.resize(m_scale_space_levels_num);
	m_props_hessian_ratio_of_proj.resize(m_scale_space_levels_num);
	m_hessian_map_calculators.resize(m_scale_space_levels_num);
	//calculate gradients and hessians fro each level of each property
#pragma omp parallel for
	for (int curr_level = 0; curr_level < m_scale_space_levels_num; ++curr_level)
	{
		m_props_hessian_ratio[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_grad_dx[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);
		m_grad_dy[curr_level].resize(ISingularPointsFinder::SURF_PROPS_NUM);

		for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++curr_prop)
		{
			m_props_hessian_ratio[curr_level][curr_prop].SetGraph(vertices_graph);
			m_hessian_map_calculators[curr_level].Process(vertices_graph, coord_map, m_output_scale_space[curr_level][curr_prop],
				m_vert_vert_dist, 1.0, m_tangent_basis_map, true, m_props_hessian_ratio[curr_level][curr_prop]);
			m_grad_dx[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_dx;
			m_grad_dy[curr_level][curr_prop] = m_hessian_map_calculators[curr_level].m_dy;
		}

		m_projected_grad_x[curr_level].SetGraph(vertices_graph);
		m_projected_grad_y[curr_level].SetGraph(vertices_graph);
		m_projected_grad_norm[curr_level].SetGraph(vertices_graph);
		m_props_hessian_ratio_of_proj[curr_level].SetGraph(vertices_graph);
		//find projections
		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			m_projected_grad_x[curr_level][*curr_vertice] = 
				ProjectPCADiff(m_grad_dx[curr_level], coord_map, m_scale_space_projecter[curr_level], *curr_vertice);
			m_projected_grad_y[curr_level][*curr_vertice] = 
				ProjectPCADiff(m_grad_dy[curr_level], coord_map, m_scale_space_projecter[curr_level], *curr_vertice);
			m_projected_grad_norm[curr_level][*curr_vertice] = 
				sqrt(Sqr(m_projected_grad_x[curr_level][*curr_vertice]) + Sqr(m_projected_grad_y[curr_level][*curr_vertice]));
		}
		//find hessian of projections
		m_hessian_map_calculators[curr_level].ProcessFromGradient(vertices_graph, coord_map, 
			m_projected_grad_x[curr_level],m_projected_grad_y[curr_level],
			m_vert_vert_dist, 1.0, m_tangent_basis_map, true, m_props_hessian_ratio_of_proj[curr_level]);
		
	}
}

void SngPtsFinderScaleSpace::CalcHessianOfProjectedLog(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
		get(boost::vertex_info_3d, vertices_graph);
	CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());
	//find hessian of LOG projections
	m_props_hessian_ratio_of_proj_LOG.resize(m_scale_space_levels_num - 1);

	for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
	{
		m_props_hessian_ratio_of_proj_LOG[curr_level].SetGraph(vertices_graph);
		m_hessian_map_calculators[0].Process(vertices_graph, coord_map, m_detector_function_projected[curr_level],
			m_vert_vert_dist, 1.0, m_tangent_basis_map, true, m_props_hessian_ratio_of_proj_LOG[curr_level]);
	}
}

void SngPtsFinderScaleSpace::CalculateDistanceMaps(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	//mean edge length
	double mean_dist = 0.0;
	double edges_num = 0.0;
	std::vector<double> edges_lengths;
	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
			get(boost::vertex_info_3d, vertices_graph);
		CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());

		for (auto curr_neighb = adjacent_vertices(*curr_vertice, vertices_graph).first, 
			end_neighb = adjacent_vertices(*curr_vertice, vertices_graph).second; 
			curr_neighb != end_neighb; ++curr_neighb)
		{
			const double curr_dist = cv::norm(coord_map[*curr_vertice] - coord_map[*curr_neighb]);
			edges_lengths.push_back(curr_dist);
			mean_dist += curr_dist;
			++edges_num;
		}
	}	
	std::nth_element(edges_lengths.begin(), edges_lengths.begin() + Round(0.9 * edges_lengths.size()), edges_lengths.end());
	std::cout << mean_dist / edges_num << " " << edges_lengths[Round(0.9 * edges_lengths.size())] << "\n";

	//distance maps

	const TrianglesGraph& triangles_graph = mesh.triangles;
	boost::property_map<const TrianglesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map_tr = 
		get(boost::vertex_info_3d, triangles_graph);
	const auto coord_map_tr = GetProxyPropMap(coord_3d_map_tr, GetCoord<MeshTriangle>());

	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
		get(boost::vertex_info_3d, vertices_graph);
	CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());
	//calculate distances
	m_vert_vert_dist.create(num_vertices(vertices_graph), num_vertices(vertices_graph));
	m_vert_tr_dist.create(num_vertices(vertices_graph), num_vertices(triangles_graph));

	if (m_use_euclid_distance)
	{
		CalcDistBetweenGraphs(vertices_graph, coord_map, vertices_graph, coord_map, m_vert_vert_dist);
		CalcDistBetweenGraphs(vertices_graph, coord_map, triangles_graph, coord_map_tr, m_vert_tr_dist);
	}
	else
	{//USE GEODESIC DISTANCE
		DijkstraDistMapCalculator<VerticesGraph, double> djikstra_dist;
		djikstra_dist.Calc(vertices_graph, m_vert_vert_dist);

		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			for (auto curr_tr = vertices(triangles_graph).first, end_tr = vertices(triangles_graph).second; 
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

void SngPtsFinderScaleSpace::CalcSingPtsFromCurvatureScales(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	const TrianglesGraph& triangles_graph = mesh.triangles;

	//DATA PREPARATION
	m_input_prop_map.resize(ISingularPointsFinder::SURF_PROPS_NUM);
	m_input_prop_map[ISingularPointsFinder::GAUSS_CURV] = m_gaussian_curvature;
	m_input_prop_map[ISingularPointsFinder::MEAN_CURV] = m_mean_curvature;
	m_input_prop_map[ISingularPointsFinder::ELECTR_POTENT] = m_vertex_charge_map;
	m_input_prop_map[ISingularPointsFinder::STERIC_POTENT] = m_vertex_lennard_jones_map;
	
	if (!m_mean_and_sigma.empty())
	{//rescale functions
		for (int curr_prop = ISingularPointsFinder::FIRST_SURF_PROP; curr_prop < ISingularPointsFinder::SURF_PROPS_NUM; ++ curr_prop)
		{
			const double min_val = m_mean_and_sigma[curr_prop][0] - 3.0 * m_mean_and_sigma[curr_prop][1];
			const double max_val = m_mean_and_sigma[curr_prop][0] + 3.0 * m_mean_and_sigma[curr_prop][1];

			for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
				curr_vertice != end_vertices; ++curr_vertice)
			{
				m_input_prop_map[curr_prop][*curr_vertice] = std::min(m_input_prop_map[curr_prop][*curr_vertice], max_val);
				m_input_prop_map[curr_prop][*curr_vertice] = std::max(m_input_prop_map[curr_prop][*curr_vertice], min_val);
				m_input_prop_map[curr_prop][*curr_vertice] -= min_val;
				m_input_prop_map[curr_prop][*curr_vertice] /= m_mean_and_sigma[curr_prop][1];
			}
		}
	}
	//PREPAIR FOR MAKING SCALE SPACE
	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
		get(boost::vertex_info_3d, vertices_graph);
	CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());

	
	//MAKE SCALE-SPACE
	const bool use_post_filter = !m_detect_blobs;
	std::cout << "scale-space" << std::endl;
	m_scale_space_blurrer.MakeScaleSpace(vertices_graph, m_vert_vert_dist, triangles_graph, m_vert_tr_dist, m_input_prop_map, 
		m_scale_space_levels_num, use_post_filter, m_output_scale_space);
	//DETECT TYPES OF VERTICES FOR EACH LEVEL, ACCORDING TO THE CURVATURE	
	m_vertex_curv_type_mesh_levels.assign(m_scale_space_levels_num, VetrticesCurvMap());
	for (size_t curr_level = 0; curr_level < m_scale_space_levels_num; ++curr_level)
	{
		m_vertex_curv_type_mesh_levels[curr_level].SetGraph(vertices_graph);
		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
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
	m_scale_space_projecter.resize(m_scale_space_levels_num);

	if (m_detect_blobs)
	{
		m_output_scale_space_diff.clear();
		m_output_scale_space_diff.resize(m_scale_space_levels_num - 1);

		std::cout << "LOG calculation" << std::endl;

		if (m_use_DOG_as_LOG_approximation)
		{
			//ScaleSpaceBlurrer<VerticesGraph, CoordMap, GaussianKernel<cv::Point3d, double>> scale_space_blurrer_shifted;//for calculation of derivative in t
			//const double d_mult_sigma = pow(m_sigma_max / m_init_curv_sigma, 1.0 / (m_sing_pts_levels_num));
			//scale_space_blurrer_shifted.Init(m_init_curv_sigma, d_mult_sigma, 0.1, false);
			//std::vector<std::vector<DoubleVertGraphProp>> output_scale_space_shifted;
			//scale_space_blurrer_shifted.MakeScaleSpace(vertices_graph, m_vert_vert_dist, triangles_graph, m_vert_tr_dist, m_input_prop_map, 
			//	m_scale_space_levels_num, use_post_filter, output_scale_space_shifted);
			//Calculate scales difference which approximates laplace operator
			for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
			{
				/*const double coeff = (scale_space_blurrer_shifted.GetSigma(curr_level) + m_scale_space_blurrer.GetSigma(curr_level)) / 2.0 
					/ (scale_space_blurrer_shifted.GetSigma(curr_level) - m_scale_space_blurrer.GetSigma(curr_level));*/

				const double coeff = (m_scale_space_blurrer.GetSigma(curr_level + 1) + m_scale_space_blurrer.GetSigma(curr_level)) / 2.0 
					/ (m_scale_space_blurrer.GetSigma(curr_level + 1) - m_scale_space_blurrer.GetSigma(curr_level));

				m_output_scale_space_diff[curr_level].resize(m_output_scale_space[curr_level].size());
				const auto& shifted_scale_space_lev = m_output_scale_space[curr_level + 1]/*output_scale_space_shifted[curr_level]*/;
				for (size_t curr_prop = 0; curr_prop < m_output_scale_space_diff[curr_level].size(); ++curr_prop)
				{
					m_output_scale_space_diff[curr_level][curr_prop].SetGraph(vertices_graph);
					for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
						curr_vertice != end_vertices; ++curr_vertice)
					{
						m_output_scale_space_diff[curr_level][curr_prop][*curr_vertice] = coeff * 				 
							(shifted_scale_space_lev[curr_prop][*curr_vertice] - m_output_scale_space[curr_level][curr_prop][*curr_vertice]);

					}	
				}
			}
		}
		else
		{
			for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
			{
				const double coeff = m_scale_space_blurrer.GetSigma(curr_level);
				m_output_scale_space_diff[curr_level].resize(m_output_scale_space[curr_level].size());

				for (size_t curr_prop = 0; curr_prop < m_output_scale_space_diff[curr_level].size(); ++curr_prop)
				{
					m_output_scale_space_diff[curr_level][curr_prop].SetGraph(vertices_graph);
					/*SimpleLaplacian(vertices_graph, m_output_scale_space[curr_level][curr_prop], 
						output_scale_space_diff[curr_level][curr_prop]);*/
				}

				//LaplaceBeltramiKernelWeightedDist<double, SignedDistFunc<double>> kernel(0.5/*(m_scale_space_blurrer.GetSigma(curr_level) + 1.7) / 4.0*//*0.3*/, SignedDistFunc<double>());
				LoGKernelWeightedDist<double> kernel(m_scale_space_blurrer.GetSigma(curr_level));
				FilterMeshWeightedFunc(vertices_graph, triangles_graph, 
					m_vert_vert_dist, m_vert_tr_dist, m_input_prop_map, true,
					kernel, m_output_scale_space_diff[curr_level]);
				/*FilterMeshWeightedFunc(vertices_graph, triangles_graph, 
				vert_vert_dist, vert_tr_dist, m_output_scale_space[curr_level], false,
				kernel, output_scale_space_diff[curr_level]);*/
					
				for (size_t curr_prop = 0; curr_prop < m_output_scale_space_diff[curr_level].size(); ++curr_prop)
				{
					for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
						curr_vertice != end_vertices; ++curr_vertice)
					{
						m_output_scale_space_diff[curr_level][curr_prop][*curr_vertice] = m_output_scale_space_diff[curr_level][curr_prop][*curr_vertice] * coeff;
					}	
				}
			}

			for (size_t curr_prop = 0; curr_prop < m_output_scale_space_diff[0].size(); ++curr_prop)
			{
				for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
					curr_vertice != end_vertices; ++curr_vertice)
				{
					for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
					{
						//std::cout << output_scale_space_diff[curr_level][curr_prop][*curr_vertice] << " ";
					}
					//std::cout << "\n";
				}	
			}
		}
		//END LOG CALCULATION

		//recalculate scales values as it is between scale
		/*for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
		{
			for (size_t curr_prop = 0; curr_prop < output_scale_space_diff[curr_level].size(); ++curr_prop)
			{
				for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
					curr_vertice != end_vertices; ++curr_vertice)
				{
					m_output_scale_space[curr_level][curr_prop][*curr_vertice] = 0.5 *				 
						(m_output_scale_space[curr_level + 1][curr_prop][*curr_vertice] + 
						m_output_scale_space[curr_level][curr_prop][*curr_vertice]);

				}	
			}
		}*/
		//CALCULATE MEAN IN NEIGHB
		std::cout << "find PCA" << std::endl;
		m_scale_space_projecter.resize(m_output_scale_space.size());
		for (int lev = 0; lev < m_output_scale_space.size(); ++lev)
		{
			m_scale_space_projecter[lev].SetGraph(vertices_graph);
			FindVectorsForProjection(vertices_graph, coord_map, m_output_scale_space[lev], 
				m_vert_vert_dist, 2.0 * m_scale_space_blurrer.GetSigma(lev), m_scale_space_projecter[lev]);
		}
		//CALCULATE PROJECTED VECTORS
		m_detector_function_projected.clear();
		m_detector_function_projected.resize(m_output_scale_space_diff.size());
		m_scale_space_projected.clear();
		m_scale_space_projected.resize(m_output_scale_space_diff.size());
		m_projecters_coords.clear();
		m_projecters_coords.resize(m_output_scale_space_diff.size());

		for (int lev = 0; lev < m_output_scale_space_diff.size(); ++lev)
		{
			m_detector_function_projected[lev].SetGraph(vertices_graph);
			m_scale_space_projected[lev].SetGraph(vertices_graph);
			const size_t kPcaPropsNum = ISingularPointsFinder::LAST_PCA_PROP - ISingularPointsFinder::FIRST_PCA_PROP + 1;
			m_projecters_coords[lev].resize(kPcaPropsNum);
			for (size_t prop_num = 0; prop_num < kPcaPropsNum; ++prop_num)
			{
				m_projecters_coords[lev][prop_num].SetGraph(vertices_graph);
			}
			for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
				curr_vertice != end_vertices; ++curr_vertice)
			{
				m_detector_function_projected[lev][*curr_vertice] = 
					ProjectPCADiff(m_output_scale_space_diff[lev], coord_map, m_scale_space_projecter[lev], *curr_vertice);
				m_scale_space_projected[lev][*curr_vertice] = 
					ProjectPCA(m_output_scale_space[lev], coord_map, m_scale_space_projecter[lev], *curr_vertice);

				for (size_t prop_num = 0; prop_num < kPcaPropsNum; ++prop_num)
				{
					m_projecters_coords[lev][prop_num][*curr_vertice]
						= std::abs(m_scale_space_projecter[lev][*curr_vertice].pca.eigenvectors.at<double>(3, 3 + prop_num));
				}
			}			
		}

		std::cout << "hessian" << std::endl;
		CalcScaleSpacePropsHessianRatio(mesh);
		CalcHessianOfProjectedLog(mesh);
		//////////////////////////////////////////////////////////////////////////
		//DETECT POINTS
		m_maximums_with_levels.clear();

		if (m_combine_channels)
		{
			std::vector<DoubleVertGraphProp> gaussian_curvature_scales_diff(m_scale_space_levels_num - 1, DoubleVertGraphProp(vertices_graph));
			std::vector<std::vector<VertexDescriptor>> maximums_gauss(m_scale_space_levels_num);
			m_maximums_with_levels.resize(m_scale_space_levels_num);
			for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
			{
				gaussian_curvature_scales_diff[curr_level] = m_output_scale_space_diff[curr_level][STERIC_POTENT];
			
				FindLocalMaxAndMin(vertices_graph, m_output_scale_space[curr_level][STERIC_POTENT], maximums_gauss[curr_level], 
						std::greater<double>(), std::less<double>());	
				/*FindLocalMaximumsOfAbsVal(vertices_graph, gaussian_curvature_scales_diff[curr_level], m_maximums_with_levels[curr_level], 
					std::greater<double>(), std::less<double>());*/	
			}

			m_vertex_curv_type_mesh_levels.resize(m_scale_space_levels_num - 1);
		/*	FindLocalMaximumsOnLevels(vertices_graph, gaussian_curvature_scales_diff, m_maximums_with_levels, 
				std::greater<double>(), std::less<double>(), m_one_ring_neighb);*/
			m_maximums_with_levels.clear();
			FindLocalMaximumsOnLevelsVect(vertices_graph,coord_map, m_output_scale_space_diff, m_scale_space_projecter,  
				m_vert_vert_dist, 1.0, std::greater<double>(), std::less<double>(), false, m_maximums_with_levels);

			for (size_t curr_level = 0; curr_level < m_maximums_with_levels.size(); ++curr_level)
			{
				std::cout << m_maximums_with_levels[curr_level].size() << " ";
			}
			std::cout << "\n";
			for (size_t curr_level = 0; curr_level < m_maximums_with_levels.size(); ++curr_level)
			{
				std::vector<VertexDescriptor> not_near_ridge;

				for (int ind1 = 0; ind1 < m_maximums_with_levels[curr_level].size();++ind1)
				{
					bool is_near = false;

					std::vector<VertexDescriptor> nearest_vertices;
					GetVerticesWithinDistPlusAdjacent(m_maximums_with_levels[curr_level][ind1], vertices_graph, m_vert_vert_dist, 0.4, nearest_vertices);

					for (auto neighb_it = nearest_vertices.begin(),	end_neighb = nearest_vertices.end(); neighb_it != end_neighb; ++neighb_it)					
					{
						if (m_props_hessian_ratio_of_proj[curr_level + m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls][*neighb_it] >= 5.0)
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
			/*for (size_t curr_level = 0; curr_level < m_maximums_with_levels.size(); ++curr_level)
			{
				std::vector<VertexDescriptor> with_max;
				std::cout << maximums_gauss[curr_level + 1].size() << " ";

				for (int ind1 = 0; ind1 < m_maximums_with_levels[curr_level].size();++ind1)
				{
					bool is_near = false;
					VertexDescriptor matched_point;
					if (std::find(maximums_gauss[curr_level + 1].begin(), maximums_gauss[curr_level + 1].end(), m_maximums_with_levels[curr_level][ind1]) 
						!= maximums_gauss[curr_level + 1].end())
					{
						is_near = true;
						matched_point = m_maximums_with_levels[curr_level][ind1];
					}	
					else
					{
						for (auto neighb_it = adjacent_vertices(m_maximums_with_levels[curr_level][ind1], vertices_graph).first,
							end_neighb = adjacent_vertices(m_maximums_with_levels[curr_level][ind1], vertices_graph).second; 
							neighb_it != end_neighb; ++neighb_it)
						{
							if (std::find(maximums_gauss[curr_level + 1].begin(), maximums_gauss[curr_level + 1].end(), *neighb_it) 
								!= maximums_gauss[curr_level + 1].end())
							{
								is_near = true;
								matched_point = *neighb_it;
							}	
						}
					}

					if (is_near)
					{
						with_max.push_back(matched_point);
					}
				}*/

				/*std::sort(with_max.begin(), with_max.end());
				with_max.erase(std::unique(with_max.begin(), with_max.end()), with_max.end());*/
				//m_maximums_with_levels[curr_level] = with_max;
				/*m_maximums_with_levels[curr_level].insert(m_maximums_with_levels[curr_level].end(), 
					maximums_gauss[curr_level].begin(), maximums_gauss[curr_level].end());*/
			/*}*/
			std::cout << "\n";

			m_maximums_with_levels.resize(m_sing_pts_levels_num);
			for (size_t curr_level = 1; curr_level < m_maximums_with_levels.size(); ++curr_level)
			{
				m_maximums_with_levels[0].insert(m_maximums_with_levels[0].end(), 
					m_maximums_with_levels[curr_level].begin(), m_maximums_with_levels[curr_level].end());
			}
			
		}
		else
		{
			std::vector<DoubleVertGraphProp> gaussian_curvature_scales_diff(m_scale_space_levels_num - 1, DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> mean_curvature_scales_diff(m_scale_space_levels_num - 1 , DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> electrostatic_scale_space_diff(m_scale_space_levels_num - 1, DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> steric_scale_space_diff(m_scale_space_levels_num - 1 , DoubleVertGraphProp(vertices_graph));

			for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
			{
				gaussian_curvature_scales_diff[curr_level] = m_output_scale_space_diff[curr_level][GAUSS_CURV];
				mean_curvature_scales_diff[curr_level] = m_output_scale_space_diff[curr_level][MEAN_CURV];
				electrostatic_scale_space_diff[curr_level] = m_output_scale_space_diff[curr_level][ELECTR_POTENT];
				steric_scale_space_diff[curr_level] = m_output_scale_space_diff[curr_level][STERIC_POTENT];
			}
			m_vertex_curv_type_mesh_levels.resize(m_scale_space_levels_num - 1);
			std::vector<std::vector<VertexDescriptor>> sing_pts_gaussian, sing_pts_mean, sing_pts_electric, sing_pts_steric;
			FindLocalMaximumsOnLevels(vertices_graph, gaussian_curvature_scales_diff, sing_pts_gaussian, std::greater<double>(), std::less<double>(), m_one_ring_neighb);
			FindLocalMaximumsOnLevels(vertices_graph, mean_curvature_scales_diff, sing_pts_mean, std::greater<double>(), std::less<double>(), m_one_ring_neighb);
			FindLocalMaximumsOnLevels(vertices_graph, electrostatic_scale_space_diff, sing_pts_electric, std::greater<double>(), std::less<double>(), m_one_ring_neighb);
			FindLocalMaximumsOnLevels(vertices_graph, steric_scale_space_diff, sing_pts_steric, std::greater<double>(), std::less<double>(), m_one_ring_neighb);

			for (size_t curr_level = 0; curr_level < m_sing_pts_levels_num; ++curr_level)
			{
				PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, vertices_graph);
				points_keeper_curr_scale.AddPoints(sing_pts_gaussian[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
				points_keeper_curr_scale.AddPoints(sing_pts_mean[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
				std::cout <<points_keeper_curr_scale.Points().size() << " ";
				points_keeper_curr_scale.AddPoints(sing_pts_electric[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
				points_keeper_curr_scale.AddPoints(sing_pts_steric[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
				m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
				std::cout <<points_keeper_curr_scale.Points().size() << "\n";
			}
		}
	}
	else
	{
		m_maximums_with_levels.clear();

		m_scale_space_projecter.resize(m_output_scale_space.size());
		for (int lev = 0; lev < m_output_scale_space.size(); ++lev)
		{
			m_scale_space_projecter[lev].SetGraph(vertices_graph);
			FindVectorsForProjection(vertices_graph, coord_map, m_output_scale_space[lev], 
				m_vert_vert_dist, 2.0 * m_scale_space_blurrer.GetSigma(lev), m_scale_space_projecter[lev]);
		}

		if (m_combine_channels)
		{
			std::vector<DoubleVertGraphProp> scale_space_projected(m_scale_space_levels_num);

			for (size_t curr_level = 0; curr_level < m_scale_space_levels_num; ++curr_level)
			{
				scale_space_projected[curr_level].SetGraph(vertices_graph);

				for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
					curr_vertice != end_vertices; ++curr_vertice)
				{
					scale_space_projected[curr_level][*curr_vertice] = 0.0;
					double norm = 0;
					for (int prop = 0; prop < SURF_PROPS_NUM; ++prop)
					{
						scale_space_projected[curr_level][*curr_vertice] += 
							m_output_scale_space[curr_level][prop][*curr_vertice];
							///*scale_space_vect_to_project*/m_output_scale_space[curr_level][prop][*curr_vertice];
						/*norm += scale_space_vect_to_project[curr_level][prop][*curr_vertice] * 
							scale_space_vect_to_project[curr_level][prop][*curr_vertice];*/
					}
					//scale_space_projected[curr_level][*curr_vertice] /= sqrt(norm);
				}
			}

			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, scale_space_projected,
				m_vertex_curv_type_mesh_levels, 1.0, m_maximums_with_levels, true);
			m_maximums_with_levels.resize(m_sing_pts_levels_num);
		}
		else
		{
			std::vector<DoubleVertGraphProp> gaussian_curvature_scales(m_scale_space_levels_num, DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> mean_curvature_scales(m_scale_space_levels_num , DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> electrostatic_scale_space(m_scale_space_levels_num, DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> steric_scale_space(m_scale_space_levels_num , DoubleVertGraphProp(vertices_graph));

			for (size_t curr_level = 0; curr_level < m_scale_space_levels_num; ++curr_level)
			{
				gaussian_curvature_scales[curr_level] = m_output_scale_space[curr_level][GAUSS_CURV];
				mean_curvature_scales[curr_level] = m_output_scale_space[curr_level][MEAN_CURV];
				electrostatic_scale_space[curr_level] = m_output_scale_space[curr_level][ELECTR_POTENT];
				steric_scale_space[curr_level] = m_output_scale_space[curr_level][STERIC_POTENT];
			}

			std::vector<std::vector<VertexDescriptor>> sing_pts_gaussian, sing_pts_mean, sing_pts_electric, sing_pts_steric;
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, gaussian_curvature_scales,
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_gaussian, true);
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, mean_curvature_scales,
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_mean, true);
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, electrostatic_scale_space,
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_electric, true);
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, steric_scale_space,
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_steric, true);

			for (size_t curr_level = 0; curr_level < m_sing_pts_levels_num; ++curr_level)
			{
				PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, vertices_graph);
				points_keeper_curr_scale.AddPoints(sing_pts_gaussian[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
				points_keeper_curr_scale.AddPoints(sing_pts_mean[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
				std::cout <<points_keeper_curr_scale.Points().size() << " ";
				points_keeper_curr_scale.AddPoints(sing_pts_electric[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
				points_keeper_curr_scale.AddPoints(sing_pts_steric[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
				m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
				std::cout <<points_keeper_curr_scale.Points().size() << "\n";
			}
		}
	}
	///////
	PointsKeeper<VerticesGraph> points_keeper(1.0, vertices_graph);

	for (size_t curr_level = 0; curr_level < m_sing_pts_levels_num; ++curr_level)
	{
		std::cout << m_maximums_with_levels[curr_level].size() << " ";	

		points_keeper.AddPoints(m_maximums_with_levels[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
			+ m_diff_btwn_sng_pts_lvls_and_scl_spc_lvls]);
	}
	m_maximums = points_keeper.Points();

	std::cout << "maximums " << m_maximums.size() << " ";	
}

//void SngPtsFinderScaleSpace::CalcSingPtsFromCurvatureScales(const Mesh& mesh, const int levels_num)
//{
//	const VerticesGraph& vertices_graph = mesh.vertices;
//	typedef ProxyPropMap<
//		boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type, GetCoord<Vertice>> CoordMap;
//	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
//		get(boost::vertex_info_3d, vertices_graph);
//	CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());
//	ScaleSpaceBlurrer<VerticesGraph, CoordMap, GaussianKernel<cv::Point3d, double>> scale_space_blurrer;
//	const size_t ext_levels_num = levels_num + 2;
//	const bool kAdditive = false;
//	const double d_mult_sigma = pow(kSigmaMax / kInitCurvSigma, 1.0 / ext_levels_num);
//	std::cout << "d_mult_sigma " << d_mult_sigma << std::endl;
//	scale_space_blurrer.Init(kInitCurvSigma, d_mult_sigma, ext_levels_num, kAdditive);
//	typedef ProxyPropMap<DoubleVertGraphProp, IdenticalTransformFunc> DouplePropMapRef;
//	std::vector<DoubleVertGraphProp> input_prop_map(ISingularPointsFinder::PROPS_NUM);
//	input_prop_map[ISingularPointsFinder::GAUSS_CURV] = m_gaussian_curvature;
//	input_prop_map[ISingularPointsFinder::MEAN_CURV] = m_mean_curvature;
//	input_prop_map[ISingularPointsFinder::ELECTR_POTENT] = m_vertex_charge_map;
//	input_prop_map[ISingularPointsFinder::STERIC_POTENT] = m_vertex_lennard_jones_map;
//
//	if (!m_mean_and_sigma.empty())
//	{//rescale functions
//		for (int curr_prop = ISingularPointsFinder::FIRST_PROP; curr_prop < ISingularPointsFinder::PROPS_NUM; ++ curr_prop)
//		{
//			for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
//				curr_vertice != end_vertices; ++curr_vertice)
//			{
//				input_prop_map[curr_prop][*curr_vertice] = std::min(input_prop_map[curr_prop][*curr_vertice], m_mean_and_sigma[curr_prop][0] + 3.0 * m_mean_and_sigma[curr_prop][1]);
//				input_prop_map[curr_prop][*curr_vertice] = std::max(input_prop_map[curr_prop][*curr_vertice], m_mean_and_sigma[curr_prop][0] - 3.0 * m_mean_and_sigma[curr_prop][1]);
//				input_prop_map[curr_prop][*curr_vertice] /= m_mean_and_sigma[curr_prop][1];
//			}
//		}
//	}
//	std::vector<std::vector<DoubleVertGraphProp>> output_scale_space;
//	scale_space_blurrer.MakeScaleSpace(vertices_graph, coord_map, input_prop_map, ext_levels_num + 1, false/*true*/, output_scale_space);
//	//попробуем определить радиус фич
//	std::vector<std::vector<DoubleVertGraphProp>> output_scale_space_diff(ext_levels_num);
//	std::vector<std::vector<DoubleVertGraphProp>> output_scale_space_mean(ext_levels_num);
//	DoubleVertGraphProp temp_prop_map(vertices_graph);
//
//	for (size_t curr_level = 0; curr_level < ext_levels_num; ++curr_level)
//	{
//		const double coeff = /*log*/(/*1.0 + */(scale_space_blurrer.GetSigma(curr_level + 1) + scale_space_blurrer.GetSigma(curr_level)) / 2.0) 
//			/** log(1.0 + (scale_space_blurrer.GetSigma(curr_level + 1) + scale_space_blurrer.GetSigma(curr_level)) / 2.0) */
//			/ (scale_space_blurrer.GetSigma(curr_level + 1) - scale_space_blurrer.GetSigma(curr_level));
//
//		output_scale_space_diff[curr_level].resize(output_scale_space[curr_level].size());
//		output_scale_space_mean[curr_level].resize(output_scale_space[curr_level].size());
//		for (size_t curr_prop = 0; curr_prop < output_scale_space_diff[curr_level].size(); ++curr_prop)
//		{
//			output_scale_space_diff[curr_level][curr_prop].SetGraph(vertices_graph);
//			output_scale_space_mean[curr_level][curr_prop].SetGraph(vertices_graph);
//			for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
//				curr_vertice != end_vertices; ++curr_vertice)
//			{
//				temp_prop_map[*curr_vertice] = coeff * 				 
//					(output_scale_space[curr_level + 1][curr_prop][*curr_vertice] - output_scale_space[curr_level][curr_prop][*curr_vertice]);
//				output_scale_space_diff[curr_level][curr_prop][*curr_vertice] = coeff * 				 
//					(output_scale_space[curr_level + 1][curr_prop][*curr_vertice] - output_scale_space[curr_level][curr_prop][*curr_vertice]);
//
//			}	
//			typedef MedianKernel<double, size_t> MedianFilter;
//			typedef AverageKernel<double, size_t> AverageFilter;
//			const size_t kMedianRadius = 1;
//			MedianFilter med_filter(kMedianRadius);
//			AverageFilter av_filter(1);
//			//FilterGraphEdgeDist(med_filter, vertices_graph, temp_prop_map, output_scale_space_diff[curr_level][curr_prop]);
//			FilterGraphEdgeDist(av_filter, vertices_graph, output_scale_space[curr_level][curr_prop], output_scale_space_mean[curr_level][curr_prop]);
//		}
//	}
//	///////
//	//DETECT TYPES OF VERTICES FOR EACH LEVEL, ACCORDING TO THE CURVATURE
//	m_gaussian_curvature_scales.resize(ext_levels_num);
//	m_mean_curvature_scales.resize(ext_levels_num);
//	std::vector<DoubleVertGraphProp> gaussian_curvature_scales_diff(ext_levels_num, DoubleVertGraphProp(vertices_graph));
//	std::vector<DoubleVertGraphProp> mean_curvature_scales_diff(ext_levels_num , DoubleVertGraphProp(vertices_graph));
//	std::vector<DoubleVertGraphProp> electrostatic_scale_space_diff(ext_levels_num, DoubleVertGraphProp(vertices_graph));
//	std::vector<DoubleVertGraphProp> steric_scale_space_diff(ext_levels_num , DoubleVertGraphProp(vertices_graph));
//
//	for (size_t curr_level = 0; curr_level < ext_levels_num; ++curr_level)
//	{
//		m_gaussian_curvature_scales[curr_level] = output_scale_space[curr_level][0];
//		m_mean_curvature_scales[curr_level] = output_scale_space[curr_level][1];
//
//		gaussian_curvature_scales_diff[curr_level] = output_scale_space_diff[curr_level][0];
//		mean_curvature_scales_diff[curr_level] = output_scale_space_diff[curr_level][1];
//		electrostatic_scale_space_diff[curr_level] = output_scale_space_diff[curr_level][2];
//		steric_scale_space_diff[curr_level] = output_scale_space_diff[curr_level][3];
//	}
//
//	m_vertex_curv_type_mesh_levels.assign(ext_levels_num, VetrticesCurvMap());
//
//	for (size_t curr_level = 0; curr_level < ext_levels_num; ++curr_level)
//	{
//		m_vertex_curv_type_mesh_levels[curr_level].SetGraph(vertices_graph);
//		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
//			curr_vertice != end_vertices; ++curr_vertice)
//		{
//			m_vertex_curv_type_mesh_levels[curr_level][*curr_vertice] = GetCurvatureType(
//				m_gaussian_curvature_scales[curr_level][*curr_vertice], m_mean_curvature_scales[curr_level][*curr_vertice]);
//		}	
//	}	
//	const bool kDetectSeparately = false;
//	PointsKeeper<VerticesGraph> points_keeper(1.0, vertices_graph);
//	m_maximums_with_levels.clear();
//
//	//if (kDetectSeparately)
//	{
//		/*std::vector<std::vector<VertexDescriptor>> sing_pts_gaussian, sing_pts_mean;
//		m_gaussian_curvature_scales.resize(ext_levels_num);
//		m_mean_curvature_scales.resize(ext_levels_num);
//
//		FindScaleSingularPointsOnFunc(vertices_graph, coord_map, m_gaussian_curvature_scales, 
//		m_vertex_curv_type_mesh_levels, 1.0, sing_pts_gaussian, false, true);
//		FindScaleSingularPointsOnFunc(vertices_graph, coord_map, m_mean_curvature_scales, 
//		m_vertex_curv_type_mesh_levels, 1.0, sing_pts_mean, false, true);*/
//		//m_maximums_with_levels.clear();
//
//		//for (size_t curr_level = 0; curr_level < levels_num; ++curr_level)
//		//{
//		//	PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, vertices_graph);
//		//	points_keeper_curr_scale.AddPoints(sing_pts_gaussian[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//		//	points_keeper_curr_scale.AddPoints(sing_pts_mean[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//		//	m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
//		//	points_keeper.AddPoints(points_keeper_curr_scale.Points(), m_vertex_curv_type_mesh_levels[0]);
//		//	std::cout << m_maximums_with_levels[curr_level].size() << "\n";
//		//}
//		//m_maximums = points_keeper.Points();
//		//повторим для разностей
//		std::vector<std::vector<VertexDescriptor>> sing_pts_gaussian1, sing_pts_mean1, sing_pts_electric, sing_pts_steric;
//		FindScaleSingularPointsOnFunc(vertices_graph, coord_map, gaussian_curvature_scales_diff, m_vertex_curv_type_mesh_levels, 1.0, sing_pts_gaussian1, true, false);
//		FindScaleSingularPointsOnFunc(vertices_graph, coord_map, mean_curvature_scales_diff, m_vertex_curv_type_mesh_levels, 1.0, sing_pts_mean1, true, false);
//
//		FindScaleSingularPointsOnFunc(vertices_graph, coord_map, electrostatic_scale_space_diff, m_vertex_curv_type_mesh_levels, 1.0, sing_pts_electric, true, false);
//		FindScaleSingularPointsOnFunc(vertices_graph, coord_map, steric_scale_space_diff, m_vertex_curv_type_mesh_levels, 1.0, sing_pts_steric, true, false);
//		for (size_t curr_level = 0; curr_level < levels_num; ++curr_level)
//		{
//			PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, vertices_graph);
//			points_keeper_curr_scale.AddPoints(sing_pts_gaussian1[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//			points_keeper_curr_scale.AddPoints(sing_pts_mean1[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//			std::cout <<points_keeper_curr_scale.Points().size() << " ";
//			points_keeper_curr_scale.AddPoints(sing_pts_electric[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//			points_keeper_curr_scale.AddPoints(sing_pts_steric[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//			m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
//			points_keeper.AddPoints(points_keeper_curr_scale.Points(), m_vertex_curv_type_mesh_levels[0]);
//			std::cout <<points_keeper_curr_scale.Points().size() << "\n";
//		}
//	}
//	//else
//	{
//		points_keeper.Clear();
//		m_maximums_with_levels.clear();
//		FindLocalMaximumsOnLevelsVect(vertices_graph, output_scale_space_diff, output_scale_space_mean/*output_scale_space*/,  
//			m_maximums_with_levels, std::greater<double>(), std::less<double>());
//		for (size_t curr_level = 0; curr_level < levels_num; ++curr_level)
//		{
//			
//			std::cout <<m_maximums_with_levels[curr_level].size() << " ";
//			points_keeper.AddPoints(m_maximums_with_levels[curr_level], m_vertex_curv_type_mesh_levels[0]);
//		}
//	}
//	m_maximums = points_keeper.Points();
//
//	std::cout << "maximums " << m_maximums.size() << " ";	
//}

//void SngPtsFinderScaleSpace::CalcSingPtsFromCurvatureScalesOptimalVector(const Mesh& mesh, const int levels_num)
//{
//	const VerticesGraph& vertices_graph = mesh.vertices;
//	typedef ProxyPropMap<
//		boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type, GetCoord<Vertice>> CoordMap;
//	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
//		get(boost::vertex_info_3d, vertices_graph);
//	CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());
//	ScaleSpaceBlurrer<VerticesGraph, CoordMap, GaussianKernel<cv::Point3d, double>> scale_space_blurrer;
//	const size_t ext_levels_num = levels_num + 2;
//	const bool kAdditive = false;
//	scale_space_blurrer.Init(kInitCurvSigma, kDMultSigma, ext_levels_num, kAdditive);
//	std::vector<DoubleVertGraphProp> electrostatic_scale_space, steric_scale_space;
//	scale_space_blurrer.MakeScaleSpace(vertices_graph, coord_map, m_gaussian_curvature, ext_levels_num + 1, true, m_gaussian_curvature_scales);
//	scale_space_blurrer.MakeScaleSpace(vertices_graph, coord_map, m_mean_curvature, ext_levels_num + 1, true, m_mean_curvature_scales);
//	scale_space_blurrer.MakeScaleSpace(vertices_graph, coord_map, m_vertex_charge_map, ext_levels_num + 1, true, electrostatic_scale_space);
//	scale_space_blurrer.MakeScaleSpace(vertices_graph, coord_map, m_vertex_lennard_jones_map, ext_levels_num + 1, true, steric_scale_space);
//	//попробуем определить радиус фич
//	std::vector<DoubleVertGraphProp> gaussian_curvature_scales_diff(ext_levels_num, DoubleVertGraphProp(vertices_graph));
//	std::vector<DoubleVertGraphProp> mean_curvature_scales_diff(ext_levels_num , DoubleVertGraphProp(vertices_graph));
//
//	std::vector<DoubleVertGraphProp> electrostatic_scale_space_diff(ext_levels_num, DoubleVertGraphProp(vertices_graph));
//	std::vector<DoubleVertGraphProp> steric_scale_space_diff(ext_levels_num , DoubleVertGraphProp(vertices_graph));
//
//	for (size_t curr_level = 0; curr_level < ext_levels_num; ++curr_level)
//	{
//		gaussian_curvature_scales_diff[curr_level].SetGraph(vertices_graph);
//		mean_curvature_scales_diff[curr_level].SetGraph(vertices_graph);
//		electrostatic_scale_space_diff[curr_level].SetGraph(vertices_graph);
//		steric_scale_space_diff[curr_level].SetGraph(vertices_graph);
//
//		const double coeff = (scale_space_blurrer.GetSigma(curr_level + 1) + scale_space_blurrer.GetSigma(curr_level)) 
//			/ (scale_space_blurrer.GetSigma(curr_level + 1) - scale_space_blurrer.GetSigma(curr_level));
//
//
//		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
//			curr_vertice != end_vertices; ++curr_vertice)
//		{
//			gaussian_curvature_scales_diff[curr_level][*curr_vertice] = coeff * 				 
//				(m_gaussian_curvature_scales[curr_level + 1][*curr_vertice] - m_gaussian_curvature_scales[curr_level][*curr_vertice]);
//			mean_curvature_scales_diff[curr_level][*curr_vertice] = coeff *
//				(m_mean_curvature_scales[curr_level + 1][*curr_vertice] - m_mean_curvature_scales[curr_level][*curr_vertice]);
//			electrostatic_scale_space_diff[curr_level][*curr_vertice] = coeff *
//				(electrostatic_scale_space[curr_level + 1][*curr_vertice] - electrostatic_scale_space[curr_level][*curr_vertice]);
//			steric_scale_space_diff[curr_level][*curr_vertice] = coeff *
//				(steric_scale_space[curr_level + 1][*curr_vertice] - steric_scale_space[curr_level][*curr_vertice]);
//		}	
//	}
//	///////
//	m_vertex_curv_type_mesh_levels.assign(ext_levels_num, VetrticesCurvMap());
//
//	for (size_t curr_level = 0; curr_level < ext_levels_num; ++curr_level)
//	{
//		m_vertex_curv_type_mesh_levels[curr_level].SetGraph(vertices_graph);
//		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
//			curr_vertice != end_vertices; ++curr_vertice)
//		{
//			m_vertex_curv_type_mesh_levels[curr_level][*curr_vertice] = GetCurvatureType(
//				m_gaussian_curvature_scales[curr_level][*curr_vertice], m_mean_curvature_scales[curr_level][*curr_vertice]);
//		}	
//	}	
//
//	std::vector<std::vector<VertexDescriptor>> sing_pts_gaussian, sing_pts_mean;
//	m_gaussian_curvature_scales.resize(ext_levels_num);
//	m_mean_curvature_scales.resize(ext_levels_num);
//
//	FindScaleSingularPointsOnFunc(vertices_graph, coord_map, m_gaussian_curvature_scales, 
//		m_vertex_curv_type_mesh_levels, 1.0, sing_pts_gaussian, false, true);
//	FindScaleSingularPointsOnFunc(vertices_graph, coord_map, m_mean_curvature_scales, 
//		m_vertex_curv_type_mesh_levels, 1.0, sing_pts_mean, false, true);
//	PointsKeeper<VerticesGraph> points_keeper(1.0, vertices_graph);
//	//m_maximums_with_levels.clear();
//
//	//for (size_t curr_level = 0; curr_level < levels_num; ++curr_level)
//	//{
//	//	PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, vertices_graph);
//	//	points_keeper_curr_scale.AddPoints(sing_pts_gaussian[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//	//	points_keeper_curr_scale.AddPoints(sing_pts_mean[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//	//	m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
//	//	points_keeper.AddPoints(points_keeper_curr_scale.Points(), m_vertex_curv_type_mesh_levels[0]);
//	//	std::cout << m_maximums_with_levels[curr_level].size() << "\n";
//	//}
//	//m_maximums = points_keeper.Points();
//	//повторим для разностей
//	m_maximums_with_levels.clear();
//	std::vector<std::vector<VertexDescriptor>> sing_pts_gaussian1, sing_pts_mean1, sing_pts_electric, sing_pts_steric;
//	FindScaleSingularPointsOnFunc(vertices_graph, coord_map, gaussian_curvature_scales_diff, m_vertex_curv_type_mesh_levels, 1.0, sing_pts_gaussian1, true, false);
//	FindScaleSingularPointsOnFunc(vertices_graph, coord_map, mean_curvature_scales_diff, m_vertex_curv_type_mesh_levels, 1.0, sing_pts_mean1, true, false);
//
//	FindScaleSingularPointsOnFunc(vertices_graph, coord_map, electrostatic_scale_space_diff, m_vertex_curv_type_mesh_levels, 1.0, sing_pts_electric, true, false);
//	FindScaleSingularPointsOnFunc(vertices_graph, coord_map, steric_scale_space_diff, m_vertex_curv_type_mesh_levels, 1.0, sing_pts_steric, true, false);
//	for (size_t curr_level = 0; curr_level < levels_num; ++curr_level)
//	{
//		PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, vertices_graph);
//		points_keeper_curr_scale.AddPoints(sing_pts_gaussian1[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//		points_keeper_curr_scale.AddPoints(sing_pts_mean1[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//		std::cout <<points_keeper_curr_scale.Points().size() << " ";
//		points_keeper_curr_scale.AddPoints(sing_pts_electric[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//		points_keeper_curr_scale.AddPoints(sing_pts_steric[curr_level], m_vertex_curv_type_mesh_levels[curr_level + 1]);
//		m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
//		points_keeper.AddPoints(points_keeper_curr_scale.Points(), m_vertex_curv_type_mesh_levels[0]);
//		std::cout <<points_keeper_curr_scale.Points().size() << "\n";
//	}
//	m_maximums = points_keeper.Points();
//	std::cout << "maximums " << m_maximums.size() << " ";	
//}

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