#include "sng_pts_finder_scale_space.h"

#include <algorithm>
#include <utility>
#include "opencv2/core/core.hpp"
#include "opencv2/core/core_c.h"
#include "SingularPoints/mesh_types.h"
#include "SingularPoints\mesh_operations.h"
#include "GraphLib\connected_components_segmentator.h"
#include "GraphLib\graph_operations.h"
#include "GraphLib\graph_filter.h"
#include "GraphLib\curvature_calculator.h"
#include "GraphLib\graph_functions.h"
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
	ReadParamFromCommandLineWithDefault(argc, argv, "-combine_channels", m_combine_channels, true);

	ReadParamFromCommandLineWithDefault(argc, argv, "-curv_sigma", m_curv_sigma, 0.2);
	ReadParamFromCommandLineWithDefault(argc, argv, "-init_curv_sigma", m_init_curv_sigma, 0.3);
	ReadParamFromCommandLineWithDefault(argc, argv, "-sigma_max", m_sigma_max, 2.65);

	m_scale_space_levels_num = m_sing_pts_levels_num + (m_detect_blobs ? 3 : 2);
	m_diff_between_sing_pts_levels_and_scale_space_levels = 1;

	//INIT SCALE-SPACE
	const bool kAdditive = false;
	const double d_mult_sigma = pow(m_sigma_max / m_init_curv_sigma, 1.0 / (m_scale_space_levels_num));
	std::cout << "d_mult_sigma " << d_mult_sigma << std::endl;
	m_scale_space_blurrer.Init(m_init_curv_sigma, d_mult_sigma, kAdditive);
}

void SngPtsFinderScaleSpace::CalcOnlyProps(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
						   const std::vector<cv::Point3i>& triangles, 
						   const std::vector<std::pair<cv::Point3d, double>>& charges, 
						   const std::vector<std::pair<cv::Point3d, double>>& wdv_radii)
{
	m_mesh_keeper.ConstructMesh(vertices, normals, triangles);
	const Mesh& mesh_to_use = /*m_mesh_keeper.GetMesh();*/GetMesh();
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
	//CalculateVerticesSurfaceType(mesh_to_use);
	CalcCurvature(mesh_to_use);
	CalculateAllPotentials(charges, m_mesh_keeper.GetMesh()/*mesh_to_use*/, m_vertex_charge_map);
	CalculateLennardJonesPotentials(wdv_radii, m_mesh_keeper.GetMesh()/*mesh_to_use*/, m_vertex_lennard_jones_map);

	CalcSingPtsFromCurvatureScales(m_mesh_keeper.GetMesh());
	//CalcShiftsMaximums();
	const int kMaxSegmSize = 500;
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
 + m_diff_between_sing_pts_levels_and_scale_space_levels][curr_descr];
			curr_point.Property().Charge() = Sign(m_vertex_charge_map[curr_descr]);
			curr_point.Property().ElectricPotential() = m_output_scale_space[level][ELECTR_POTENT][curr_descr]
			* m_mean_and_sigma[ELECTR_POTENT][1]+ m_mean_and_sigma[ELECTR_POTENT][0] - 3 * m_mean_and_sigma[ELECTR_POTENT][1];
			
			curr_point.Property().LennardJones() = m_output_scale_space[level][STERIC_POTENT][curr_descr]
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
	vertices_with_types_lev.resize(m_vertex_curv_type_mesh_levels.size() - 2);

	for (size_t lev = 0; lev < vertices_with_types_lev.size() - 2; ++lev)
	{
		vertices_with_types_lev[lev].reserve(num_vertices(vertices_graph));

		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			const size_t curr_curv_type = m_vertex_curv_type_mesh_levels[lev + 1][*curr_vertice];
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

void SngPtsFinderScaleSpace::CalcShiftsMaximums()
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;
	PointsKeeper<VerticesGraph> points_keeper(1.0, vertices_graph);
	CurvatureCalculator<VerticesGraph> curvature_calculator(10);
	m_maximums.clear();
	m_maximums_with_levels.clear()/*resize(filtered_curvatures_0.size())*/;
	m_vertex_curv_type_mesh_levels.assign(m_filtered_mesh_levels.size(), VetrticesCurvMap(vertices_graph));

	for (size_t mesh_level = 0; mesh_level < m_filtered_mesh_levels.size(); ++mesh_level)
	{
		std::cout << "mesh lev " << mesh_level << "\n";
		const VerticesGraph& curr_mesh_level = m_filtered_mesh_levels[mesh_level].vertices;
		const VerticesGraph& curr_mesh_level_2 = m_filtered_mesh_levels_2[mesh_level].vertices;
		DoubleVertGraphProp curr_shifts(vertices_graph);

		for (auto curr_vertice = vertices(curr_mesh_level).first, end_vertices = vertices(curr_mesh_level).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			const Point3d diff = get(boost::vertex_info_3d, curr_mesh_level, *curr_vertice).Center() - 
				get(boost::vertex_info_3d, curr_mesh_level_2, *curr_vertice).Center();
			const Point3d norm = get(boost::vertex_info_3d, curr_mesh_level, *curr_vertice).Normal();
			curr_shifts[*curr_vertice] = abs(norm.x * diff.x + norm.y * diff.y + norm.z * diff.z);

		}
		for (auto curr_vertice = vertices(curr_mesh_level).first, end_vertices = vertices(curr_mesh_level).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			Curvature curvatures;
			curvature_calculator.CalculateCurvatureCubic(curr_mesh_level, get(boost::vertex_info_3d, curr_mesh_level), *curr_vertice, curvatures);
			m_vertex_curv_type_mesh_levels[mesh_level][*curr_vertice] = GetCurvatureType(curvatures.gaussian_curv, curvatures.mean_curv);
		}
		typedef MedianKernel<double, size_t> KernelType1;
		DoubleVertGraphProp temp0(vertices_graph);
		FilterGraphEdgeDist(KernelType1(1), vertices_graph, curr_shifts, temp0);
		DoubleVertGraphProp temp00(curr_mesh_level);
		DoubleVertGraphProp temp11(curr_mesh_level);
		 
		PointsKeeper<VerticesGraph> points_keeper_mesh_level(1.0 /** sqrt(mesh_level + 1.0)*/, vertices_graph);

			std::vector<VertexDescriptor> maximums_0;
			FindLocalMaximumsOfAbsVal(curr_mesh_level/*vertices_graph*/, temp0, maximums_0, std::greater_equal<double>(), std::less_equal<double>());
			PointsKeeper<VerticesGraph> points_keeper_level(1.0, vertices_graph);
			points_keeper_level.AddPoints(maximums_0, m_vertex_curv_type_mesh_levels[mesh_level]);
			points_keeper_mesh_level.AddPoints(points_keeper_level.Points(), m_vertex_curv_type_mesh_levels[mesh_level]);
			//m_maximums_with_levels[curr_level] = points_keeper_level.Points();
			std::cout << maximums_0.size() << "\n";
			points_keeper.AddPoints(maximums_0, m_vertex_curv_type_mesh_levels[mesh_level]);
			//m_maximums.insert(m_maximums.end(), maximums_0.begin(), maximums_0.end());
			//m_maximums.insert(m_maximums.end(), maximums_1.begin(), maximums_1.end());
		m_maximums_with_levels.push_back(points_keeper_mesh_level.Points());
		std::cout << m_maximums_with_levels.back().size() << "\n";
	}
}



void SngPtsFinderScaleSpace::CalculateVerticesSurfaceType(const Mesh& mesh)
{
	
	//std::vector<DoubleVertGraphProp> filtered_gaussian_curvatures(kCurvLevels, DoubleVertGraphProp(vertices_graph));
	//std::vector<DoubleVertGraphProp> filtered_mean_curvatures(kCurvLevels, DoubleVertGraphProp(vertices_graph));

	//m_mean_curvature.SetGraph(vertices_graph);
	//m_gaussian_curvature.SetGraph(vertices_graph);
	//const int kMaxNeighbNum = 10;
	//CurvatureCalculator<VerticesGraph> curvature_calculator(kMaxNeighbNum);
	////typedef GaussianKernel<cv::Point3d, double> KernelType;
	//const size_t kMedRad = 5;
	//DoubleVertGraphProp gauss_curv_temp(vertices_graph);
	//DoubleVertGraphProp mean_curv_temp(vertices_graph);

	//PointsKeeper<VerticesGraph> points_keeper(1.0, vertices_graph);
	//m_maximums.clear();
	//m_maximums_with_levels.clear();
	//m_vertex_curv_type_mesh_levels.assign(m_filtered_mesh_levels.size(), VetrticesCurvMap(vertices_graph));

	//for (size_t mesh_level = 0; mesh_level < m_filtered_mesh_levels.size(); ++mesh_level)
	//{
	//	std::cout << "mesh lev " << mesh_level << "\n";
	//	const VerticesGraph& curr_mesh_level = m_filtered_mesh_levels[mesh_level].vertices;
	//	for (auto curr_vertice = vertices(curr_mesh_level).first, end_vertices = vertices(curr_mesh_level).second; 
	//		curr_vertice != end_vertices; ++curr_vertice)
	//	{
	//		Curvature curvatures;
	//		curvature_calculator.CalculateCurvatureCubic(curr_mesh_level, get(boost::vertex_info_3d, curr_mesh_level), *curr_vertice, curvatures);
	//		m_gaussian_curvature[*curr_vertice] = curvatures.gaussian_curv;
	//		m_mean_curvature[*curr_vertice] = curvatures.mean_curv;

	//	}
	//	typedef MedianKernel<double, size_t> KernelType1;
	//	FilterGraphEdgeDist(KernelType1(1), curr_mesh_level, m_gaussian_curvature, gauss_curv_temp);
	//	FilterGraphEdgeDist(KernelType1(1), curr_mesh_level, m_mean_curvature, mean_curv_temp);
	//	DoubleVertGraphProp gauss_curv_temp_1(curr_mesh_level);
	//	DoubleVertGraphProp mean_curv_temp_1(curr_mesh_level);

	//	for (size_t curr_level = 0; curr_level < kCurvLevels; ++curr_level)
	//	{
	//		typedef GaussianKernel<cv::Point3d, double> KernelType;
	//		FilterGraphDist(KernelType(curr_level * 0.2 + 0.4), curr_mesh_level/*vertices_graph*/, GetProxyPropMap(get(boost::vertex_info_3d, curr_mesh_level/*vertices_graph*/), GetCoord<Vertice>()), 
	//			/*m_curvature_0*/ gauss_curv_temp, gauss_curv_temp_1/*filtered_curvatures_0[curr_level]*/);///////////
	//		FilterGraphDist(KernelType(curr_level * 0.2 + 0.4), curr_mesh_level/*vertices_graph*/, GetProxyPropMap(get(boost::vertex_info_3d, curr_mesh_level/*vertices_graph*/), GetCoord<Vertice>()),
	//			/*m_curvature_1*/ mean_curv_temp, mean_curv_temp_1/*filtered_curvatures_1[curr_level]*/);	
	//		FilterGraphEdgeDist(KernelType1(1), curr_mesh_level/*vertices_graph*/, gauss_curv_temp_1, filtered_gaussian_curvatures[curr_level]/*, temp0*/);
	//		FilterGraphEdgeDist(KernelType1(1), curr_mesh_level/*vertices_graph*/, mean_curv_temp_1, filtered_mean_curvatures[curr_level]/*, temp1*/);

	//		if (curr_level == 0)
	//		{
	//			for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
	//				curr_vertice != end_vertices; ++curr_vertice)
	//			{
	//				m_vertex_curv_type_mesh_levels[mesh_level][*curr_vertice] = GetCurvatureType(
	//					filtered_gaussian_curvatures[curr_level][*curr_vertice], filtered_mean_curvatures[curr_level][*curr_vertice]);
	//			}
	//		}
	//	
	//	}
	//	//m_curvature_1 = filtered_curvatures_1[0];
	//	//m_curvature_0 = filtered_curvatures_2[0];
	//	PointsKeeper<VerticesGraph> points_keeper_mesh_level(1.0 * sqrt(mesh_level + 1.0), vertices_graph);

	//	for (size_t curr_level = 0; curr_level < filtered_gaussian_curvatures.size(); ++curr_level)
	//	{
	//		std::vector<VertexDescriptor> gauss_curv_maximums, mean_curv_maximums;
	//		FindLocalMaximumsOfAbsVal(curr_mesh_level/*vertices_graph*/, filtered_gaussian_curvatures[curr_level], gauss_curv_maximums, std::greater<double>());
	//		FindLocalMaximumsOfAbsVal(curr_mesh_level/*vertices_graph*/, filtered_mean_curvatures[curr_level], mean_curv_maximums, std::greater<double>());
	//		PointsKeeper<VerticesGraph> points_keeper_level(1.0, vertices_graph);
	//		points_keeper_level.AddPoints(gauss_curv_maximums, m_vertex_curv_type_mesh_levels[mesh_level]);
	//		points_keeper_level.AddPoints(mean_curv_maximums, m_vertex_curv_type_mesh_levels[mesh_level]);
	//		points_keeper_mesh_level.AddPoints(points_keeper_level.Points(), m_vertex_curv_type_mesh_levels[mesh_level]);
	//		//m_maximums_with_levels[curr_level] = points_keeper_level.Points();
	//		std::cout << gauss_curv_maximums.size() << " " << mean_curv_maximums.size() << "\n";
	//		points_keeper.AddPoints(gauss_curv_maximums, m_vertex_curv_type_mesh_levels[mesh_level]);
	//		points_keeper.AddPoints(mean_curv_maximums, m_vertex_curv_type_mesh_levels[mesh_level]);
	//		//m_maximums.insert(m_maximums.end(), maximums_0.begin(), maximums_0.end());
	//		//m_maximums.insert(m_maximums.end(), maximums_1.begin(), maximums_1.end());
	//	}
	//	m_maximums_with_levels.push_back(points_keeper_mesh_level.Points());
	//	std::cout << m_maximums_with_levels.back().size() << "\n";
	//}
	///*for (size_t curr_level = 0; curr_level < filtered_curvatures_0.size(); ++curr_level)
	//{
	//	PointsKeeper<VerticesGraph> points_keeper_level(0.5, vertices_graph);
	//	for (size_t next_levels = curr_level; next_levels < filtered_curvatures_0.size(); ++next_levels)
	//	{
	//		points_keeper_level.AddPoints(m_maximums_with_levels[next_levels]);
	//	}
	//	m_maximums_with_levels[curr_level] = points_keeper_level.Points();
	//}*/
	//m_maximums = points_keeper.Points();
	/*std::vector<DoubleVertGraphProp> diff_0(filtered_curvatures_0.size() - 1, DoubleVertGraphProp(vertices_graph));
	std::vector<DoubleVertGraphProp> diff_1(filtered_curvatures_1.size() - 1, DoubleVertGraphProp(vertices_graph));

	for (size_t curr_level = 0; curr_level < filtered_curvatures_0.size() - 1; ++curr_level)
	{
		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			diff_0[curr_level][*curr_vertice] = abs( filtered_curvatures_0[curr_level + 1][*curr_vertice]  - 
				filtered_curvatures_0[curr_level][*curr_vertice] ) / (filtered_curvatures_0[curr_level][*curr_vertice] + 0.001);
			diff_1[curr_level][*curr_vertice] = abs( filtered_curvatures_1[curr_level + 1][*curr_vertice]  -
				filtered_curvatures_1[curr_level][*curr_vertice] ) / (filtered_curvatures_1[curr_level][*curr_vertice] + 0.001);		
		}
	}
	FindLocalMaximumsOnLevels(vertices_graph, diff_0, m_maximums);
	std::vector<VertexDescriptor> maximums_1;
	FindLocalMaximumsOnLevels(vertices_graph, diff_1, maximums_1);
	m_maximums.insert(m_maximums.end(), maximums_1.begin(), maximums_1.end());*/
	std::cout << "maximums " << m_maximums.size() << " ";
	//filter curvature

	//ContPropMap<VerticesGraph, std::vector<double>, VERTEX> filtered_curvature_1(vertices_graph);
	//ContPropMap<VerticesGraph, std::vector<double>, VERTEX> filtered_curvature_0(vertices_graph);
	/*typedef MedianKernel<double, size_t> KernelType;
	const size_t kMedRad = 1;
	FilterGraphEdgeDist(KernelType(kMedRad), vertices_graph, m_curvature_1, filtered_curvature_1);
	FilterGraphEdgeDist(KernelType(kMedRad), vertices_graph, m_curvature_0, filtered_curvature_0);*/
	//m_curvature_1 = filtered_curvature_1;
	//m_curvature_0 = filtered_curvature_0;
	/*typedef GaussianKernel<cv::Point3d, double> KernelType; 

	GraphDistFilter<VerticesGraph, ProxyCoordMapGraph, GaussianKernel<cv::Point3d, double>> gauss_filter_1(GaussianKernel<cv::Point3d, double>(0.5));
	ContPropMap<VerticesGraph, std::vector<double>, VERTEX> filtered_curvature_1(vertices_graph);
	ContPropMap<VerticesGraph, std::vector<double>, VERTEX> filtered_curvature_0(vertices_graph);
	gauss_filter_1.Filter(vertices_graph, map_3d_coord, curvature_1, filtered_curvature_1);
	gauss_filter_1.Filter(vertices_graph, map_3d_coord, curvature_0, filtered_curvature_0);*/
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
		curvature_calculator.CalculateCurvatureCubic(vertices_graph, get(boost::vertex_info_3d, vertices_graph), *curr_vertice, curvatures);
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

void SngPtsFinderScaleSpace::CalcSingPtsFromCurvatureScales(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;	

	//DATA PREPARATION
	std::vector<DoubleVertGraphProp> input_prop_map(ISingularPointsFinder::PROPS_NUM);
	input_prop_map[ISingularPointsFinder::GAUSS_CURV] = m_gaussian_curvature;
	input_prop_map[ISingularPointsFinder::MEAN_CURV] = m_mean_curvature;
	input_prop_map[ISingularPointsFinder::ELECTR_POTENT] = m_vertex_charge_map;
	input_prop_map[ISingularPointsFinder::STERIC_POTENT] = m_vertex_lennard_jones_map;

	if (!m_mean_and_sigma.empty())
	{//rescale functions
		for (int curr_prop = ISingularPointsFinder::FIRST_PROP; curr_prop < ISingularPointsFinder::PROPS_NUM; ++ curr_prop)
		{
			const double min_val = m_mean_and_sigma[curr_prop][0] - 3.0 * m_mean_and_sigma[curr_prop][1];
			const double max_val = m_mean_and_sigma[curr_prop][0] + 3.0 * m_mean_and_sigma[curr_prop][1];

			for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
				curr_vertice != end_vertices; ++curr_vertice)
			{
				input_prop_map[curr_prop][*curr_vertice] = std::min(input_prop_map[curr_prop][*curr_vertice], max_val);
				input_prop_map[curr_prop][*curr_vertice] = std::max(input_prop_map[curr_prop][*curr_vertice], min_val);
				input_prop_map[curr_prop][*curr_vertice] -= min_val;
				input_prop_map[curr_prop][*curr_vertice] /= m_mean_and_sigma[curr_prop][1];
			}
		}
	}
	//PREPAIR FOR MAKING SCALE SPACE
	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
		get(boost::vertex_info_3d, vertices_graph);
	CoordMap coord_map = GetProxyPropMap(coord_3d_map, GetCoord<Vertice>());
	//MAKE SCALE-SPACE
	const bool use_post_filter = !m_detect_blobs;
	m_scale_space_blurrer.MakeScaleSpace(vertices_graph, coord_map, input_prop_map, 
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
	//CALCULATE MEAN IN NEIGHB
	std::vector<std::vector<DoubleVertGraphProp>> scale_space_vect_to_project(m_scale_space_levels_num);

	struct DataFiller
	{
		void operator()(const CoordMap& coord, const double coord_scale, 
			const std::vector<DoubleVertGraphProp>& props, const double props_scale,
			const VertexDescriptor curr_vert, const int data_row, Mat_<double>& data) const
		{
			data(data_row, 0) = coord_scale * coord[curr_vert].x;
			data(data_row, 1) = coord_scale * coord[curr_vert].y;
			data(data_row, 2) = coord_scale * coord[curr_vert].z;

			for (int col = 3; col < 3 + PROPS_NUM; ++col)
			{
				data(data_row, col) = props_scale * props[col - 3][curr_vert];
			}
		}
	};
	for (size_t curr_level = 0; curr_level < m_scale_space_levels_num; ++curr_level)
	{
		scale_space_vect_to_project[curr_level].resize(m_output_scale_space[curr_level].size());
		for (size_t curr_prop = 0; curr_prop < m_output_scale_space[curr_level].size(); ++curr_prop)
		{
			scale_space_vect_to_project[curr_level][curr_prop].SetGraph(vertices_graph);
		}

		for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
			curr_vertice != end_vertices; ++curr_vertice)
		{
			bool success = false;
			const int neighb_num = adjacent_vertices(*curr_vertice, vertices_graph).second - adjacent_vertices(*curr_vertice, vertices_graph).first;
			cv::Mat_<double> data(neighb_num + 1, 3 + PROPS_NUM);
			DataFiller()(coord_map, 1.0, m_output_scale_space[curr_level], 1.0, *curr_vertice, 0, data);

			for (auto neighb_it = adjacent_vertices(*curr_vertice, vertices_graph).first,
				end_neighb = adjacent_vertices(*curr_vertice, vertices_graph).second; neighb_it != end_neighb; ++neighb_it)
			{
				DataFiller()(coord_map, 1.0, m_output_scale_space[curr_level], 1.0, *neighb_it, end_neighb - neighb_it, data);
			}
				
			cv::PCA pca(data, cv::noArray(), CV_PCA_DATA_AS_ROW, 4);
			cv::Mat_<double> data_basis(1, 3 + PROPS_NUM), proj_basis;
			DataFiller()(coord_map, 1.0, m_output_scale_space[curr_level], 0.0, *curr_vertice, 0, data_basis);
			cv::Mat_<double> curr_data(1, 3 + PROPS_NUM), proj_curr;
			DataFiller()(coord_map, 1.0, m_output_scale_space[curr_level], 1.0, *curr_vertice, 0, curr_data);
			pca.project(data_basis, proj_basis);
			pca.project(curr_data, proj_curr);
			cv::Mat_<double> vect_to_project;
			pca.backProject(proj_curr - proj_basis, vect_to_project);

			for (int prop = 0; prop < PROPS_NUM; ++prop)
			{
				scale_space_vect_to_project[curr_level][prop][*curr_vertice] = vect_to_project(0, prop + 3);
			}
		}
	}
	//DETECT FEATURES
	if (m_detect_blobs)
	{
		std::vector<std::vector<DoubleVertGraphProp>> output_scale_space_diff(m_scale_space_levels_num - 1);

		for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
		{
			const double coeff = (m_scale_space_blurrer.GetSigma(curr_level + 1) + m_scale_space_blurrer.GetSigma(curr_level)) / 2.0 
				/ (m_scale_space_blurrer.GetSigma(curr_level + 1) - m_scale_space_blurrer.GetSigma(curr_level));

			output_scale_space_diff[curr_level].resize(m_output_scale_space[curr_level].size());

			for (size_t curr_prop = 0; curr_prop < output_scale_space_diff[curr_level].size(); ++curr_prop)
			{
				output_scale_space_diff[curr_level][curr_prop].SetGraph(vertices_graph);
				for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
					curr_vertice != end_vertices; ++curr_vertice)
				{
					output_scale_space_diff[curr_level][curr_prop][*curr_vertice] = coeff * 				 
						(m_output_scale_space[curr_level + 1][curr_prop][*curr_vertice] - m_output_scale_space[curr_level][curr_prop][*curr_vertice]);

				}	
			}
		}

		m_maximums_with_levels.clear();

		if (m_combine_channels)
		{
			FindLocalMaximumsOnLevelsVect(vertices_graph, output_scale_space_diff, scale_space_vect_to_project/*output_scale_space*/,  
				m_maximums_with_levels, std::greater<double>(), std::less<double>(), true);
		}
		else
		{
			std::vector<DoubleVertGraphProp> gaussian_curvature_scales_diff(m_scale_space_levels_num - 1, DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> mean_curvature_scales_diff(m_scale_space_levels_num - 1 , DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> electrostatic_scale_space_diff(m_scale_space_levels_num - 1, DoubleVertGraphProp(vertices_graph));
			std::vector<DoubleVertGraphProp> steric_scale_space_diff(m_scale_space_levels_num - 1 , DoubleVertGraphProp(vertices_graph));

			for (size_t curr_level = 0; curr_level < m_scale_space_levels_num - 1; ++curr_level)
			{
				gaussian_curvature_scales_diff[curr_level] = output_scale_space_diff[curr_level][GAUSS_CURV];
				mean_curvature_scales_diff[curr_level] = output_scale_space_diff[curr_level][MEAN_CURV];
				electrostatic_scale_space_diff[curr_level] = output_scale_space_diff[curr_level][ELECTR_POTENT];
				steric_scale_space_diff[curr_level] = output_scale_space_diff[curr_level][STERIC_POTENT];
			}
			m_vertex_curv_type_mesh_levels.resize(m_scale_space_levels_num - 1);
			std::vector<std::vector<VertexDescriptor>> sing_pts_gaussian, sing_pts_mean, sing_pts_electric, sing_pts_steric;
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, gaussian_curvature_scales_diff, 
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_gaussian, true, false);
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, mean_curvature_scales_diff, 
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_mean, true, false);

			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, electrostatic_scale_space_diff, 
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_electric, true, false);
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, steric_scale_space_diff, 
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_steric, true, false);

			for (size_t curr_level = 0; curr_level < m_sing_pts_levels_num; ++curr_level)
			{
				PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, vertices_graph);
				points_keeper_curr_scale.AddPoints(sing_pts_gaussian[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
				points_keeper_curr_scale.AddPoints(sing_pts_mean[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
				std::cout <<points_keeper_curr_scale.Points().size() << " ";
				points_keeper_curr_scale.AddPoints(sing_pts_electric[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
				points_keeper_curr_scale.AddPoints(sing_pts_steric[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
				m_maximums_with_levels.push_back(points_keeper_curr_scale.Points());
				std::cout <<points_keeper_curr_scale.Points().size() << "\n";
			}
		}
	}
	else
	{
		m_maximums_with_levels.clear();

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
					for (int prop = 0; prop < PROPS_NUM; ++prop)
					{
						scale_space_projected[curr_level][*curr_vertice] += 
							m_output_scale_space[curr_level][prop][*curr_vertice];
							///*scale_space_vect_to_project*/m_output_scale_space[curr_level][prop][*curr_vertice];
						norm += scale_space_vect_to_project[curr_level][prop][*curr_vertice] * 
							scale_space_vect_to_project[curr_level][prop][*curr_vertice];
					}
					//scale_space_projected[curr_level][*curr_vertice] /= sqrt(norm);
				}
			}

			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, scale_space_projected,
				m_vertex_curv_type_mesh_levels, 1.0, m_maximums_with_levels, false, true);
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
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_gaussian, false, true);
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, mean_curvature_scales,
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_mean, false, true);
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, electrostatic_scale_space,
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_electric, false, true);
			FindScaleSingularPointsOnFunc(vertices_graph, coord_map, steric_scale_space,
				m_vertex_curv_type_mesh_levels, 1.0, sing_pts_steric, false, true);

			for (size_t curr_level = 0; curr_level < m_sing_pts_levels_num; ++curr_level)
			{
				PointsKeeper<VerticesGraph> points_keeper_curr_scale(1.0, vertices_graph);
				points_keeper_curr_scale.AddPoints(sing_pts_gaussian[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
				points_keeper_curr_scale.AddPoints(sing_pts_mean[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
				std::cout <<points_keeper_curr_scale.Points().size() << " ";
				points_keeper_curr_scale.AddPoints(sing_pts_electric[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
				points_keeper_curr_scale.AddPoints(sing_pts_steric[curr_level], m_vertex_curv_type_mesh_levels[curr_level 
					+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
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
			+ m_diff_between_sing_pts_levels_and_scale_space_levels]);
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