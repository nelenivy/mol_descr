#include "singular_points_finder.h"

#include <algorithm>
#include <utility>
#include "SingularPoints/mesh_types.h"
#include "SingularPoints\mesh_operations.h"
#include "GraphLib\connected_components_segmentator.h"
#include "GraphLib\graph_operations.h"
#include "GraphLib\graph_filter.h"
#include "GraphLib\curvature_calculator.h"
#include "GraphLib\proxy_property_map.h"
#include "GraphLib\graph_functions.h"
#include "CommonUtilities\common_functions.h"
#include "GraphLib\coordinates_transform.h"
#include "GraphLib\graph_dist_calculate.h"
#include "points_keeper.h"
#include "properties_calculator.h"

namespace molecule_descriptor
{

using namespace std;
using namespace cv;

void SngPtsFinderSegmentation::CalcSingPtsFromCalculatedProperties(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
	const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges, 
	const std::vector<std::pair<cv::Point3d, double>>& wdv_radii,
	const bool calc_prop_as_average)
{
	m_mesh_keeper.ConstructMesh(vertices, normals, triangles);
	const Mesh& mesh_to_use = /*m_mesh_keeper.GetMesh();*/GetMesh();
	CalculateDistanceMaps(m_mesh_keeper.GetMesh());
	CalculateCurvature(m_mesh_keeper.GetMesh());
	const int kMaxSegmSize = 500;
	SegmentMolecularSurface(kMaxSegmSize, mesh_to_use);
	FindSegmentsGraphAndCenters(mesh_to_use);
	CalculateAllPotentials(charges, m_mesh_keeper.GetMesh()/*mesh_to_use*/, m_vertex_charge_map);
	CalculateLennardJonesPotentials(wdv_radii, m_mesh_keeper.GetMesh()/*mesh_to_use*/, m_vertex_lennard_jones_map);
	CalculatePropsInSingPts(calc_prop_as_average);
}

void SngPtsFinderSegmentation::CalculatePropsInSingPts(const bool calc_prop_as_average)
{
	CalcPropInSingPts(m_vertex_charge_map, calc_prop_as_average, m_sing_pts_potential);
	CalcPropInSingPts(m_vertex_lennard_jones_map, calc_prop_as_average, m_sing_pts_lennard_jones);
	CalcSegmentsArea();
}

void SngPtsFinderSegmentation::CalcSegmentsArea(/*const Mesh& mesh*/)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;
	const TrianglesGraph& triangles_graph = m_mesh_keeper.GetMesh().triangles;

	//find segments number
	int max_segm = - 1;
	for (auto tr_iter = vertices(triangles_graph).first, end_tr = vertices(triangles_graph).second; tr_iter != end_tr; ++tr_iter)
	{
		CV_Assert(m_triangles_segm[*tr_iter] >= 0);
		max_segm = std::max(max_segm, static_cast<int>(m_triangles_segm[*tr_iter]));
	}
	//calc area for all segments
	std::vector<double> areas(max_segm + 1, 0.0);
	for (auto tr_iter = vertices(triangles_graph).first, end_tr = vertices(triangles_graph).second; tr_iter != end_tr; ++tr_iter)
	{
		const int curr_segm = static_cast<int>(m_triangles_segm[*tr_iter]);
		const MeshTriangle& curr_triangle = get(boost::vertex_info_3d, triangles_graph, *tr_iter);
		areas[curr_segm] += curr_triangle.Area();
	}

	m_sing_pts_segm_area.SetGraph(m_singular_points_graph);
	for (auto vert_iter = vertices(m_singular_points_graph).first, 
		end_vert = vertices(m_singular_points_graph).second; vert_iter != end_vert; ++vert_iter)
	{
		const VertexDescriptor curr_vert = get(boost::vertex_parent, m_singular_points_graph, *vert_iter);
		const int curr_segm = static_cast<int>(m_vertex_segm[*vert_iter]);
		m_sing_pts_segm_area[*vert_iter] = areas[curr_segm];
	}
}
template<typename GraphPropMap>
void SngPtsFinderSegmentation::CalcPropInSingPts(const GraphPropMap& graph_prop_map, const bool calc_prop_as_average, SingPtsDoublePropMap& sing_pts_prop_map)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;
	sing_pts_prop_map.SetGraph(m_singular_points_graph);

	if (calc_prop_as_average)
	{
		std::vector<std::pair<double, bool>> average_prop;
		CalculateAverageProp(vertices_graph, m_vertex_segm, graph_prop_map, average_prop);

		for (auto vert_iter = vertices(m_singular_points_graph).first, 
			end_vert = vertices(m_singular_points_graph).second; vert_iter != end_vert; ++vert_iter)
		{
			const VertexDescriptor curr_vert = get(boost::vertex_parent, m_singular_points_graph, *vert_iter);
			const size_t curr_segm = m_vertex_segm[curr_vert];
			CV_Assert(average_prop[curr_segm].second);
			sing_pts_prop_map[*vert_iter] = average_prop[curr_segm].first;
		}
	}
	else
	{
		for (auto vert_iter = vertices(m_singular_points_graph).first, 
			end_vert = vertices(m_singular_points_graph).second; vert_iter != end_vert; ++vert_iter)
		{
			const VertexDescriptor curr_vert = get(boost::vertex_parent, m_singular_points_graph, *vert_iter);
			sing_pts_prop_map[*vert_iter] = graph_prop_map[curr_vert];
		}
	}
}

void SngPtsFinderSegmentation::CalculateDistanceMaps(const Mesh& mesh)
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
		CoordMap coord_map = GetProxyPropMapVal(coord_3d_map, GetCoord<Vertice>());

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
	const auto coord_map_tr = GetProxyPropMapVal(coord_3d_map_tr, GetCoord<MeshTriangle>());

	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
		get(boost::vertex_info_3d, vertices_graph);
	CoordMap coord_map = GetProxyPropMapVal(coord_3d_map, GetCoord<Vertice>());
	//calculate distances
	m_vert_vert_dist.create(num_vertices(vertices_graph), num_vertices(vertices_graph));
	m_vert_tr_dist.create(num_vertices(vertices_graph), num_vertices(triangles_graph));

	if (0/*m_use_euclid_distance*/)
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
void SngPtsFinderSegmentation::GetNonMarkedSingularPointsLevels(std::pair<std::vector<std::vector<NonMarkedSingularPoint>>, std::vector<std::vector<size_t>>>& non_marked_singular_points)
{
	
}

void SngPtsFinderSegmentation::GetNonMarkedSingularPoints(std::pair<std::vector<NonMarkedSingularPoint>, std::vector<size_t>>& non_marked_singular_points)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;
	const size_t sing_pts_num = num_vertices(m_singular_points_graph);
	non_marked_singular_points.first.resize(sing_pts_num);
	non_marked_singular_points.second.resize(sing_pts_num);

	for (auto vert_iter = vertices(m_singular_points_graph).first, 
		end_vert = vertices(m_singular_points_graph).second; vert_iter != end_vert; ++vert_iter)
	{
		const size_t ind = get(boost::vertex_index, m_singular_points_graph, *vert_iter);
		auto& curr_point = non_marked_singular_points.first[ind];
		
		const auto curr_graph_vertice =  get(boost::vertex_parent, m_singular_points_graph, *vert_iter);
		curr_point.Property().SurfaceType() = m_vertex_curv_type[curr_graph_vertice];
		curr_point.Property().Charge() = Sign(m_sing_pts_potential[*vert_iter]);
		curr_point.Property().ElectricPotential() = m_sing_pts_potential[*vert_iter];
		curr_point.Property().LennardJones() = m_sing_pts_lennard_jones[*vert_iter];
		curr_point.Property().Area() = m_sing_pts_segm_area[*vert_iter];
		curr_point.Coord() = get(boost::vertex_info_3d, vertices_graph, curr_graph_vertice).Center();
		//non_marked_singular_points.second[ind] = m_vertices_type_map[curr_vertice];
	}
}

void SngPtsFinderSegmentation::GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers)
{
	const VerticesGraph& vertices_graph = GetMesh().vertices;
	vertices_with_segm_numbers.clear();
	vertices_with_segm_numbers.reserve(num_vertices(vertices_graph));
	const auto map_3d = get(boost::vertex_info_3d, vertices_graph);
	
	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		const size_t curr_curv = 10.0 * (std::min(1.0, std::max(-1.0, abs(m_gaussian_curvature[*curr_vertice]))) /* + 1.0*/ );
		size_t curr_segm = m_vertex_segm[*curr_vertice];
		const cv::Point3d curr_coord = map_3d[*curr_vertice].Center();
		vertices_with_segm_numbers.push_back(std::make_pair(curr_coord, curr_segm));
	}
}

void SngPtsFinderSegmentation::GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types)
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

void SngPtsFinderSegmentation::Clear()
{
	m_singular_points_graph = SingularPointsGraph();

	//m_vertices_type_map.Clear();
	m_vertex_charge_map.Clear();
	m_vertex_curv_type.Clear();
	m_vertex_segm.Clear();
	m_triangle_curv_type.Clear();
	m_triangles_segm.Clear();
}

void SngPtsFinderSegmentation::CalcTangentBasis(const Mesh& mesh)
{
	//calculate tangent basis
	boost::property_map<const VerticesGraph, boost::vertex_info_3d_t>::const_type coord_3d_map = 
		get(boost::vertex_info_3d, mesh.vertices);
	CoordMap coord_map = GetProxyPropMapVal(coord_3d_map, GetCoord<Vertice>());
	NormalMap norm_map = GetProxyPropMapVal(coord_3d_map, GetNormal<Vertice>());
	m_tangent_basis_map.SetGraph(mesh.vertices);
	CalcTangentCoordSystemMap(mesh.vertices, coord_map, norm_map, 0.01, m_tangent_basis_map);
}


void SngPtsFinderSegmentation::CalculateCurvature(const Mesh& mesh)
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
		if (m_mean_curvature[*curr_vertice] > 0.0)
		{
			m_vertex_curv_type[*curr_vertice] = kConcaveType;
		}
		else if (m_gaussian_curvature[*curr_vertice] < 0.0)
		{
			m_vertex_curv_type[*curr_vertice] = kSaddleType;
		}
		else
		{
			m_vertex_curv_type[*curr_vertice] = kConvexType ;
		}
	}

	//assign triangles labels
	const TrianglesGraph& triangles_graph = mesh.triangles;
	m_triangle_curv_type.SetGraph(triangles_graph);

	for (auto curr_tr_iter = vertices(triangles_graph).first, last_tr = vertices(triangles_graph).second; 
		curr_tr_iter != last_tr; ++curr_tr_iter)
	{
		const MeshTriangle curr_triangle = get(boost::vertex_info_3d, triangles_graph, *curr_tr_iter);
		int8_t types[kSaddleType+1] = {};
		types[m_vertex_curv_type[curr_triangle.GetA()]]++;
		types[m_vertex_curv_type[curr_triangle.GetB()]]++;
		types[m_vertex_curv_type[curr_triangle.GetC()]]++;

		int max_count = -1;
		int max_type = -1;
		for (int curr_type = kConvexType; curr_type <= kSaddleType; ++curr_type)
		{
			if (types[curr_type] > max_count)
			{
				max_count = types[curr_type];
				max_type = curr_type;
			}
		}
		if (max_count >= 2)
		{			
			m_triangle_curv_type[*curr_tr_iter] = max_type;
		}
		else
		{
			m_triangle_curv_type[*curr_tr_iter] = kSaddleType;
		}
	}
}
const size_t kMinElementsInSegm = 5;

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
void SngPtsFinderSegmentation::SegmentMolecularSurface(const size_t max_segm_size, const Mesh& mesh)
{	
	const VerticesGraph& vertices_graph = mesh.vertices;
	const TrianglesGraph& triangles_graph = mesh.triangles;
	m_triangles_segm.SetGraph(triangles_graph);
	m_vertex_segm.SetGraph(vertices_graph);

	//segment with connected components algorithm
	const unsigned char types[] = {kConvexType, kConcaveType};
	size_t all_segm_num = 0;
	ConnectedComponentsSegmentator<TrianglesGraph, TrianglesCurvMap, TrianglesSegmMap> conn_comp_segmentator;

	for (size_t type_ind = 0; type_ind < 2; type_ind++)
	{
		auto& curr_temp_segm = m_type_triangle_segm[type_ind];
		curr_temp_segm.SetGraph(triangles_graph);
		const unsigned char curr_type = types[type_ind];
		size_t curr_segm_num = 0;
		conn_comp_segmentator.SegmentImageValue(triangles_graph, m_triangle_curv_type, curr_type, kMinElementsInSegm, curr_temp_segm, curr_segm_num);

		if (curr_segm_num == 0)
		{
			continue;
		}
		//remap segments number and write to common segments map
		vector<size_t> segments_sizes(curr_segm_num + 1, 0);
		for (auto tr_iter = vertices(triangles_graph).first, last_tr = vertices(triangles_graph).second; 
			tr_iter != last_tr; ++tr_iter)
		{
			if (m_triangle_curv_type[*tr_iter] == curr_type)
			{
				const size_t curr_segm = curr_temp_segm[*tr_iter];
				m_triangles_segm[*tr_iter] = curr_segm + all_segm_num;
				segments_sizes[curr_segm]++;
			}
		}

		//split too big segments
		//find sizes deviation
		double mean_elems_in_segm = 0, dev = 0;
		const std::vector<size_t>::iterator beg_next = segments_sizes.begin() + 1, size_end = segments_sizes.end();
		CalcMeanAndDev(beg_next, size_end, mean_elems_in_segm, dev);
		//split big segments
		const double diff_thresh = 1.75 * dev;
		//const int max_segm_size = static_cast<int>(std::min(2.0 * mean_elems_in_segm, mean_elems_in_segm + 1.75 * dev));

		const size_t segm_num_before_split = curr_segm_num;
		curr_segm_num += all_segm_num;
		for (size_t curr_segm = 1; curr_segm <= segm_num_before_split; ++curr_segm)
		{
			if (segments_sizes[curr_segm] > max_segm_size)
			{
				const size_t parts_of_segm_num = static_cast<int>(ceil(segments_sizes[curr_segm] / static_cast<double>(mean_elems_in_segm)));
				SplitSegment(triangles_graph, curr_segm + all_segm_num, parts_of_segm_num, m_triangles_segm, curr_segm_num);
			}
		}

		all_segm_num = curr_segm_num;
	}

	all_segm_num++;
	//segment saddle vertices with k-means
	int segmented_triangles_num = 0;
	const int saddle_type_ini_segm_num = all_segm_num;
	for (auto tr_iter = vertices(triangles_graph).first, end_iter = vertices(triangles_graph).second; 
		tr_iter != end_iter; ++tr_iter)
	{		
		if (m_triangle_curv_type[*tr_iter] == kSaddleType)
		{
			m_triangles_segm[*tr_iter] = saddle_type_ini_segm_num;//mark for segmenting with k-means
		}
		else
		{
			++segmented_triangles_num;
		}
	}
	//calculate segments num for saddle regions
	const int non_segmented_triangles = num_vertices(triangles_graph) - segmented_triangles_num;
	const size_t parts_num = std::max(Round(all_segm_num * static_cast<double>(non_segmented_triangles) / segmented_triangles_num), 1);
	const size_t segm_num_to_split = saddle_type_ini_segm_num;
	SplitSegment(triangles_graph, segm_num_to_split, parts_num, m_triangles_segm, all_segm_num);
	std::cout << all_segm_num << "\n";
	//write segments_num to vertices
	m_vertex_score_map.SetGraph(vertices_graph);
	for (auto tr_iter = vertices(triangles_graph).first, end_iter = vertices(triangles_graph).second; 
		tr_iter != end_iter; ++tr_iter)
	{
		const MeshTriangle& curr_triangle = get(boost::vertex_info_3d, triangles_graph, *tr_iter);
		const size_t curr_segm = m_triangles_segm[*tr_iter];
		m_vertex_score_map[curr_triangle.GetA()][curr_segm]++;
		m_vertex_score_map[curr_triangle.GetB()][curr_segm]++;
		m_vertex_score_map[curr_triangle.GetC()][curr_segm]++;
	}
	//write values to the vertices segments map
	for (auto vert_iter = vertices(vertices_graph).first, end_iter = vertices(vertices_graph).second; 
		vert_iter != end_iter; ++vert_iter)
	{
		const auto& curr_map = m_vertex_score_map[*vert_iter];
		CV_Assert(!curr_map.empty());
		size_t max_segm = 0;
		uint8_t max_score = 0;
		for (auto segm_iter = curr_map.begin(), end_segm = curr_map.end(); segm_iter != end_segm; ++segm_iter)
		{
			if (segm_iter->second > max_score)
			{
				max_score = segm_iter->second;
				max_segm = segm_iter->first;
			}
		}

		m_vertex_segm[*vert_iter] = max_segm;
	}
}

void SngPtsFinderSegmentation::FindSegmentsGraphAndCenters(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	CalculateSegmentGraphAndCenters(vertices_graph, m_vertex_segm, m_singular_points_graph);
}

}