#include "singular_points_finder.h"

#include <algorithm>

#include "SingularPoints/mesh_types.h"
#include "SingularPoints\mesh_operations.h"
#include "GraphLib\connected_components_segmentator.h"
#include "GraphLib\graph_operations.h"
#include "GraphLib\graph_filter.h"
#include "GraphLib\curvature_calculator.h"
#include "GraphLib\proxy_property_map.h"
#include "GraphLib\graph_functions.h"
#include "CommonUtilities\common_functions.h"

namespace molecule_descriptor
{

using namespace std;
using namespace cv;

std::shared_ptr<ISingularPointsFinder> CreateSingularPointsFinder()
{
	ISingularPointsFinder* instance = new SingularPointsFinder;
	return std::shared_ptr<ISingularPointsFinder>(instance, ReleaseDeleter());
}

cv::Point3d CalculatePotential(const cv::Point3d& point, const std::vector<std::pair<cv::Point3d, double>>& charges)
{
	cv::Point3d force = 0.0;

	for (auto charges_iter = charges.begin(); charges_iter != charges.end(); ++charges_iter)
	{
		const cv::Point3d curr_vect = charges_iter->first - point;
		const double inverse_dist = 1.0 / (cv::norm(curr_vect) + 0.00000001);
		const cv::Point3d curr_direction = curr_vect * inverse_dist;
		force += curr_direction * charges_iter->second * (inverse_dist * inverse_dist);
	}

	return force;
}

double CalculateLennardJonesPotential(const cv::Point3d& point, const std::vector<std::pair<cv::Point3d, double>>& wdv_radii, const double probe_radius)
{
	double potential = 0.0;

	for (auto radius_iter = wdv_radii.begin(); radius_iter != wdv_radii.end(); ++radius_iter)
	{
		const cv::Point3d curr_vect = radius_iter->first - point;
		const double inverse_dist = 1.0 / (cv::norm(curr_vect) + 0.00000001);
		const double sigma = (radius_iter->second + probe_radius) / 2.0;
		potential += pow(sigma * inverse_dist, 12.0) - pow(sigma * inverse_dist, 6.0);
	}

	return potential;
}

void SingularPointsFinder::CalculateAllPotentials(const std::vector<std::pair<cv::Point3d, double>>& charges, const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	m_vertex_charge_map.SetGraph(vertices_graph);
	const double kProbeRadius = 1.53;
	for (auto vertex_iter = vertices(vertices_graph).first, 
		end_iter = vertices(vertices_graph).second; vertex_iter != end_iter; ++vertex_iter)
	{
		const Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Center(); 
		const Point3d curr_norm = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Normal();
		const cv::Point3d curr_force = CalculatePotential(curr_coord + curr_norm * kProbeRadius, charges);
		m_vertex_charge_map[*vertex_iter] = curr_force.x * curr_norm.x + curr_force.y * curr_norm.y + curr_force.z * curr_norm.z;
	}
}

void SingularPointsFinder::CalculateLennardJonesPotentials(const std::vector<std::pair<cv::Point3d, double>>& wdv_radii, const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	m_vertex_lennard_jones_map.SetGraph(vertices_graph);
	const double kProbeRadius = 1.53;

	for (auto vertex_iter = vertices(vertices_graph).first, 
		end_iter = vertices(vertices_graph).second; vertex_iter != end_iter; ++vertex_iter)
	{
		const Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Center(); 
		const Point3d curr_norm = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Normal();
		const cv::Point3d curr_force = CalculateLennardJonesPotential(curr_coord + curr_norm * kProbeRadius, wdv_radii, kProbeRadius);
		m_vertex_lennard_jones_map[*vertex_iter] = curr_force.x * curr_norm.x + curr_force.y * curr_norm.y + curr_force.z * curr_norm.z;
	}
}

void SingularPointsFinder::Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
	const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges, 
	const std::vector<std::pair<cv::Point3d, double>>& wdv_radii,
	const bool calc_prop_as_average)
{
	m_mesh_keeper.ConstructMesh(vertices, normals, triangles);
	FilterMesh(m_mesh_keeper.GetMesh(), GaussianKernel<cv::Point3d, cv::Point3d>(1.0), m_filtered_mesh);
	const Mesh& mesh_to_use = /*m_mesh_keeper.GetMesh();*/GetMesh();
	CalculateVerticesSurfaceType(mesh_to_use);
	const int kMaxSegmSize = 500;
	SegmentMolecularSurface(kMaxSegmSize, mesh_to_use);
	FindSegmentsGraphAndCenters(mesh_to_use);
	CalculateAllPotentials(charges, m_mesh_keeper.GetMesh()/*mesh_to_use*/);
	CalculateLennardJonesPotentials(wdv_radii, m_mesh_keeper.GetMesh()/*mesh_to_use*/);
	CalculateSingularPointsHistograms();
	CalculatePropsInSingPts(calc_prop_as_average);
}

void SingularPointsFinder::CalculatePropsInSingPts(const bool calc_prop_as_average)
{
	CalcPropInSingPts(m_vertex_charge_map, calc_prop_as_average, m_sing_pts_potential);
	CalcPropInSingPts(m_vertex_lennard_jones_map, calc_prop_as_average, m_sing_pts_lennard_jones);
}

template<typename GraphPropMap>
void SingularPointsFinder::CalcPropInSingPts(const GraphPropMap& graph_prop_map, const bool calc_prop_as_average, SingPtsDoublePropMap& sing_pts_prop_map)
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

const size_t kElectricSignMax = 2;
const size_t kSurfaceTypeMax = kSaddleType;

void SingularPointsFinder::GetMarkedSingularPoints(std::vector<MarkedSingularPoint>& marked_singular_points)
{
	marked_singular_points.clear();
	marked_singular_points.reserve(m_singular_points.size());
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;

	/*std::vector<VertexDescriptor> centers;
	FindLocalMaximums(vertices_graph, m_curvature_1, centers);

	for (size_t ind = 0; ind < centers.size(); ++ind)
	{
		const auto curr_vertice = centers[ind];
		const size_t curr_type = 1;
		const cv::Point3d& curr_coord = get(boost::vertex_info_3d, vertices_graph, curr_vertice).Center();
		marked_singular_points.emplace_back(curr_coord, curr_type);
	}*/

	/*for (size_t ind = 0; ind < m_singular_points.size(); ++ind)
	{
		const auto curr_vertice = m_singular_points[ind];
		const size_t curr_type = m_vertices_type_map[curr_vertice];
		const cv::Point3d& curr_coord = get(boost::vertex_info_3d, vertices_graph, curr_vertice).Center();
		marked_singular_points.emplace_back(curr_coord, curr_type);
	}*/
}

void SingularPointsFinder::GetNonMarkedSingularPoints(std::pair<std::vector<NonMarkedSingularPoint>, std::vector<size_t>>& non_marked_singular_points)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;
	const size_t sing_pts_num = num_vertices(m_singular_points_graph);
	non_marked_singular_points.first.resize(sing_pts_num);
	non_marked_singular_points.second.resize(sing_pts_num);

	for (auto vert_iter = vertices(m_singular_points_graph).first, 
		end_vert = vertices(m_singular_points_graph).second; vert_iter != end_vert; ++vert_iter)
	{
		const int ind = get(boost::vertex_index, m_singular_points_graph, *vert_iter);
		auto& curr_point = non_marked_singular_points.first[ind];
		
		const auto curr_graph_vertice =  get(boost::vertex_parent, m_singular_points_graph, *vert_iter);
		curr_point.Property().SurfaceType() = m_vertex_curv_type[curr_graph_vertice];
		curr_point.Property().Charge() = Sign(m_sing_pts_potential[*vert_iter]);
		curr_point.Property().ElectricPotential() = m_sing_pts_potential[*vert_iter];
		curr_point.Property().LennardJones() = m_sing_pts_lennard_jones[*vert_iter];
		curr_point.Coord() = get(boost::vertex_info_3d, vertices_graph, curr_graph_vertice).Center();
		//non_marked_singular_points.second[ind] = m_vertices_type_map[curr_vertice];
	}
}

void SingularPointsFinder::GetSingularPointsHisto(std::vector<HistogramSingularPoint<kHistSize>>& singular_points_hist)
{
	singular_points_hist.resize(num_vertices(m_singular_points_graph));
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMesh().vertices;
	const auto coord_map = get(boost::vertex_info_3d, vertices_graph);
	int ind = 0;
	for (auto curr_center = vertices(m_singular_points_graph).first, end_centers = vertices(m_singular_points_graph).second; 
		curr_center != end_centers; ++curr_center, ++ind)
	{
		auto& curr_point = singular_points_hist[ind];
		curr_point.Property() = m_sing_pts_histo[*curr_center];
		curr_point.Coord() = coord_map[get(boost::vertex_parent, m_singular_points_graph, *curr_center)].Center();
	}
}

void SingularPointsFinder::GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers)
{
	const VerticesGraph& vertices_graph = GetMesh().vertices;
	vertices_with_segm_numbers.clear();
	vertices_with_segm_numbers.reserve(num_vertices(vertices_graph));
	const auto map_3d = get(boost::vertex_info_3d, vertices_graph);
	
	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		const size_t curr_segm = m_vertex_segm[*curr_vertice];
		//const size_t curr_segm = 10.0 * (std::min(1.0, std::max(-1.0, abs(m_curvature_0[*curr_vertice]))) /* + 1.0*/ ); 
		const cv::Point3d curr_coord = map_3d[*curr_vertice].Center();
		vertices_with_segm_numbers.push_back(std::make_pair(curr_coord, curr_segm));
	}
}

void SingularPointsFinder::GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types)
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

void SingularPointsFinder::Clear()
{
	m_singular_points.clear();
	m_singular_points_graph = SingularPointsGraph();
	m_sing_pts_histo.Clear();

	//m_vertices_type_map.Clear();
	m_vertex_charge_map.Clear();
	m_vertex_curv_type.Clear();
	m_vertex_segm.Clear();
	m_triangle_curv_type.Clear();
	m_triangles_segm.Clear();
}

void SingularPointsFinder::CalculateVerticesSurfaceType(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	m_vertex_curv_type.SetGraph(vertices_graph);
	
	//calculate curvature
	m_curvature_1.SetGraph(vertices_graph);
	m_curvature_0.SetGraph(vertices_graph);
	const int kMaxNeighbNum = 10;
	CurvatureCalculator<VerticesGraph> curvature_calculator(kMaxNeighbNum);

	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		Vec2d curvatures;
		curvature_calculator.CalculateCurvatureCubic(vertices_graph, get(boost::vertex_info_3d, vertices_graph), *curr_vertice, curvatures);
		m_curvature_1[*curr_vertice] = curvatures[1];
		m_curvature_0[*curr_vertice] = curvatures[0];
	}

	//filter curvature
	ContPropMap<VerticesGraph, std::vector<double>, VERTEX> filtered_curvature_1(vertices_graph);
	ContPropMap<VerticesGraph, std::vector<double>, VERTEX> filtered_curvature_0(vertices_graph);
	typedef MedianKernel<double, size_t> KernelType;
	const size_t kMedRad = 1;
	FilterGraphEdgeDist(KernelType(kMedRad), vertices_graph, m_curvature_1, filtered_curvature_1);
	FilterGraphEdgeDist(KernelType(kMedRad), vertices_graph, m_curvature_0, filtered_curvature_0);
	//m_curvature_1 = filtered_curvature_1;
	//m_curvature_0 = filtered_curvature_0;
	/*typedef GaussianKernel<cv::Point3d, double> KernelType; 

	GraphDistFilter<VerticesGraph, ProxyCoordMapGraph, GaussianKernel<cv::Point3d, double>> gauss_filter_1(GaussianKernel<cv::Point3d, double>(0.5));
	ContPropMap<VerticesGraph, std::vector<double>, VERTEX> filtered_curvature_1(vertices_graph);
	ContPropMap<VerticesGraph, std::vector<double>, VERTEX> filtered_curvature_0(vertices_graph);
	gauss_filter_1.Filter(vertices_graph, map_3d_coord, curvature_1, filtered_curvature_1);
	gauss_filter_1.Filter(vertices_graph, map_3d_coord, curvature_0, filtered_curvature_0);*/

	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		if (/*filtered_*/m_curvature_1[*curr_vertice] > 0.0)
		{
			m_vertex_curv_type[*curr_vertice] = kConcaveType;
		}
		else if (/*filtered_*/m_curvature_0[*curr_vertice] < 0.0)
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
		types[m_vertex_curv_type[curr_triangle.a]]++;
		types[m_vertex_curv_type[curr_triangle.b]]++;
		types[m_vertex_curv_type[curr_triangle.c]]++;

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
void SingularPointsFinder::SegmentMolecularSurface(const size_t max_segm_size, const Mesh& mesh)
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
		m_vertex_score_map[curr_triangle.a][curr_segm]++;
		m_vertex_score_map[curr_triangle.b][curr_segm]++;
		m_vertex_score_map[curr_triangle.c][curr_segm]++;
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
	//segment remaining vertices
	
	/*const size_t segm_num_to_split = all_segm_num;
	const size_t parts_num = all_segm_num / 2 + 1;
	SplitSegment(vertices_graph, segm_num_to_split, parts_num, m_vertex_segm, all_segm_num);*/
}

void SingularPointsFinder::FindSegmentsGraphAndCenters(const Mesh& mesh)
{
	const VerticesGraph& vertices_graph = mesh.vertices;
	CalculateSegmentGraphAndCenters(vertices_graph, m_vertex_segm, m_singular_points, m_singular_points_graph);
}



void SingularPointsFinder::CalculateSingularPointsHistograms()
{
	m_sing_pts_histo.SetGraph(m_singular_points_graph);
	//CV_Assert(kHistSize == GetTypesNum());

	/*for (auto node_iter = vertices(m_singular_points_graph).first, end_iter = vertices(m_singular_points_graph).second; 
		node_iter != end_iter; ++node_iter)
	{
		std::array<uint8_t, kHistSize>& curr_hist = m_sing_pts_histo[*node_iter];
		std::fill(curr_hist.begin(), curr_hist.end(), 0);
		const size_t curr_type = m_vertices_type_map[get(boost::vertex_parent, m_singular_points_graph, *node_iter)];
		++curr_hist[curr_type];

		for (auto neighb_iter = adjacent_vertices(*node_iter, m_singular_points_graph).first,
			end_neighb_iter = adjacent_vertices(*node_iter, m_singular_points_graph).second; 
			neighb_iter != end_neighb_iter; ++neighb_iter)
		{
			const size_t neighb_type = m_vertices_type_map[get(boost::vertex_parent, m_singular_points_graph, *neighb_iter)];
			++curr_hist[neighb_type];
		}
	}*/
}

}