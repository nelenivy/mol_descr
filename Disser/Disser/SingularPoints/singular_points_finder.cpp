#include "singular_points_finder.h"

#include <algorithm>

#include "SingularPoints/mesh_types.h"
#include "GraphLib\connected_components_segmentator.h"
#include "GraphLib\graph_operations.h"

#include "GraphLib\curvature_calculator.h"

#include "singular_point_types.h"

namespace molecule_descriptor
{

using namespace std;
using namespace cv;

std::shared_ptr<ISingularPointsFinder> CreateSingularPointsFinder()
{
	ISingularPointsFinder* instance = new SingularPointsFinder;
	return std::shared_ptr<ISingularPointsFinder>(instance, ReleaseDeleter());
}

double CalculatePotential(const cv::Point3d& point, const std::vector<std::pair<cv::Point3d, double>>& charges)
{
	double potential = 0.0;

	for (auto charges_iter = charges.begin(); charges_iter != charges.end(); ++charges_iter)
	{
		const double curr_dist = cv::norm(point - charges_iter->first) + 0.00000001;
		potential += charges_iter->second / curr_dist;
	}

	return potential;
}

void SingularPointsFinder::CalculateAllPotentials(const std::vector<std::pair<cv::Point3d, double>>& charges)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
	m_vertex_charge_map.SetGraph(vertices_graph);

	for (auto vertex_iter = vertices(vertices_graph).first, 
		end_iter = vertices(vertices_graph).second; vertex_iter != end_iter; ++vertex_iter)
	{
		const Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *vertex_iter).Center();
		m_vertex_charge_map[*vertex_iter] = CalculatePotential(curr_coord, charges);
	}
}

void SingularPointsFinder::Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
	const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges)
{
	m_mesh_keeper.ConstructMesh(vertices, normals, triangles);
	CalculateVerticesSurfaceType();
	const int kMaxSegmSize = 500;
	SegmentMolecularSurface(kMaxSegmSize);
	FindSegmentsGraphAndCenters();
	CalculateAllPotentials(charges);
	CalculateSingularPointsTypes();
	CalculateSingularPointsHistograms();
}

const size_t kElectricSignMax = 2;
const size_t kSurfaceTypeMax = kSaddleType;

size_t SingularPointsFinder::GetTypesNum()
{
	return (kElectricSignMax + 1) * kSurfaceTypeMax;
}

void SingularPointsFinder::GetMarkedSingularPoints(std::vector<MarkedSingularPoint>& marked_singular_points)
{
	marked_singular_points.clear();
	marked_singular_points.reserve(m_singular_points.size());
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();

	for (size_t ind = 0; ind < m_singular_points.size(); ++ind)
	{
		const auto curr_vertice = m_singular_points[ind];
		const size_t curr_type = m_vertices_type_map[curr_vertice];
		const cv::Point3d& curr_coord = get(boost::vertex_info_3d, vertices_graph, curr_vertice).Center();
		marked_singular_points.emplace_back(curr_coord, curr_type);
	}
}

void SingularPointsFinder::GetNonMarkedSingularPoints(std::vector<NonMarkedSingularPoint>& non_marked_singular_points)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
	non_marked_singular_points.resize(m_singular_points.size());

	for (size_t ind = 0; ind < m_singular_points.size(); ++ind)
	{
		auto& curr_point = non_marked_singular_points[ind];
		const auto curr_vertice = m_singular_points[ind];
		curr_point.Property().SurfaceType() = m_vertex_curv_type[curr_vertice];
		curr_point.Property().Charge() = Sign(m_vertex_charge_map[curr_vertice]);
		curr_point.Property().ElectricPotential() = m_vertex_charge_map[curr_vertice];
		curr_point.Coord() = get(boost::vertex_info_3d, vertices_graph, curr_vertice).Center();
	}
}

void SingularPointsFinder::GetSingularPointsHisto(std::vector<HistogramSingularPoint<kHistSize>>& singular_points_hist)
{
	singular_points_hist.resize(num_vertices(m_singular_points_graph));
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
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
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
	vertices_with_segm_numbers.clear();
	vertices_with_segm_numbers.reserve(num_vertices(vertices_graph));

	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		const size_t curr_segm = m_vertex_segm[*curr_vertice];
		const cv::Point3d curr_coord = get(boost::vertex_info_3d, vertices_graph, *curr_vertice).Center();
		vertices_with_segm_numbers.push_back(std::make_pair(curr_coord, curr_segm));
	}
}

void SingularPointsFinder::GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types)
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
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

	m_vertices_type_map.Clear();
	m_vertex_charge_map.Clear();
	m_vertex_curv_type.Clear();
	m_vertex_segm.Clear();
	m_triangle_curv_type.Clear();
	m_triangles_segm.Clear();
}

void SingularPointsFinder::CalculateVerticesSurfaceType()
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
	m_vertex_curv_type.SetGraph(vertices_graph);
	const int kMaxNeighbNum = 10;
	CurvatureCalculator<VerticesGraph> curvature_calculator(kMaxNeighbNum);
	
	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		Vec2d curvatures;
		curvature_calculator.CalculateCurvatureCubic(vertices_graph, *curr_vertice, curvatures);
		
		if (curvatures[1] > 0.0)
		{
			m_vertex_curv_type[*curr_vertice] = kConcaveType;
		}
		else if (curvatures[0] < 0.0)
		{
			m_vertex_curv_type[*curr_vertice] = kSaddleType;
		}
		else
		{
			m_vertex_curv_type[*curr_vertice] = kConvexType ;
		}
	}

	//filter surface types labels
	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		const size_t curr_type = m_vertex_curv_type[*curr_vertice];
		size_t differences_num = 0;

		for (auto curr_neighbour = adjacent_vertices(*curr_vertice, vertices_graph).first,
			end_neighb = adjacent_vertices(*curr_vertice, vertices_graph).second;
			curr_neighbour != end_neighb; ++curr_neighbour)
		{
			const size_t neighb_type = m_vertex_curv_type[*curr_neighbour];
			if (neighb_type != curr_type)
			{
				differences_num++;
			}
		}

		if (differences_num == in_degree(*curr_vertice, vertices_graph))
		{
			auto first_neighbour = adjacent_vertices(*curr_vertice, vertices_graph).first;
			m_vertex_curv_type[*curr_vertice] = m_vertex_curv_type[*first_neighbour];
		}
	}

	//assign triangles labels
	const TrianglesGraph& triangles_graph = m_mesh_keeper.GetMeshTriangles();
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

void SingularPointsFinder::SegmentMolecularSurface(const size_t max_segm_size)
{	
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
	const TrianglesGraph& triangles_graph = m_mesh_keeper.GetMeshTriangles();
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
		const size_t segm_num_before_split = curr_segm_num;
		curr_segm_num += all_segm_num;

		for (size_t curr_segm = 1; curr_segm <= segm_num_before_split; ++curr_segm)
		{
			if (segments_sizes[curr_segm] > max_segm_size)
			{
				const size_t parts_of_segm_num = static_cast<int>(ceil(segments_sizes[curr_segm] / static_cast<double>(max_segm_size)));
				SplitSegment(triangles_graph, curr_segm + all_segm_num, parts_of_segm_num, m_triangles_segm, curr_segm_num);
			}
		}

		all_segm_num = curr_segm_num;
	}

	//write segments_num to vertices
	all_segm_num++;
	m_vertex_score_map.SetGraph(vertices_graph);
	for (auto tr_iter = vertices(triangles_graph).first, end_iter = vertices(triangles_graph).second; 
		tr_iter != end_iter; ++tr_iter)
	{
		
		if (m_triangle_curv_type[*tr_iter] == kSaddleType)
		{
			m_triangles_segm[*tr_iter] = all_segm_num;//mark for segmenting with k-means
		}

		const MeshTriangle curr_triangle = get(boost::vertex_info_3d, triangles_graph, *tr_iter);
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
	
	const size_t segm_num_to_split = all_segm_num;
	const size_t parts_num = all_segm_num / 2 + 1;
	SplitSegment(vertices_graph, segm_num_to_split, parts_num, m_vertex_segm, all_segm_num);
}

void SingularPointsFinder::FindSegmentsGraphAndCenters()
{
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
	CalculateSegmentGraphAndCenters(vertices_graph, m_vertex_segm, m_singular_points, m_singular_points_graph);
}

struct TypeWithMax
{
	size_t type;
	IntRange range;

	int Order() const
	{
		return (type - range.min_val) / range.step;
	}
};

size_t CalculateSingularPointType(const std::vector<TypeWithMax>& type_with_max)
{
	size_t type = 0;

	for (size_t curr_type = 0; curr_type < type_with_max.size(); ++curr_type)
	{
		if (curr_type > 0)
		{
			type *= type_with_max[curr_type - 1].range.ElemsInRange();
		}

		type += type_with_max[curr_type].Order();
	}

	return type;
}

void SingularPointsFinder::CalculateSingularPointsTypes()
{
	IntRange curvature_range;
	curvature_range.min_val = kConvexType;
	curvature_range.max_val = kSaddleType;
	curvature_range.step = 1;

	IntRange charge_range;
	charge_range.min_val = -1;
	charge_range.max_val = 1;
	charge_range.step = 1; 

	std::vector<TypeWithMax> type_with_max(2);
	enum {SIGN, CURVATURE};
	type_with_max[SIGN].range = charge_range;
	type_with_max[CURVATURE].range = curvature_range;
	
	const VerticesGraph& vertices_graph = m_mesh_keeper.GetMeshVertices();
	m_vertices_type_map.SetGraph(vertices_graph);
	for (auto iter = vertices(vertices_graph).first, 
		end = vertices(vertices_graph).second; iter != end; ++iter)
	{
		const size_t curr_sign = Sign(m_vertex_charge_map[*iter]);
		type_with_max[SIGN].type = curr_sign;
		type_with_max[CURVATURE].type = m_vertex_curv_type[*iter];
		const size_t curr_type = CalculateSingularPointType(type_with_max);
		m_vertices_type_map[*iter] = curr_type;
	}
}

void SingularPointsFinder::CalculateSingularPointsHistograms()
{
	m_sing_pts_histo.SetGraph(m_singular_points_graph);
	CV_Assert(kHistSize == GetTypesNum());

	for (auto node_iter = vertices(m_singular_points_graph).first, end_iter = vertices(m_singular_points_graph).second; 
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
	}
}

}