#include "singular_points_finder.h"

#include <algorithm>

#include "GraphLib\connected_components_segmentator.h"
#include "GraphLib\segmentation_types.h"
#include "CommonUtilities/attributes_container.h"
#include "GraphLib\graph_operations.h"

#include "GraphLib\curvature_calculator.h"

#include "singular_point_types.h"
#include "points_marker.h"

namespace molecule_descriptor
{

using namespace std;
using namespace cv;

std::shared_ptr<ISingularPointsFinder> CreateSingularPointsFinder()
{
	ISingularPointsFinder* instance = new SingularPointsFinder;
	return std::shared_ptr<ISingularPointsFinder>(instance, ReleaseDeleter());
}

void SingularPointsFinder::Process(const std::vector<cv::Point3d>& vertices, const std::vector<cv::Point3d>& normals, 
	const std::vector<cv::Point3i>& triangles, const std::vector<std::pair<cv::Point3d, double>>& charges)
{
	m_mesh_keeper.ConstructMesh(vertices, normals, triangles);
	CalculateVerticesSurfaceType();
	const int kMaxSegmSize = 500;
	SegmentMolecularSurface(kMaxSegmSize);
	FindSegmentsGraphAndCenters();
	CalculatePotential(m_singular_points, charges);
	CalculateSingularPointsTypes();
	CalculateSingularPointsHistograms();
}

size_t SingularPointsFinder::GetTypesNum()
{
	return GetRange<ElectricSign>() * GetRange<SurfaceTypeProp>();
}

void SingularPointsFinder::GetMarkedSingularPoints(std::vector<MarkedSingularPoint>& marked_singular_points)
{
	marked_singular_points.reserve(m_singular_points.size());

	for (size_t ind = 0; ind < m_singular_points.size(); ++ind)
	{
		const size_t curr_type = m_singular_points[ind].attr.Get<SingularPointType>();
		const cv::Point3d& curr_coord = m_singular_points[ind].Center();
		marked_singular_points.emplace_back(curr_coord, curr_type);
	}
}

void SingularPointsFinder::GetNonMarkedSingularPoints(std::vector<NonMarkedSingularPoint>& non_marked_singular_points)
{
	non_marked_singular_points.resize(m_singular_points.size());

	for (size_t ind = 0; ind < m_singular_points.size(); ++ind)
	{
		auto& curr_point = non_marked_singular_points[ind];
		curr_point.Property().SurfaceType() = m_singular_points[ind].attr.Get<SurfaceTypeProp>();
		curr_point.Property().Charge() = m_singular_points[ind].attr.Get<ElectricSign>();
		curr_point.Property().ElectricPotential() = m_singular_points[ind].attr.Get<ElectricPotential>();
		curr_point.Coord() = m_singular_points[ind].Center();
	}
}

void SingularPointsFinder::GetSingularPointsHisto(std::vector<HistogramSingularPoint<kHistSize>>& singular_points_hist)
{
	singular_points_hist.resize(m_sing_pts_histo.size());

	for (size_t ind = 0; ind < m_sing_pts_histo.size(); ++ind)
	{
		auto& curr_point = singular_points_hist[ind];
		curr_point.Property() = m_sing_pts_histo[ind];
		curr_point.Coord() = m_singular_points[ind].Center();
	}
}

void SingularPointsFinder::GetSegmentedVertices(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_segm_numbers)
{
	const std::vector<MeshVertice>& vertices = *(m_mesh_keeper.GetMeshVertices());
	vertices_with_segm_numbers.resize(vertices.size());

	for (size_t ind = 0; ind < vertices.size(); ++ind)
	{
		const size_t curr_segm = vertices[ind].attr.Get<SegmentNumProp>();
		const cv::Point3d curr_coord = vertices[ind].element->Center();
		vertices_with_segm_numbers[ind] = std::make_pair(curr_coord, curr_segm);
	}
}

void SingularPointsFinder::GetVerticesWithTypes(std::vector<std::pair<cv::Point3d, size_t>>& vertices_with_types)
{
	const std::vector<MeshVertice>& vertices = *(m_mesh_keeper.GetMeshVertices());
	vertices_with_types.resize(vertices.size());

	for (size_t ind = 0; ind < vertices.size(); ++ind)
	{
		const size_t curr_segm = vertices[ind].attr.Get<SurfaceTypeProp>();
		const cv::Point3d curr_coord = vertices[ind].element->Center();
		vertices_with_types[ind] = std::make_pair(curr_coord, curr_segm);
	}
}

void SingularPointsFinder::Clear()
{
	m_singular_points.clear();
	m_singular_points_graph.clear();
	m_sing_pts_histo.clear();
}

void SingularPointsFinder::CalculateVerticesSurfaceType()
{
	const std::vector<GraphNode<Vertice>>& vertices = *(m_mesh_keeper.GetMeshVertices());
	const std::vector<GraphNode<Triangle>>& triangles = *(m_mesh_keeper.GetMeshTriangles());
	const int kMaxNeighbNum = 10;
	CurvatureCalculator<Vertice> curvature_calculator(kMaxNeighbNum);

	for (auto vert_iter = vertices.begin(); vert_iter != vertices.end(); ++vert_iter)
	{
		Vec2d curvatures;
		curvature_calculator.CalculateCurvatureCubic(*vert_iter, curvatures);
		vert_iter->element->attr.Add<SurfaceTypeProp>();
		SurfaceTypeProp::value_type surface_type;

		if (curvatures[1] > 0.0)
		{
			surface_type = kConvexType;
		}
		else if (curvatures[0] < 0.0)
		{
			surface_type = kConcaveType;
		}
		else
		{
			surface_type = kSaddleType;
		}

		vert_iter->element->attr.Get<SurfaceTypeProp>() = surface_type;
	}

	//filter surface types labels
	for (auto vert_iter = vertices.begin(); vert_iter != vertices.end(); ++vert_iter)
	{
		const size_t curr_type = vert_iter->element->attr.Get<SurfaceTypeProp>();
		size_t differences_num = 0;

		for (auto neighb_iter = vert_iter->neighbours.begin(); neighb_iter != vert_iter->neighbours.end(); ++neighb_iter)
		{
			const size_t neighb_type = (*neighb_iter)->element->attr.Get<SurfaceTypeProp>();
			if (neighb_type != curr_type)
			{
				differences_num++;
			}
		}

		if (differences_num == vert_iter->neighbours.size())
		{
			auto first_neighb = *(vert_iter->neighbours.begin());
			vert_iter->element->attr.Get<SurfaceTypeProp>() = first_neighb->element->attr.Get<SurfaceTypeProp>();
		}

		CopyPropFromElemToNode<SurfaceTypeProp>()(*vert_iter);
	}

	//assign triangles labels
	for (auto tr_iter = triangles.begin(); tr_iter != triangles.end(); ++tr_iter)
	{
		const Triangle& curr_triangle = *(tr_iter->element);
		curr_triangle.attr.Add<SurfaceTypeProp>();

		if (curr_triangle.a->attr.Get<SurfaceTypeProp>() == curr_triangle.b->attr.Get<SurfaceTypeProp>() &&
			curr_triangle.a->attr.Get<SurfaceTypeProp>() == curr_triangle.c->attr.Get<SurfaceTypeProp>())
		{
			curr_triangle.attr.Get<SurfaceTypeProp>() = curr_triangle.a->attr.Get<SurfaceTypeProp>();
		}
		else
		{
			curr_triangle.attr.Get<SurfaceTypeProp>() = kSaddleType;
		}

		CopyPropFromElemToNode<SurfaceTypeProp>()(*tr_iter);
	}
}

const size_t kMinElementsInSegm = 5;

void SingularPointsFinder::SegmentMolecularSurface(const size_t max_segm_size)
{	
	const std::vector<GraphNode<Vertice>>& vertices = *(m_mesh_keeper.GetMeshVertices());
	const std::vector<GraphNode<Triangle>>& triangles = *(m_mesh_keeper.GetMeshTriangles());	

	//segment with connected components algorithm
	const unsigned char types[] = {kConvexType, kConcaveType};
	size_t all_segm_num = 0;
	ConnectedComponentsSegmentator<Triangle> conn_comp_segmentator;

	for (size_t type_ind = 0; type_ind < 2; type_ind++)
	{
		const unsigned char curr_type = types[type_ind];
		size_t curr_segm_num = 0;
		conn_comp_segmentator.SegmentImageValue<SurfaceTypeProp>(triangles, curr_type, kMinElementsInSegm, curr_segm_num);
		//remap segments number
		vector<size_t> segments_sizes(curr_segm_num + 1, 0);
		for (auto tr_iter = triangles.begin(); tr_iter != triangles.end(); ++tr_iter)
		{
			if (tr_iter->element->attr.Get<SurfaceTypeProp>() == curr_type)
			{
				size_t& curr_segm = tr_iter->attr.Get<SegmentNumProp>();
				segments_sizes[curr_segm]++;
				curr_segm += all_segm_num;
			}
		}
		//split too big segments
		const size_t segm_num_before_split = curr_segm_num;

		for (size_t curr_segm = 1; curr_segm <= segm_num_before_split; ++curr_segm)
		{
			if (segments_sizes[curr_segm] > max_segm_size)
			{
				const size_t new_segm_num = static_cast<int>(ceil(segments_sizes[curr_segm] / static_cast<double>(max_segm_size)));
				SplitSegment(triangles, curr_segm + all_segm_num, new_segm_num, curr_segm_num);
			}
		}

		all_segm_num += curr_segm_num;
	}

	std::for_each(triangles.begin(), triangles.end(), CopyPropFromNodeToElem<SegmentNumProp>());

	//write segments_num to vertices
	all_segm_num++;
	for (auto tr_iter = triangles.begin(); tr_iter != triangles.end(); ++tr_iter)
	{
		const Triangle& curr_triangle = *(tr_iter->element);
		
		if (tr_iter->attr.Get<SegmentNumProp>() == 0)
		{
			tr_iter->attr.Get<SegmentNumProp>() = all_segm_num;//mark for segmenting with k-means
			CopyPropFromNodeToElem<SegmentNumProp>()(*tr_iter);
		}

		curr_triangle.a->attr.Add<SegmentNumProp>();
		curr_triangle.b->attr.Add<SegmentNumProp>();
		curr_triangle.c->attr.Add<SegmentNumProp>();
		curr_triangle.a->attr.Get<SegmentNumProp>() = curr_triangle.b->attr.Get<SegmentNumProp>() =
			curr_triangle.c->attr.Get<SegmentNumProp>() = curr_triangle.attr.Get<SegmentNumProp>();		
	}

	std::for_each(vertices.begin(), vertices.end(), CopyPropFromElemToNode<SegmentNumProp>());
	//segment remaining vertices
	
	const size_t segm_num_to_split = all_segm_num;
	const size_t parts_num = all_segm_num / 2 + 1;
	SplitSegment(vertices, segm_num_to_split, parts_num, all_segm_num);
}

void SingularPointsFinder::FindSegmentsGraphAndCenters()
{
	const std::vector<GraphNode<Vertice>>& vertices_graph = *(m_mesh_keeper.GetMeshVertices());
	CalculateSegmentGraphAndCenters(vertices_graph, m_singular_points, m_singular_points_graph);
}

void SingularPointsFinder::CalculateSingularPointsTypes()
{
	for (auto iter = m_singular_points.begin(); iter != m_singular_points.end(); ++iter)
	{
		const size_t curr_type = CalculateSingularPointType<SurfaceTypeProp, ElectricSign>(iter->attr);
		iter->attr.Add<SingularPointType>();
		iter->attr.Get<SingularPointType>() = curr_type;
	}
}

void SingularPointsFinder::CalculateSingularPointsHistograms()
{
	m_sing_pts_histo.resize(m_singular_points_graph.size());
	CV_Assert(kHistSize == GetTypesNum());

	for (auto node_iter = m_singular_points_graph.begin(); node_iter != m_singular_points_graph.end(); ++node_iter)
	{
		std::array<uint8_t, kHistSize>& curr_hist = m_sing_pts_histo[node_iter - m_singular_points_graph.begin()];
		std::fill(curr_hist.begin(), curr_hist.end(), 0);
		const size_t curr_type = node_iter->element->attr.Get<SingularPointType>();
		++curr_hist[curr_type];

		for (auto neighb_iter = node_iter->neighbours.begin(); neighb_iter != node_iter->neighbours.end(); ++neighb_iter)
		{
			const size_t neighb_type = (*neighb_iter)->element->attr.Get<SingularPointType>();
			++curr_hist[neighb_type];
		}
	}
}

}