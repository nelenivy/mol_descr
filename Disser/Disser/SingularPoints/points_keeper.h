#pragma once

#include <vector>
#include <functional>
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/properties.hpp"
#include "boost/tuple/tuple.hpp"
#include "boost/property_map/property_map.hpp"
#include "opencv2/core/core.hpp"
#include "mesh_types.h"

namespace molecule_descriptor
{
template <typename ElemT>
struct ContainerImitator
{
	ContainerImitator(const ElemT& elem) : elem(elem) {}
	size_t size() const
	{
		return 1;
	}
	const ElemT& operator[](const size_t /*ind*/) const
	{
		return elem;
	}
	const ElemT& elem;
};
//if new point lies too close to existing it isn't added
template <class GraphT>
class PointsKeeper
{
public:
	typedef typename boost::graph_traits<GraphT>::vertex_descriptor vertex_descriptor;
	PointsKeeper(const double dist_thresh, const GraphT& graph)
		: m_dist_thresh(dist_thresh), m_graph(graph) { }
	void Clear()
	{
		m_points.clear();
	}
	template <class PointsCont, class ContOfTypeCont>
	void AddSeveralPointsSeveralTypes(const std::vector<PointsCont>& point_cont, const ContOfTypeCont& types_cont)
	{
		for (size_t ind = 0; ind < point_cont.size(); ++ind)
		{
			AddPointsSeveralTypes(point_cont[ind], types_cont);
		}
	}
	template <class PointsCont, class ContOfTypeCont>
	void AddPointsSeveralTypes(const PointsCont& point_cont, const ContOfTypeCont& types_cont)
	{
		for (size_t ind = 0; ind < point_cont.size(); ++ind)
		{
			AddPointSeveralTypes(point_cont[ind], types_cont);
		}
	}
	template <class PointsCont, class CurvatureTypeCont>
	void AddPoints(const PointsCont& point_cont, const CurvatureTypeCont& types)
	{
		for (size_t ind = 0; ind < point_cont.size(); ++ind)
		{
			AddPoint(point_cont[ind], types);
		}
	}
	template <class CurvatureTypeCont>
	void AddPoint(vertex_descriptor new_point, const CurvatureTypeCont& types)
	{//look if there is any point close to this
		AddPointSeveralTypes(new_point, ContainerImitator<CurvatureTypeCont>(types));
	}
	template <class ContOfTypeCont>
	void AddPointSeveralTypes(vertex_descriptor new_point, const ContOfTypeCont& types_cont)
	{//look if there is any point close to this
		static_assert(
			std::is_integral<typename std::decay<decltype(types_cont[0][0])>::type>::value, "");
		bool is_close_to_something = false;
		const auto info_3d_map = get(boost::vertex_info_3d, m_graph.get());
		const cv::Point3d new_coord = info_3d_map[new_point].Center();

		for (size_t ind = 0; ind < m_points.size(); ++ind)
		{
			const cv::Point3d curr_coord = info_3d_map[m_points[ind]].Center();
			//check distance to other points
			if (cv::norm(curr_coord - new_coord) < m_dist_thresh)
			{
				bool same_type = true;
				//if it's close check points types
				for (size_t type_ind = 0; type_ind < types_cont.size(); ++type_ind)
				{
					if (types_cont[type_ind][m_points[ind]] != types_cont[type_ind][new_point])
					{
						same_type = false;
						break;
					}
				}

				if (same_type)
				{
					is_close_to_something = true;
					break;
				}
			}
		}

		if (!is_close_to_something)
		{
			m_points.push_back(new_point);
		}
	}
	const std::vector<vertex_descriptor>& Points()
	{
		return m_points;
	}
private:
	std::vector<vertex_descriptor> m_points;
	double m_dist_thresh;
	std::reference_wrapper<const GraphT> m_graph;
};
}