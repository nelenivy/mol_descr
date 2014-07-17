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
	template <class PointsCont, class CurvatureTypeCont>
	void AddPoints(const PointsCont& point_cont, const CurvatureTypeCont& curv_cont)
	{
		for (size_t ind = 0; ind < point_cont.size(); ++ind)
		{
			AddPoint(point_cont[ind], curv_cont);
		}
	}
	template <class CurvatureTypeCont>
	void AddPoint(vertex_descriptor new_point, const CurvatureTypeCont& curv_cont)
	{//look if there is any point close to this
		static_assert(
			std::is_integral<typename std::decay<decltype(curv_cont[0])>::type>::value, "");
		bool is_close_to_something = false;
		const auto info_3d_map = get(boost::vertex_info_3d, m_graph.get());
		const cv::Point3d new_coord = info_3d_map[new_point].Center();

		for (size_t ind = 0; ind < m_points.size(); ++ind)
		{
			const cv::Point3d curr_coord = info_3d_map[m_points[ind]].Center();
			if (cv::norm(curr_coord - new_coord) < m_dist_thresh && curv_cont[m_points[ind]] == curv_cont[new_point])
			{
				is_close_to_something = true;
				break;
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