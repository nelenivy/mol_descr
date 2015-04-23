#pragma once

#include <vector>
#include <boost/graph/graph_traits.hpp>
#include "opencv2/core/core.hpp"
#include "GraphLib/array_property_map.h"
#include "SingularPoints/points_keeper.h"
#include "GraphLib\graph_filter.h"

namespace molecule_descriptor
{

template <class GraphT, class CoordMapT, class ConvolutionKernel>
class ScaleSpaceBlurrer
{
public:
	typedef typename boost::graph_traits<GraphT>::vertex_descriptor vertex_descriptor;
	typedef ContPropMap<GraphT, std::vector<double>, VERTEX> DoubleVertGraphProp;

	ScaleSpaceBlurrer() : m_init_sigma(0),
	m_sigma_diff(0),
	m_additive(false)
	{ }

	void Init(const double init_sigma, const double sigma_diff, const bool additive)
	{
		m_init_sigma = init_sigma;
		m_sigma_diff = sigma_diff;
		m_additive = additive;
	}

	double GetSigma(const int level)
	{
		if (m_additive)
		{
			return m_init_sigma + level * m_sigma_diff;
		}
		else
		{
			return m_init_sigma * pow(m_sigma_diff, level);
		}
	}
	template <class PropMapT>
	double GetValInVertex(const GraphT& graph, const PropMapT& prop_map, const CoordMapT& coord_map, const vertex_descriptor vert, const int level)
	{
		ConvolutionKernel av_kernel(GetSigma(level));
		return m_av_finder.GetAverageInRad(graph,prop_map,coord_map,av_kernel,vert);
	}

	template <class CoordMapT, class PropMapT>
	void MakeScaleSpace(const GraphT& vertices_graph, const CoordMapT& coord_map, const std::vector<PropMapT>& prop_map_vect, 
		const int levels_num, const bool post_filter,
		std::vector<std::vector<DoubleVertGraphProp>>& filtered_prop_map_vect)
	{
		//calculate curvature
		filtered_prop_map_vect.assign(levels_num, std::vector<DoubleVertGraphProp>(prop_map_vect.size(), DoubleVertGraphProp(vertices_graph)));
		
		if (post_filter) 
		{
			m_filtered_before_postproc.assign(prop_map_vect.size(), DoubleVertGraphProp(vertices_graph));
		}

		for (size_t curr_level = 0; curr_level < levels_num; ++curr_level)
		{
			const double curr_sigma = GetSigma(curr_level);
			ConvolutionKernel av_kernel(curr_sigma);

			std::vector<DoubleVertGraphProp>& filtered_output = post_filter ? m_filtered_before_postproc : filtered_prop_map_vect[curr_level];
			FilterGraphDist(av_kernel, vertices_graph, coord_map, prop_map_vect, filtered_output);///////////

			if (post_filter)
			{
				typedef MedianKernel<double, size_t> MedianFilter;
				const size_t kMedianRadius = 1;
				MedianFilter med_filter(kMedianRadius);
				FilterGraphEdgeDist(med_filter, vertices_graph, filtered_output, filtered_prop_map_vect[curr_level]);
			}
		}			
	}
	template <class CoordMapT, class PropMapT>
	void MakeScaleSpace(const GraphT& vertices_graph, const CoordMapT& coord_map, const PropMapT& prop_map, 
		const int levels_num, const bool post_filter,
		std::vector<DoubleVertGraphProp>& filtered_prop_map)
	{
		std::vector<ProxyPropMap<const PropMapT>, IdenticalTransformFunc> prop_map_ref(1, prop_map);
		std::vector<std::vector<DoubleVertGraphProp>> filtered_vect;

		MakeScaleSpace(vertices_graph, coord_map, prop_map_ref, 
			levels_num, post_filter,filtered_vect);

		filtered_prop_map.resize(filtered_vect.size());

		for (size_t curr_level = 0; curr_level < filtered_vect.size(); ++curr_level)
		{
			filtered_prop_map[curr_level] = filtered_vect[curr_level][0];	
		}
	}
private:
	double m_init_sigma;
	double m_sigma_diff;
	bool m_additive;//true - additive, false - multiplicative
	AverageFinder<GraphT, CoordMapT> m_av_finder;
	std::vector<DoubleVertGraphProp> m_filtered_before_postproc;
};

template <class GraphT, class CoordMapT, class PropMapT, class PropTypesMap>
void FindScaleSingularPointsOnFunc(const GraphT& vertices_graph, const CoordMapT& coord_map, std::vector<PropMapT>& scales_prop_map, 
								   const std::vector<PropTypesMap>& scales_prop_map_types, const double dist_thresh,
								   std::vector<std::vector<typename boost::graph_traits<GraphT>::vertex_descriptor>>& sing_points,
								   const bool maximum_in_scale_space, const bool combine_from_diff_levels)
{
	CV_Assert(scales_prop_map.size() > 2);
	CV_Assert(scales_prop_map.size() == scales_prop_map_types.size());

	const size_t levels_num = scales_prop_map.size();
	typedef std::vector<typename boost::graph_traits<GraphT>::vertex_descriptor> SingPtsCont;
	std::vector<SingPtsCont> sing_pts_extended(levels_num);
	if (maximum_in_scale_space)
	{
		FindLocalMaximumsOnLevels(vertices_graph, scales_prop_map, sing_pts_extended, std::greater<double>());
	}
	else
	{
		for (size_t curr_level = 0; curr_level < levels_num; ++curr_level)
		{
			FindLocalMaximumsOfAbsVal(vertices_graph, scales_prop_map[curr_level], sing_pts_extended[curr_level], 
				std::greater<double>(), std::less<double>());
		}	
	}

	if (combine_from_diff_levels)
	{
		sing_points.clear();
		for (size_t curr_level = 1; curr_level < levels_num - 1; ++curr_level)
		{
			PointsKeeper<VerticesGraph> points_keeper_level(dist_thresh, vertices_graph);
			points_keeper_level.AddPoints(sing_pts_extended[curr_level ], scales_prop_map_types[curr_level]);
			points_keeper_level.AddPoints(sing_pts_extended[curr_level + 1], scales_prop_map_types[curr_level]);
			points_keeper_level.AddPoints(sing_pts_extended[curr_level - 1], scales_prop_map_types[curr_level]);
			sing_points.push_back(points_keeper_level.Points());
			//std::cout << points_keeper_level.Points().size() << "\n";
		}
	}
	else
	{
		sing_points = sing_pts_extended;
	}
}

}