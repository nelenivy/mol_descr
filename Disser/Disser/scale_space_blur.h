#pragma once

#include <vector>
#include <memory>
#include <boost/graph/graph_traits.hpp>
#include "opencv2/core/core.hpp"
#include "GraphLib/array_property_map.h"
#include "SingularPoints/points_keeper.h"
#include "GraphLib\graph_filter.h"
#include "mesh_weighted_filter.h"
#include "GraphLib/graph_dist_calculate.h"
namespace molecule_descriptor
{

template <class VerticesGraph, class CoordMapT, class ConvolutionKernel>
class ScaleSpaceBlurrer
{
public:
	typedef typename boost::graph_traits<VerticesGraph>::vertex_descriptor vertex_descriptor;
	typedef ContPropMap<VerticesGraph, std::vector<double>, VERTEX> DoubleVertGraphProp;
	typedef ContPropMap<VerticesGraph, std::vector<cv::Point3d>, VERTEX> Point3dVertGraphProp;

	ScaleSpaceBlurrer() : m_init_sigma(0),
	m_sigma_diff(0),
	m_offset(0.0),
	m_additive(false),
	m_exp_approx(std::make_shared<ExpApproxForGauss>(1000))
	{ }

	void Init(const double init_sigma, const double sigma_diff, const double offset, const bool additive, const bool recursive, const double sigma_delta)
	{
		m_init_sigma = init_sigma;
		m_sigma_diff = sigma_diff;
		m_additive = additive;
		m_offset = offset;
		m_recursive = recursive;
		m_recursive_sigma_delta = sigma_delta;
	}

	double GetSigma(const int level)
	{
		if (m_additive)
		{
			return m_init_sigma + level * m_sigma_diff + m_offset;
		}
		else
		{
			return m_init_sigma * pow(m_sigma_diff, level) + m_offset;
		}
	}
	template <class PropMapT>
	double GetValInVertex(const VerticesGraph& graph, const PropMapT& prop_map, const CoordMapT& coord_map, const vertex_descriptor vert, const int level)
	{
		ConvolutionKernel av_kernel(GetSigma(level));
		return m_av_finder.GetAverageInRad(graph,prop_map,coord_map,av_kernel,vert);
	}

	template <class TrianglesGraph, class PropMapT, class VertVertDistMap, class VertTrDistMap>
	void MakeScaleSpaceNonRecursive(const VerticesGraph& vertices_graph, const VertVertDistMap& vert_vert_dist,
		const TrianglesGraph& tr_graph, const VertTrDistMap& vert_tr_dist,
		const std::vector<PropMapT>& prop_map_vect, 
		const int levels_num, const bool post_filter,
		std::vector<std::vector<DoubleVertGraphProp>>& filtered_prop_map_vect)
	{
	}
	template <class TrianglesGraph, class PropMapT, class VertVertDistMap, class VertTrDistMap>
	void MakeScaleSpace(const VerticesGraph& vertices_graph, const VertVertDistMap& vert_vert_dist,
		const TrianglesGraph& tr_graph, const VertTrDistMap& vert_tr_dist,
		const std::vector<PropMapT>& prop_map_vect, 
		const int levels_num, const bool post_filter,
		std::vector<std::vector<DoubleVertGraphProp>>& filtered_prop_map_vect)
	{
		filtered_prop_map_vect.assign(levels_num, std::vector<DoubleVertGraphProp>(prop_map_vect.size(), DoubleVertGraphProp(vertices_graph)));
		
		if (post_filter) 
		{
			m_filtered_before_postproc.resize(levels_num);
		}
//#pragma omp parallel for
		double true_sigma_for_recursive = 0.0;

		for (int curr_level = 0; curr_level < levels_num; ++curr_level)
		{
			if (post_filter) 
			{
				m_filtered_before_postproc[curr_level].assign(prop_map_vect.size(), DoubleVertGraphProp(vertices_graph));
			}
			const double curr_sigma = GetSigma(curr_level);
			//ConvolutionKernel av_kernel(curr_sigma);
			std::vector<DoubleVertGraphProp>& filtered_output = post_filter ? 
				m_filtered_before_postproc[curr_level] : filtered_prop_map_vect[curr_level];
			
			//GaussianKernelWeightedDistTable<double> av_kernel(curr_sigma, m_exp_approx);
			if (m_recursive)
			{
				m_recursive_prop_map = curr_level > 0 ? filtered_prop_map_vect[curr_level - 1] : prop_map_vect;
				GaussianKernelWeightedDist<double> av_kernel(m_recursive_sigma_delta);
				m_recursive_prop_map_res.assign(prop_map_vect.size(), DoubleVertGraphProp(vertices_graph));
				for (;true_sigma_for_recursive < curr_sigma; true_sigma_for_recursive += m_recursive_sigma_delta)
				{
					std::cout << true_sigma_for_recursive << "\n";
					FilterMeshWeightedFuncMultiThread(vertices_graph, tr_graph, vert_vert_dist, vert_tr_dist, m_recursive_prop_map, true, 8,
						av_kernel, m_recursive_prop_map_res);
					m_recursive_prop_map = m_recursive_prop_map_res;
				}
				true_sigma_for_recursive -= m_recursive_sigma_delta;
				filtered_output = m_recursive_prop_map_res;
			}
			else
			{
				GaussianKernelWeightedDist<double> av_kernel(curr_sigma);
				FilterMeshWeightedFuncMultiThread(vertices_graph, tr_graph, vert_vert_dist, vert_tr_dist, prop_map_vect, true, 8,
					av_kernel, filtered_output);
			}
			//FilterGraphDist(av_kernel, vertices_graph, coord_map, prop_map_vect, filtered_output);///////////

			if (post_filter)
			{
				typedef MedianKernel<double, size_t> MedianFilter;
				const size_t kMedianRadius = 1;
				MedianFilter med_filter(kMedianRadius);
				FilterGraphEdgeDist(med_filter, vertices_graph, filtered_output, filtered_prop_map_vect[curr_level]);
			}

		}			
	}

	template <class TrianglesGraph, class VertVertDistMap, class VertTrDistMap>
	void MakeScaleSpaceSpherical (const VerticesGraph& vertices_graph, const VertVertDistMap& vert_vert_dist,
		const TrianglesGraph& tr_graph, const VertTrDistMap& vert_tr_dist,
		const std::vector<Point3dVertGraphProp>& prop_map_vect, 
		const int levels_num, 
		std::vector<std::vector<Point3dVertGraphProp>>& filtered_prop_map_vect)
	{
		filtered_prop_map_vect.assign(levels_num, std::vector<Point3dVertGraphProp>(prop_map_vect.size(), Point3dVertGraphProp(vertices_graph)));
		double true_sigma_for_recursive = 0.0;

		for (int curr_level = 0; curr_level < levels_num; ++curr_level)
		{
			const double curr_sigma = GetSigma(curr_level);
			std::vector<Point3dVertGraphProp>& filtered_output = filtered_prop_map_vect[curr_level];

			//for spherical data only recursive variant of smoothing because it gives true solution of heat equation
			m_recursive_spherical_prop_map = curr_level > 0 ? filtered_prop_map_vect[curr_level - 1] : prop_map_vect;
			GaussianKernelWeightedDist<cv::Point3d> av_kernel(m_recursive_sigma_delta);
			m_recursive_spherical_prop_map_res.assign(prop_map_vect.size(), Point3dVertGraphProp(vertices_graph));
			for (;true_sigma_for_recursive < curr_sigma; true_sigma_for_recursive += m_recursive_sigma_delta)
			{
				std::cout << true_sigma_for_recursive << "\n";
				FilterMeshWeightedFuncMultiThread(vertices_graph, tr_graph, vert_vert_dist, vert_tr_dist, m_recursive_spherical_prop_map, true, 8,
					av_kernel, m_recursive_spherical_prop_map_res);
				//reproject 3d data on sphere
				for (size_t curr_prop = 0; curr_prop < m_recursive_spherical_prop_map_res.size(); ++curr_prop)
				{
					for (auto curr_vert = vertices(vertices_graph).first, end_vert = vertices(vertices_graph).second;
						curr_vert != end_vert; ++curr_vert)
					{
						m_recursive_spherical_prop_map_res[curr_prop][*curr_vert].x /= cv::norm(m_recursive_spherical_prop_map_res[curr_prop][*curr_vert]);
						m_recursive_spherical_prop_map_res[curr_prop][*curr_vert].y /= cv::norm(m_recursive_spherical_prop_map_res[curr_prop][*curr_vert]);
						m_recursive_spherical_prop_map_res[curr_prop][*curr_vert].z /= cv::norm(m_recursive_spherical_prop_map_res[curr_prop][*curr_vert]);
					}
				}
				m_recursive_spherical_prop_map = m_recursive_spherical_prop_map_res;
			}
			true_sigma_for_recursive -= m_recursive_sigma_delta;
			filtered_prop_map_vect[curr_level] = m_recursive_spherical_prop_map_res;
		}			
	}

	template <class TrianglesGraph, class PropMapT, class VertVertDistMap, class VertTrDistMap>
	void MakeScaleSpace(const VerticesGraph& vertices_graph, const VertVertDistMap& vert_vert_dist,
		const TrianglesGraph& tr_graph, const VertTrDistMap& vert_tr_dist,
		const std::vector<PropMapT>& prop_map_vect, 
		const int levels_num, const bool post_filter,
		std::vector<DoubleVertGraphProp>& filtered_prop_map)
	{
		std::vector<ProxyPropMap<const PropMapT>, IdenticalTransformFunc> prop_map_ref(1, prop_map);
		std::vector<std::vector<DoubleVertGraphProp>> filtered_vect;

		MakeScaleSpace(vertices_graph, vert_vert_dist,
			tr_graph, vert_tr_dist,prop_map_vect, 
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
	double m_offset;
	bool m_additive;//true - additive, false - multiplicative
	bool m_recursive;//true - additive, false - multiplicative
	double m_recursive_sigma_delta;
	AverageFinder<VerticesGraph, CoordMapT> m_av_finder;
	std::vector<std::vector<DoubleVertGraphProp>> m_filtered_before_postproc;
	std::vector<DoubleVertGraphProp> m_recursive_prop_map;
	std::vector<DoubleVertGraphProp> m_recursive_prop_map_res;

	std::vector<Point3dVertGraphProp> m_recursive_spherical_prop_map;
	std::vector<Point3dVertGraphProp> m_recursive_spherical_prop_map_res;

	std::shared_ptr<ExpApproxForGauss> m_exp_approx;

};

template <class VerticesGraph, class CoordMapT, class PropMapT, class PropTypesMap>
void FindScaleSingularPointsOnFunc(const VerticesGraph& vertices_graph, const CoordMapT& coord_map, std::vector<PropMapT>& scales_prop_map, 
								   const std::vector<PropTypesMap>& scales_prop_map_types, const double dist_thresh,
								   std::vector<std::vector<typename boost::graph_traits<VerticesGraph>::vertex_descriptor>>& sing_points,
								 const bool combine_from_diff_levels)
{
	CV_Assert(scales_prop_map.size() > 2);
	CV_Assert(scales_prop_map.size() == scales_prop_map_types.size());

	const size_t levels_num = scales_prop_map.size();
	typedef std::vector<typename boost::graph_traits<VerticesGraph>::vertex_descriptor> SingPtsCont;
	std::vector<SingPtsCont> sing_pts_extended(levels_num);
	
	for (size_t curr_level = 0; curr_level < levels_num; ++curr_level)
	{
		FindLocalMaxAndMin(vertices_graph, scales_prop_map[curr_level], sing_pts_extended[curr_level], 
			std::greater<double>(), std::less<double>());
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