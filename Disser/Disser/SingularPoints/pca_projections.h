#pragma once

#include <vector>
#include "boost/graph/properties.hpp"
#include "boost/graph/graph_traits.hpp"
#include "opencv2/core/core.hpp"

namespace molecule_descriptor
{

template <typename Graph, typename CoordMap, typename DoublePropMap, typename DistVertMap, typename CVPCAPropMap>
void FindVectorsForProjection(const Graph& vertices_graph, const CoordMap& coord_map, 
	const std::vector<DoublePropMap>& functions, 
	const DistVertMap& dist_vert_map, const double dist_thresh,
	CVPCAPropMap& scale_space_projecter)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;
	
	struct DataFiller
	{
		void operator()(const CoordMap& coord, const double coord_scale, 
			const std::vector<DoublePropMap>& props, const double props_scale,
			const VertexDescriptor curr_vert, const int data_row, Mat_<double>& data) const
		{
			data(data_row, 0) = coord_scale * coord[curr_vert].x;
			data(data_row, 1) = coord_scale * coord[curr_vert].y;
			data(data_row, 2) = coord_scale * coord[curr_vert].z;

			for (int col = 3; col < 3 + ISingularPointsFinder::SURF_PROPS_NUM; ++col)
			{
				data(data_row, col) = props_scale * props[col - 3][curr_vert];
			}
		}
	};

	std::vector<VertexDescriptor> vertices_withit_dist;
	scale_space_projecter.SetGraph(vertices_graph);
		
	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		//find vertices that we need to process
		GetVerticesWithinDistPlusAdjacent(*curr_vertice, vertices_graph, dist_vert_map, dist_thresh, vertices_withit_dist);
		//fill data
		const double kCoordMult = 1000.0;
		const double kPropMult = 1.0;
		cv::Mat_<double> data(vertices_withit_dist.size(), 3 + ISingularPointsFinder::SURF_PROPS_NUM);
		size_t curr_vert_in_mask = 0;
		for (auto neighb_vertice = vertices_withit_dist.begin(), end_neighb_vertices =  vertices_withit_dist.end(); 
			neighb_vertice != end_neighb_vertices; ++neighb_vertice)
		{
			DataFiller()(coord_map, kCoordMult, functions, kPropMult, 
				*neighb_vertice, curr_vert_in_mask, data);
			++curr_vert_in_mask;
		}	
		//calculate coordiantes
		cv::PCA pca(data, cv::Mat()/*data.row(0)*/, CV_PCA_DATA_AS_ROW, 4);
		//find images of coordiante vectors basis
		cv::Mat_<double> coord_base[3];	//basis vectors for 3d coordinated
		cv::Mat_<double> proj_coord_base[3];//projections for basis vectors	
		for (int coord_ind = 0; coord_ind < 3; ++coord_ind)
		{
			coord_base[coord_ind].create(1, 3 + ISingularPointsFinder::SURF_PROPS_NUM);
			coord_base[coord_ind].setTo(0);
			coord_base[coord_ind](0, coord_ind) = 1;
			//parallel shift of vector to the origin of the new coordinate system
			coord_base[coord_ind] += pca.mean;
			pca.project(coord_base[coord_ind], proj_coord_base[coord_ind]);
			const double proj_norm = cv::norm(proj_coord_base[coord_ind]);

			if (proj_norm > 0.0)
			{
				proj_coord_base[coord_ind] *= 1.0/ proj_norm;
			}				
		}
		//find projection of properties basis vectors
		cv::Mat_<double> prop_base[ISingularPointsFinder::SURF_PROPS_NUM];	//basis vectors for 3d coordinated
		cv::Mat_<double> proj_prop_base[ISingularPointsFinder::SURF_PROPS_NUM];//projections for basis vectors	
		for (int prop_ind = 0; prop_ind <  ISingularPointsFinder::SURF_PROPS_NUM; ++prop_ind)
		{
			prop_base[prop_ind].create(1, 3 + ISingularPointsFinder::SURF_PROPS_NUM);
			//parallel shift of vector to the origin of the new coordinate system
			pca.mean.copyTo(prop_base[prop_ind]);
			prop_base[prop_ind](0, prop_ind + 3) = functions[prop_ind][*curr_vertice] * kPropMult;
			pca.project(prop_base[prop_ind], proj_prop_base[prop_ind]);
			//make projection orthogonal to each of the projected coordinate vectors
			for (int coord_ind = 0; coord_ind < 3; ++coord_ind)
			{
				const double proj_norm = proj_prop_base[prop_ind].dot(proj_coord_base[coord_ind]);
				proj_prop_base[prop_ind] -= proj_norm * proj_coord_base[coord_ind];
			}			
		}
		//calculate new property vector as sum of orthogonalized projected vectors
		cv::Mat_<double> new_prop(1, 4);
		new_prop.setTo(0);

		for (int prop_ind = 0; prop_ind <  ISingularPointsFinder::SURF_PROPS_NUM; ++prop_ind)
		{
			new_prop += proj_prop_base[prop_ind];			
		}
		
		cv::Mat_<double> vect_to_project;
		scale_space_projecter[*curr_vertice].pca = pca;
		scale_space_projecter[*curr_vertice].vect_to_project_on = new_prop;
	}
}


template <typename Graph, typename DoublePropMap, typename DistVertMap, typename CVPCAPropMap>
void SetPCAAsMeanOfFunction(const Graph& vertices_graph, const std::vector<DoublePropMap>& functions, 
							  const DistVertMap& dist_vert_map, const double dist_thresh,
							  CVPCAPropMap& scale_space_projecter)
{
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VertexDescriptor;
	struct DataFiller
	{
		void operator()(const std::vector<DoublePropMap>& props,
			const VertexDescriptor curr_vert, const int data_row, Mat_<double>& data) const
		{
			for (int col = 0; col < ISingularPointsFinder::SURF_PROPS_NUM; ++col)
			{
				data(data_row, col) = props[col][curr_vert];
			}
		}
	};
	std::vector<VertexDescriptor> vertices_withit_dist;
	scale_space_projecter.SetGraph(vertices_graph);

	for (auto curr_vertice = vertices(vertices_graph).first, end_vertices = vertices(vertices_graph).second; 
		curr_vertice != end_vertices; ++curr_vertice)
	{
		//find vertices that we need to process
		GetVerticesWithinDistPlusAdjacent(*curr_vertice, vertices_graph, dist_vert_map, dist_thresh, vertices_withit_dist);
		//fill data
		cv::Mat_<double> data(vertices_withit_dist.size(), ISingularPointsFinder::SURF_PROPS_NUM);
		size_t curr_vert_in_mask = 0;
		for (auto neighb_vertice = vertices_withit_dist.begin(), end_neighb_vertices =  vertices_withit_dist.end(); 
			neighb_vertice != end_neighb_vertices; ++neighb_vertice)
		{
			DataFiller()(functions, *neighb_vertice, curr_vert_in_mask, data);
			++curr_vert_in_mask;
		}	
		//calculate coordiantes
		cv::Mat_<double> mean(1, ISingularPointsFinder::SURF_PROPS_NUM);
		mean.setTo(0);
		cv::PCA pca(data, mean, CV_PCA_DATA_AS_ROW, 1);
		cv::Mat_<double> vect_to_project_on(1, 1);
		vect_to_project_on.setTo(1);
		scale_space_projecter[*curr_vertice].pca = pca;
		scale_space_projecter[*curr_vertice].vect_to_project_on = vect_to_project_on;
	}
}

//functions for PCA projections
//this is not for original vectors but for the difference of original vector
template <class PropMap, class CoordMap, class PCAProjecter, class VertDescr>
double ProjectPCADiffDifferentVert(const std::vector<PropMap>& prop_map, const CoordMap& coord_map, const VertDescr vert,
								   const PCAProjecter& projecter)
{
	if (projecter.pca.mean.total() == prop_map.size() + 3)
	{
		cv::Mat_<double> prop_vect(1, prop_map.size() + 3);
		prop_vect(0, 0)=coord_map[vert].x;
		prop_vect(0, 1)=coord_map[vert].y;
		prop_vect(0, 2)=coord_map[vert].z;
		for (int ind = 3; ind < 3 + prop_map.size(); ++ind)
		{
			prop_vect(0, ind) = prop_map[ind - 3][vert];
		}
		prop_vect += projecter.pca.mean;
		cv::Mat_<double> proj_prop_vect;
		projecter.pca.project(prop_vect, proj_prop_vect);
		double res = proj_prop_vect.dot(projecter.vect_to_project_on)/ cv::norm(projecter.vect_to_project_on);
		return res;
	}
	else if (projecter.pca.mean.total() == prop_map.size())
	{
		cv::Mat_<double> prop_vect(1, prop_map.size());

		for (int ind = 0; ind < prop_map.size(); ++ind)
		{
			prop_vect(0, ind) = prop_map[ind][vert];
		}
		prop_vect += projecter.pca.mean;
		cv::Mat_<double> proj_prop_vect;
		projecter.pca.project(prop_vect, proj_prop_vect);
		double res = proj_prop_vect.dot(projecter.vect_to_project_on)/ cv::norm(projecter.vect_to_project_on);
		return res;
	}
	CV_Assert(0);
}

template <class PropMap, class CoordMap, class PCAProjMap, class VertDescr>
double ProjectPCADiff(const std::vector<PropMap>& prop_map, const CoordMap& coord_map, const PCAProjMap& proj_map, const VertDescr vert)
{
	return ProjectPCADiffDifferentVert(prop_map, coord_map, vert, proj_map[vert]);
}

template <class PropMap, class CoordMap, class PCAProjecter,  class VertDescr>
double ProjectPCADifferentVert(
	const std::vector<PropMap>& prop_map, const CoordMap& coord_map, const VertDescr vert, const PCAProjecter& projecter)
{
	if (projecter.pca.mean.total() == prop_map.size() + 3)
	{
		cv::Mat_<double> prop_vect(1, prop_map.size() + 3);
		prop_vect(0, 0)=coord_map[vert].x;
		prop_vect(0, 1)=coord_map[vert].y;
		prop_vect(0, 2)=coord_map[vert].z;
		for (int ind = 3; ind < 3 + prop_map.size(); ++ind)
		{
			prop_vect(0, ind) = prop_map[ind - 3][vert];
		}
		cv::Mat_<double> proj_prop_vect;
		projecter.pca.project(prop_vect, proj_prop_vect);
		double res = proj_prop_vect.dot(projecter.vect_to_project_on)/ cv::norm(projecter.vect_to_project_on);
		return res;
	}
	else if (projecter.pca.mean.total() == prop_map.size())
	{
		cv::Mat_<double> prop_vect(1, prop_map.size());

		for (int ind = 0; ind < prop_map.size(); ++ind)
		{
			prop_vect(0, ind) = prop_map[ind][vert];
		}
		cv::Mat_<double> proj_prop_vect;
		projecter.pca.project(prop_vect, proj_prop_vect);
		double res = proj_prop_vect.dot(projecter.vect_to_project_on)/ cv::norm(projecter.vect_to_project_on);
		return res;
	}
	CV_Assert(0);
}

template <class PropMap, class CoordMap, class PCAProjMap, class VertDescr>
double ProjectPCA(const std::vector<PropMap>& prop_map, const CoordMap& coord_map, const PCAProjMap& proj_map, const VertDescr vert)
{
	return ProjectPCADifferentVert(prop_map, coord_map, vert, proj_map[vert]);
}
}