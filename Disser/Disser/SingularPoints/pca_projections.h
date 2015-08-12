#pragma once

#include <vector>
#include "boost/graph/properties.hpp"
#include "boost/graph/graph_traits.hpp"
#include "opencv2/core/core.hpp"

namespace molecule_descriptor
{

template <typename Graph, typename CoordMap, typename DoublePropMap, typename DistVertMap, typename CVPCAPropMap>
void FindVectorsForProjection(const Graph& vertices_graph, const CoordMap& coord_map, 
	const std::vector<DoublePropMap>& scale_space, 
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
			DataFiller()(coord_map, kCoordMult, scale_space, kPropMult, 
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
			prop_base[prop_ind](0, prop_ind + 3) = scale_space[prop_ind][*curr_vertice] * kPropMult;
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
		/*for (int prop_coord_ind = 0; prop_coord_ind < 4; ++prop_coord_ind)
		{
			std::cout << proj_prop_base[0](0, prop_coord_ind) << " ";
		}
		std::cout << std::endl;*/
		for (int prop_ind = 0; prop_ind <  ISingularPointsFinder::SURF_PROPS_NUM; ++prop_ind)
		{
			new_prop += proj_prop_base[prop_ind];
			for (int prop_coord_ind = 0; prop_coord_ind < 4; ++prop_coord_ind)
			{
				//std::cout << proj_prop_base[prop_ind](0, prop_coord_ind) / proj_prop_base[0](0, prop_coord_ind) << " ";
			}
			//std::cout << std::endl;
		}
		/*cv::Mat_<double> data_basis(1, 3 + ISingularPointsFinder::PROPS_NUM), proj_basis;
		DataFiller()(coord_map, 1.0, scale_space[curr_level], 0.0, *curr_vertice, 0, data_basis);
		cv::Mat_<double> curr_data(1, 3 + ISingularPointsFinder::PROPS_NUM), proj_curr;
		DataFiller()(coord_map, 1.0, scale_space[curr_level], 1.0, *curr_vertice, 0, curr_data);
		pca.project(data_basis, proj_basis);
		pca.project(curr_data, proj_curr);
		cv::Mat_<double> new_prop = proj_curr - proj_basis;
		*/
		/*std::cout << pca.eigenvectors << "\n" << pca.mean << "\n" << new_prop << "\n";
		for (int coord_ind = 0; coord_ind < 3; ++coord_ind)
		{
			std::cout << proj_coord_base[coord_ind] << "\n";
		}	
		std::cout << "\n";*/

		cv::Mat_<double> vect_to_project;
		scale_space_projecter[*curr_vertice].pca = pca;
		scale_space_projecter[*curr_vertice].vect_to_project_on = new_prop;
		//pca.backProject(new_prop, vect_to_project);
		//vect_to_project -= pca.mean;

		//for (int prop = 0; prop < ISingularPointsFinder::PROPS_NUM; ++prop)
		//{
		//	scale_space_vect_to_project[curr_level][prop][*curr_vertice] = vect_to_project(0, prop + 3);
		//}

		////check
		//cv::Mat_<double> vect_data(1, 3 + ISingularPointsFinder::PROPS_NUM);
		//for (auto neighb_it = adjacent_vertices(*curr_vertice, vertices_graph).first,
		//	end_neighb = adjacent_vertices(*curr_vertice, vertices_graph).second; neighb_it != end_neighb; ++neighb_it)
		//{
		//	DataFiller()(coord_map, kCoordMult, scale_space[curr_level], kPropMult, *neighb_it, 0, vect_data);
		//	cv::Mat_<double> vect_projected_data;
		//	pca.project(vect_data, vect_projected_data);
		//	double prod1 = 0;
		//	double norm1 = 0;
		//	for (int coord = 0; coord < 4; ++coord)
		//	{
		//		prod1 += vect_projected_data(0, coord) * new_prop(0, coord);
		//		norm1 += new_prop(0, coord) * new_prop(0, coord);
		//	}
		//	prod1 /= norm1;
		//	vect_data -= pca.mean;
		//	double prod = 0;
		//	double norm = 0;
		//	for (int coord = 0; coord < ISingularPointsFinder::PROPS_NUM; ++coord)
		//	{
		//		prod += vect_data(0, coord + 3) * vect_to_project(0, coord + 3);
		//		norm += vect_to_project(0, coord + 3) * vect_to_project(0, coord + 3);
		//	}
		//	prod/= sqrt(norm);
		//	std::cout << prod << " " << prod1 << "\n";
		//}
	}
}

}