#pragma once
#include <utility>
#include <vector>

#include "opencv2/core/core.hpp"
#include "shark/Models/Kernels/KernelExpansion.h"
//namespace boost
//{
//	namespace serialization
//	{
//		template <class Archive/*, typename T*/>
//		void serialize(Archive& ar, /*const*/ cv::Point3d/*_<T>*/& point, const unsigned int version)
//		{
//			ar & point.x;
//			ar & point.y;
//			ar & point.z;
//		}
//
//		template <class Archive, typename ElemType/*, typename T*/>
//		void serialize(Archive& ar, /*const*/ /*_<T>*/molecule_descriptor::ElemWithIndexAndID<ElemType>& point, const unsigned int version)
//		{
//			ar & point.Elem();
//		}
//
//		template <class Archive, typename ElemType/*, typename T*/>
//		void serialize(Archive& ar, /*const*/ /*_<T>*/molecule_descriptor::SingularPoint<ElemType>& point, const unsigned int version)
//		{
//			ar & point.GetAsPair();
//		}
//
//		template <class Archive/*, typename T*/>
//		void serialize(Archive& ar, /*const*/ /*_<T>*/molecule_descriptor::PropertiesSet& point, const unsigned int version)
//		{
//			ar & point.SurfaceType();
//			ar & point.Charge();
//			ar & point.ElectricPotential();
//		}
//	}
//}

namespace molecule_descriptor
{
	using std::vector;
	using std::string;
	using shark::RealVector;

	class DescriptorSVMTrainerManager
	{
	public:
		typedef size_t singular_point;
		typedef /*std::vector<singular_point>*/RealVector sing_pts_seq;
		typedef sing_pts_seq/*ElemWithIndexAndID<sing_pts_seq>*/ sing_pts_seq_wth_ind;

		DescriptorSVMTrainerManager() 
			: m_kernel_expansion(true), 
			m_curr_dataset_id(0)
		{	}

		void SetData(const cv::Mat_<size_t>& data,
			const vector<unsigned int>& labels);//Convert input to shark format
		void Train(const std::string& file_name);
		void Write(const std::string& file_name);
	private:
		void EvaluateTrained();

		int m_curr_dataset_id;
		shark::LabeledData<sing_pts_seq_wth_ind, unsigned int> m_labeled_data;
		shark::KernelExpansion<sing_pts_seq_wth_ind> m_kernel_expansion;
		shark::Data<shark::blas::vector<double>> m_eval_result;
		double m_training_error;
	};
}