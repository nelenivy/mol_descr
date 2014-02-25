#pragma once
#include <utility>
#include <vector>
#include <fstream>
#include <array>
#include "boost/archive/polymorphic_text_oarchive.hpp"

#include "opencv2/core/core.hpp"
#include "shark/Algorithms/Trainers/CSvmTrainer.h"
#include "shark/Models/Kernels/KernelExpansion.h"
#include "shark/ObjectiveFunctions/Loss/ZeroOneLoss.h"
#include <shark/Data/CVDatasetTools.h>
#include <shark/Algorithms/DirectSearch/GridSearch.h>
#include <shark/ObjectiveFunctions/CrossValidationError.h>
#include "shark/Data/Csv.h"

#include "../Common/singular_point.h"
#include "Kernels/pharm_kernel_seq_calc.h"
#include "element_with_index_and_id.h"
#include "modified_grid_search.h"
#include "prepare_data_for_svm.h"
#include "zero_one_loss_fuzzy.h"

namespace boost
{
	namespace serialization
	{
		template <class Archive/*, typename T*/>
		void serialize(Archive& ar, /*const*/ cv::Point3d/*_<T>*/& point, const unsigned int version)
		{
			ar & point.x;
			ar & point.y;
			ar & point.z;
		}

		template <class Archive, typename ElemType/*, typename T*/>
		void serialize(Archive& ar, /*const*/ /*_<T>*/molecule_descriptor::ElemWithIndexAndID<ElemType>& point, const unsigned int version)
		{
			ar & point.Elem();
		}

		template <class Archive, typename ElemType/*, typename T*/>
		void serialize(Archive& ar, /*const*/ /*_<T>*/molecule_descriptor::SingularPoint<ElemType>& point, const unsigned int version)
		{
			ar & point.GetAsPair();
		}

		template <class Archive/*, typename T*/>
		void serialize(Archive& ar, /*const*/ /*_<T>*/molecule_descriptor::PropertiesSet& point, const unsigned int version)
		{
			ar & point.SurfaceType();
			ar & point.Charge();
			ar & point.ElectricPotential();
		}

		template <class Archive, typename T, size_t kSize>
		void serialize(Archive& ar, /*const*/ /*_<T>*/std::array<T, kSize>& point, const unsigned int version)
		{
			for (size_t ind = 0; ind < point.size(); ++ind)
			{
				ar & point[ind];
			}
		}
	}
}

namespace molecule_descriptor
{

using std::vector;
using std::string;
using  shark::RealVector;
//SingularPoint<> is a pair of 3d point and its PropType
//PropType is used for kernel PropKernel
template <typename PropType, class DistKernel, class PropKernel>
class KernelSVMTrainerManager
{
public:
	typedef SingularPoint<PropType> singular_point;
	typedef std::vector<singular_point> sing_pts_seq;
	typedef ElemWithIndexAndID<sing_pts_seq> sing_pts_seq_wth_ind;

	KernelSVMTrainerManager() 
		: m_kernel_expansion(true), 
		m_curr_dataset_id(0)
	{	}

	void SetData(const vector<vector<singular_point>>& data,
		const vector<unsigned int>& labels);//Convert input to shark format
	void SetKernels(DistKernel dist_kernel, PropKernel prop_kernel);
	void Train(const std::string& file_name);
	void Write(const std::string& file_name);
private:
	void EvaluateTrained();

	int m_curr_dataset_id;
	shark::LabeledData<sing_pts_seq_wth_ind, unsigned int> m_labeled_data;
	PharmSequenceKernelPairs<PropType, DistKernel, PropKernel> m_pairs_kernel;	
	shark::KernelExpansion<sing_pts_seq_wth_ind> m_kernel_expansion;

	shark::Data<shark::blas::vector<double>> m_eval_result;
	double m_training_error;
};

template <typename PropType, class DistKernel, class PropKernel>
void KernelSVMTrainerManager<PropType, DistKernel, PropKernel>::SetData(const vector<vector<singular_point>>& data, const vector<unsigned int>& labels)
{
	//convert data to shark format
	std::vector<sing_pts_seq_wth_ind> data_wth_ind;
	CreateSeqWithIndexAndIdFromSeq(data.begin(), data.end(), m_curr_dataset_id, data_wth_ind);
	++m_curr_dataset_id;
	m_labeled_data = PrepareDataForSVM(data_wth_ind.begin(), data_wth_ind.end(), labels.begin(), labels.end());
	auto s = m_labeled_data.inputs().elements();
	/*for (auto it = s.begin(); it != s.end(); ++it)
	{
		auto d = it->ElemConst();
		for (auto it1 = d.begin(); it1 != d.end(); ++it1)
		{
			std::cout << it1->Property().SurfaceType();
		}
		std::cout << "\n";


	}*/
}

template <typename PropType, class DistKernel, class PropKernel>
void KernelSVMTrainerManager<PropType, DistKernel, PropKernel>::SetKernels(DistKernel dist_kernel, PropKernel prop_kernel)
{
	PharmSequenceKernelPairs<PropType, DistKernel, PropKernel> pairs_kernel(dist_kernel, prop_kernel);
	m_pairs_kernel = pairs_kernel;
}

template <typename PropType, class DistKernel, class PropKernel>
void KernelSVMTrainerManager<PropType, DistKernel, PropKernel>::Train(const std::string& file_name)
{
	//trainer
	const double C = 0.001;
	shark::CSvmTrainer<sing_pts_seq_wth_ind> svm_trainer(&m_pairs_kernel, C, true, false);
	std::cout << svm_trainer.parameterVector();
	std::cout << m_pairs_kernel.parameterVector();
	// cross-validation error
	const unsigned int N= 5;//m_labeled_data.elements().size();  // number of folds
	shark::ZeroOneLoss/*Fuzzy*/<unsigned int, shark::RealVector> loss;//(-1.0, 1.0);
	shark::CVFolds<shark::LabeledData<sing_pts_seq_wth_ind, unsigned int>> folds = shark::createCVSameSizeBalanced(m_labeled_data, N);
	shark::CrossValidationError<shark::KernelExpansion<sing_pts_seq_wth_ind>, unsigned int> cv_error(
		folds, &svm_trainer, &m_kernel_expansion, &svm_trainer, &loss);
	// find best parameters
	std::vector<double> min(2);
	std::vector<double> max(2);
	std::vector<size_t> sections(2);
	min[1] = -8; max[1] = -1; sections[1] = 8;  // regularization parameter C
	min[0] = 0.00001; max[0] = 0.00010; sections[0] = 5;   // kernel parameter gamma
	GridSearchFromEnd grid;
	grid.configure(min, max, sections);
	std::ofstream out(file_name + "_res_001_01_all_triangle_triangle.txt");
	grid.SetOutput(out);
	grid.step(cv_error);
	// train model on the full dataset
	svm_trainer.setParameterVector(grid.solution().point);
	std::cout << grid.solution().value << "\n";
	std::cout << grid.solution().point << "\n";
	svm_trainer.train(m_kernel_expansion, m_labeled_data);

	EvaluateTrained();
}

template <typename PropType, class DistKernel, class PropKernel>
void KernelSVMTrainerManager<PropType, DistKernel, PropKernel>::EvaluateTrained()
{
	m_eval_result = m_kernel_expansion(m_labeled_data.inputs());
	shark::ZeroOneLoss<unsigned int, shark::RealVector> loss;
	m_training_error = loss.eval(m_labeled_data.labels(), m_eval_result);
}

template <typename PropType, class DistKernel, class PropKernel>
void KernelSVMTrainerManager<PropType, DistKernel, PropKernel>::Write(const std::string& file_name)
{
	std::ofstream out(file_name + "_res_001_01_triangle_triangle.txt");
	boost::archive::polymorphic_text_oarchive oa(out);
	m_eval_result.write(oa);
	out.close();

	out.open(file_name + "_res_normal_001_01_triangle_triangle.txt");
	for (auto iter= m_eval_result.elements().begin(); iter != m_eval_result.elements().end() ; ++iter)
	{
		out << *iter << " ";
	}
	out.close();

	out.open(file_name + "_labels_in_001_01_triangle_triangle.txt");
	boost::archive::polymorphic_text_oarchive oa3(out);
	m_labeled_data.labels().write(oa);
	out.close();

	out.open(file_name + "_svm_001_01_triangle_triangle.txt");
	boost::archive::polymorphic_text_oarchive oa1(out);
	m_kernel_expansion.write(oa1);
	out.close();
	out.open(file_name + "_res_set_001_01_triangle_triangle.txt");
	out << m_training_error;
}

}