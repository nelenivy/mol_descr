#include "set_manager.h"

#include <sstream>
#include <fstream>
#include <iostream>

#include <utility>
#include <functional>
#include <limits>

#include "opencv2/core/core.hpp"
#include "shark\LinAlg\BLAS\matrix.hpp"
#include "shark/Models/Kernels/LinearKernel.h"

#include "extensions.h"
#include "InputOutput/params_reader.h"
#include "InputOutput/surface_reader.h"
#include "../Common/interval_read_write.h"
#include "CommonUtilities/classic_k_means.h"
#include "types_calculation.h"
#include "SingularPoints/i_singular_points_finder.h"

#include "../../SVM/kernel_svm_trainer_manager.h"
#include "../../SVM/descriptor_svm_trainer.h"
#include "../../SVM/cross_validation_svm_trainer.h"
#include "../../SVM/pharm_kernel_naive_configure.h"

#include "../../SVM/Kernels/basic_kernels.h"
#include "../../SVM/Kernels/kernels_for_points.h"
namespace molecule_descriptor
{

const bool kCalculate = true;
const bool kNotCalculate = false;

void SetManager::Init(int argc, char** argv)
{
	for (int ind = 1; ind < argc; ++ind)
	{
		m_command_line += argv[ind];
	}
	ReadParamsFromCommandLine(argc, argv);
	FindOutMoleculesNum();
	m_molecule_manager.SetLevelsNum(m_mesh_levels_num);
	m_molecule_manager.SetPairsLevelsOverlap(m_pairs_levels_overlap);

	if (m_use_levels)
	{
		m_molecule_manager.SetSingPtsAlgorithm(SCALE_SPACE);
	}
	else
	{
		m_molecule_manager.SetSingPtsAlgorithm(SEGMENTATION);
	}

	m_molecule_manager.Init(argc, argv);
}

template <typename Iter>
double CalcMeanAndDev(const Iter it_beg, const Iter it_end, double& mean, double& dev)
{
	mean = 0;
	int length = 0;

	for (Iter iter = it_beg; iter != it_end; ++iter)
	{
		mean += *iter;
		++length;
	}

	if (length == 0)
	{
		return 0;
	}
	mean /= static_cast<double>(length);

	dev = 0;
	double mean_diff = 0;
	for (Iter iter = it_beg; iter != it_end; ++iter)
	{
		dev += (*iter - mean) * (*iter - mean);
		mean_diff += std::abs(*iter - mean);
	}
	dev /= static_cast<double>(length);
	mean_diff /= static_cast<double>(length);
	dev = sqrt(dev);
	const double span = *(it_end - 1) - *it_beg;
	return (6.0 * dev - span) / span;
}

void SetManager::ProcessSingularPoints(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	if (calculate)
	{
		std::cout << "Collecting properties\n";
		std::vector<std::vector<double>> props;

		for (int ind = 0; ind < m_molecules_num; ++ind)
		{
			std::cout << ind << "\n";
			m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
			m_molecule_manager.CollectProperties(props);			
		}
		std::vector<std::vector<double>> prop_params(props.size());
		for (size_t ind = 0; ind < props.size(); ++ind)
		{			
			prop_params[ind].resize(2);
			//cut too high props
			int nearest_to_zer_lev = 0;
			double least_to_zero = std::numeric_limits<double>::max();
			for (int lev = 0; lev < 50; ++ lev)
			{
				const int low_thresh_ind = (0.01 * lev) * props[ind].size();
				const int high_thresh_ind = (1.0 -0.01 * lev)  * props[ind].size();
				std::nth_element(props[ind].begin(), props[ind].begin() + low_thresh_ind, props[ind].end());
				std::nth_element(props[ind].begin(), props[ind].begin() + high_thresh_ind, props[ind].end());

				const double curr_ratio_to_span = CalcMeanAndDev(props[ind].begin() + low_thresh_ind, 
					props[ind].begin() + high_thresh_ind, prop_params[ind][0], prop_params[ind][1]);

				if (abs(curr_ratio_to_span) < abs(least_to_zero))
				{
					nearest_to_zer_lev = lev;
					least_to_zero = curr_ratio_to_span;
				}
			}
			const int low_thresh_ind = (0.01 * nearest_to_zer_lev) * props[ind].size();
			const int high_thresh_ind = (1.0 -0.01 * nearest_to_zer_lev)  * props[ind].size();
			std::nth_element(props[ind].begin(), props[ind].begin() + low_thresh_ind, props[ind].end());
			std::nth_element(props[ind].begin(), props[ind].begin() + high_thresh_ind, props[ind].end());

			const double curr_ratio_to_span = CalcMeanAndDev(props[ind].begin() + low_thresh_ind, 
				props[ind].begin() + high_thresh_ind, prop_params[ind][0], prop_params[ind][1]);

			std::cout << "mean = " << prop_params[ind][0] << "dev = " << prop_params[ind][1] 
			<< " 3 sigma " << curr_ratio_to_span
				<< " lev " << nearest_to_zer_lev << std::endl;
		}
		m_molecule_manager.SetMeansAndSigma(prop_params);
	}

	std::cout << "Singular Points Processing\n";
	std::vector<double> charges, charges_thresh;//for charges classification
	std::vector<double> lennard_jones, lennard_jones_thresh;//for lennard-jones classification
	std::vector<double> areas, areas_thresh;//for areas classification

	const bool calc_as_average = false;
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.FindSingularPoints(calculate, calc_as_average);
		m_molecule_manager.AppendCharges(charges, false/*true*/);
		m_molecule_manager.AppendLennardJones(lennard_jones, false/*true*/);
		m_molecule_manager.AppendLennardJones(areas, false/*true*/);
	}

	CalculateThresholdsQuantiles(m_charges_threshes_num, charges, charges_thresh);
	CalculateThresholdsQuantiles(m_lennard_jones_threshes_num, lennard_jones, lennard_jones_thresh);
	CalculateThresholdsQuantiles(m_areas_threshes_num, areas, areas_thresh);
	m_molecule_manager.SetChargesThresholds(charges_thresh);
	m_molecule_manager.SetLennardJonesThresholds(lennard_jones_thresh);
	m_molecule_manager.SetAreaThresholds(areas_thresh);

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ClassifySingularPoints();
	}
}

void SetManager::ProcessPairsLevels(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	m_distances_levels.clear();
	m_dist_thresholds.resize(m_distances_intervals);
	std::cout << "Pairs Processing\n";

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairsLevels(calculate, m_write_pairs);
		m_molecule_manager.AppendDistancesLevels(m_distances_levels);
	}

	//CalculateDistThresholdsKMeans();
	CalculateThresholdLevels();
	m_molecule_manager.SetDistLevelsThresholds(m_dist_thresholds_levels, m_high_thresholds);
	std::vector<double> levels_scales(m_dist_thresholds_levels.size());
	levels_scales[0] = 1.0;
	for (int ind = 0; ind < levels_scales.size(); ++ind)
	{
		levels_scales[ind] = 1.0;//static_cast<double>(m_distances_levels[0].size()) / m_distances_levels[ind].size();
	}
	m_molecule_manager.SetLevelsScales(levels_scales);
	const int pairs_types_num = static_cast<int>(m_molecule_manager.GetPairsTypeNumLevels());
	m_md_matrix_double.create(m_molecules_num, pairs_types_num);
	//calculate distance types and write them to the matrix
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairsLevels(!m_write_pairs, false);
		m_molecule_manager.CalculatePairsWithTypesLevels();
		//m_molecule_manager.CalculatePairsWithTypesLevelsAllMeshSmoothedCurv();
		m_molecule_manager.GetPairsHistogrammLevels(m_md_matrix_double.row(ind));
	}

	m_processed = true;
	WriteMatrix(m_md_matrix_double, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrixLevels());
}

void SetManager::ProcessTriplesLevels(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	m_distances_levels.clear();
	m_dist_thresholds.resize(m_distances_intervals);
	std::cout << "Pairs Processing\n";

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairsLevels(calculate, m_write_pairs);
		m_molecule_manager.AppendDistancesLevels(m_distances_levels);
	}

	//CalculateDistThresholdsKMeans();
	CalculateThresholdLevels();
	m_molecule_manager.SetDistLevelsThresholds(m_dist_thresholds_levels, m_high_thresholds);
	std::vector<double> levels_scales(m_dist_thresholds_levels.size());
	levels_scales[0] = 1.0;
	for (int ind = 0; ind < levels_scales.size(); ++ind)
	{
		levels_scales[ind] = 1.0;//static_cast<double>(m_distances_levels[0].size()) / m_distances_levels[ind].size();
	}
	m_molecule_manager.SetLevelsScales(levels_scales);
	const int triples_types_num = static_cast<int>(m_molecule_manager.GetTriplesTypeNumLevels());
	std::cout << triples_types_num << "\n";
	m_for_large_matrix.resize(m_molecules_num * triples_types_num);
	m_md_matrix_float = cv::Mat_<float>(m_molecules_num, triples_types_num, m_for_large_matrix.data());
	//calculate distance types and write them to the matrix
	std::cout << "Pairs Processing\n";

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";

		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindTriplesLevels(!m_write_pairs, false);
		m_molecule_manager.CalculateTriplesWithTypesLevels();
		//m_molecule_manager.CalculatePairsWithTypesLevelsAllMeshSmoothedCurv();
		m_molecule_manager.GetTriplesHistogrammLevels(m_md_matrix_float.row(ind));
	}

	m_processed = true;
	//WriteMatrix(m_md_matrix_float, m_mol_folder + m_mol_prefix + Extensions::TriplesMDMatrixLevels());
}
void SetManager::ProcessPairs(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	m_distances.clear();
	m_dist_thresholds.resize(m_distances_intervals);
	std::cout << "Pairs Processing\n";

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairs(calculate);
		m_molecule_manager.AppendDistances(m_distances);
	}

	//CalculateDistThresholdsKMeans();
	CalculateThresholdsLevels(m_distances_intervals, 0.05, m_distances, m_dist_thresholds);
	//CalculateThresholdsQuantiles(m_distances_intervals, m_distances, m_dist_thresholds);
	m_molecule_manager.SetDistThresholds(m_dist_thresholds);
	const int pairs_types_num = static_cast<int>(m_molecule_manager.GetPairsTypeNum());
	m_md_matrix.create(m_molecules_num, pairs_types_num);
	//calculate distance types and write them to the matrix
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindPairs(kNotCalculate);
		m_molecule_manager.CalculatePairsWithTypes();
		m_molecule_manager.GetPairsHistogramm(m_md_matrix.row(ind));
	}

	m_processed = true;
	WriteMatrix(m_md_matrix, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrix());
}

void SetManager::CalculateThresholdLevels()
{
	const int levels_num = m_distances_levels.size();
	CV_Assert(levels_num > 1);
	m_high_thresholds.resize(levels_num, std::numeric_limits<double>::max());
	m_dist_thresholds_levels.resize(levels_num);
	m_high_thresholds[levels_num - 1] = std::numeric_limits<double>::max();
	CV_Assert(m_distance_quantile_for_thresh <= m_distance_interval_levels);

	for (int level = levels_num - 1; level >= 0; --level)
	{
		//delete too large distances
		const auto first_elem_higher_than_thresh = std::partition(m_distances_levels[level].begin(), m_distances_levels[level].end(),
			[&](const double dist)->bool
		{
			return dist < m_high_thresholds[level];
		});
		m_distances_levels[level].erase(first_elem_higher_than_thresh, m_distances_levels[level].end());
		std::vector<double> thresh_levels, thresh_quantiles;
		CalculateThresholdsLevels(m_distance_interval_levels, 0.02, m_distances_levels[level], thresh_levels);
		CalculateThresholdsQuantiles(m_distance_interval_levels, m_distances_levels[level], thresh_quantiles);
		const double alpha = 1.0;
		m_dist_thresholds_levels[level].resize(m_distance_interval_levels);
		for (size_t ind = 0; ind < m_distance_interval_levels; ++ind)
		{
			m_dist_thresholds_levels[level][ind] = alpha * thresh_levels[ind] + (1.0  - alpha) * thresh_quantiles[ind];
		}
		if (level > 0 /*&& level < 8*/)
		{
			m_high_thresholds[level - 1] = m_dist_thresholds_levels[level][m_distance_quantile_for_thresh];
		}
	}
}

void SetManager::ProcessTriples(const bool calculate)
{
	CV_Assert(m_inited);
	//calculate pairs and write distances
	std::cout << "Triples Processing\n";
	const int triples_types_num = static_cast<int>(m_molecule_manager.GetTriplesTypeNum());
	m_md_matrix.create(m_molecules_num, triples_types_num);

	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.FindTriples(calculate);
		m_molecule_manager.CalculateTriplesWithTypes();
		m_molecule_manager.GetTriplesHistogramm(m_md_matrix.row(ind));
	}

	m_processed = true;
	WriteMatrix(m_md_matrix, m_mol_folder + m_mol_prefix + Extensions::TriplesMDMatrix());
}

void SetManager::ProcessKernelSVMPoints()
{
	std::vector<std::vector<NonMarkedSingularPoint>> pts(m_molecules_num);

	for (int ind = 0; ind < m_molecules_num; ind++)
	{
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.GetNonMarkedSingularPoints(pts[ind]);
	}

	std::vector<int> labels;
	ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}

	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	typedef TriangleKernel/*AlwaysOneKernel<float>*/ PotentialKernel;
	typedef KernelForParametersSet<DiracKernel<size_t>, /*AlwaysOneKernel*/DiracKernel<int>, PotentialKernel> kernel_type;	
	PotentialKernel potential_kernel(0.2);//(0.020, -1.0, 1.0);
	kernel_type tuple_kernel(
		DiracKernel<size_t>(),/*AlwaysOneKernel*/DiracKernel<int>(), potential_kernel);
	KernelSVMTrainerManager<PropertiesSet, GaussianKernelOneDim, kernel_type> svm_trainer;
	svm_trainer.SetData(pts, unsigned_labels);
	std::cout << tuple_kernel.parameterVector();
	svm_trainer.SetKernels(GaussianKernelOneDim(m_dist_thresholds), tuple_kernel);
	svm_trainer.Train(m_mol_folder + m_mol_prefix);
	svm_trainer.Write(m_mol_folder + m_mol_prefix);
}

void SetManager::ProcessKernelSVMPointsWithFiltering()
{	//Calculate good descriptors
	ProcessMGUASVMClassification(1, 1);
	const size_t descriptors_arr[] = {1369, 1542, 1381};//,
	/*1295, 157, 1009, 1048, 118, 38,
	1038, 14, 128, 176, 872, 135,
	901, 37, 271, 1048, 118, 38, 
	1157, 180, 1009, 1048, 118, 38,
	185, 6, 116, 1048, 118, 38,
	180, 157, 1009, 1048, 118, 38 ,
	415, 157, 1009, 1048, 118, 38,
	1053, 14, 128, 176, 872, 135,
	1295, 258, 128, 176, 872, 135};*/
	
	const size_t kGoodSize = sizeof(descriptors_arr) / sizeof(descriptors_arr[0]);
	std::vector<bool> good_descriptors = m_mgua_trainer.GetGoodDescriptors();	
	//ProcessPairs must be called
	//read pairs
	CV_Assert(m_inited);
	//calculate all pairs
	//////////////////////////////////////////////////////////////////////////
	std::cout << "Pairs Processing\n";
	std::vector<NonMarkedSingularPoint> sing_pts;
	std::vector<std::vector<SingularPointsPair<PropertiesSet>>> sing_pts_pairs(m_molecules_num);
	std::vector<std::vector<size_t>> sing_pts_pairs_types(m_molecules_num);
	std::vector<size_t> curr_pairs_types;
	std::vector<SingularPointsPair<PropertiesSet>> curr_non_marked_pairs;
	std::vector<std::vector<size_t>> new_descr_mat_sparsed(m_molecules_num);
	std::vector<std::vector<size_t>> new_descr_mat(m_molecules_num);
	for (int ind = 0; ind < m_molecules_num; ++ind)
	{
		new_descr_mat_sparsed[ind].resize(kGoodSize, 0);
		new_descr_mat[ind].resize(m_molecule_manager.GetPairsTypeNum(), 0);
		size_t good_pairs = 0;
		size_t bad_pairs = 0;
		std::vector<SingularPointsPair<PropertiesSet>>& curr_pairs = sing_pts_pairs[ind];
		std::cout << ind << "\n";
		m_molecule_manager.SetCurrFilePrefix(MakeMoleculePrefix(ind));
		m_molecule_manager.ReadAllSingularPoints(m_mesh_levels_num);
		m_molecule_manager.GetNonMarkedSingularPoints(sing_pts);
		m_molecule_manager.FindPairs(kNotCalculate);
		m_molecule_manager.CalculatePairsWithTypes();
		m_molecule_manager.GetPairsTypes(curr_pairs_types);
		m_molecule_manager.GetNonMarkedPairs(curr_non_marked_pairs);

		for (size_t pair_num = 0; pair_num < curr_pairs_types.size(); ++pair_num)
		{
			const size_t curr_descr = curr_pairs_types[pair_num];
			const size_t* pos = std::find(descriptors_arr, descriptors_arr + kGoodSize, curr_descr);
			new_descr_mat[ind][curr_descr]++;
			if (descriptors_arr + kGoodSize == pos)
			//if (!good_descriptors[curr_descr])
			{
				bad_pairs++;
				continue;
			}
			good_pairs++;
			sing_pts_pairs[ind].push_back(curr_non_marked_pairs[pair_num]);
			new_descr_mat_sparsed[ind][pos - descriptors_arr]++;
		}

		std::cout << "molecule num " << ind << " bad_pairs " << bad_pairs << " good_pairs " << good_pairs << "\n";
	}
	//////////////////////////////////////////////////////////////////////////
	//get trainer
	typedef /*AlwaysOneKernel<float>*/TriangleKernel PotentialKernel;
	typedef KernelForParametersSet<DiracKernel<size_t>, /*AlwaysOneKernel*/DiracKernel<int>, PotentialKernel> kernel_type;	
	PotentialKernel potential_kernel(0.03);//, -1.0, 1.0);
	kernel_type tuple_kernel(
		DiracKernel<size_t>(),/*AlwaysOneKernel*/DiracKernel<int>(), potential_kernel);
	auto trainer = GetNaiveKernelSVMTrainerForClassification<PropertiesSet>(tuple_kernel, GaussianKernelOneDim(0.3/*m_dist_thresholds*/));
	std::vector<double> min(2);
	std::vector<double> max(2);
	std::vector<size_t> sections(2);
	min[1] = -7; max[1] = 7; sections[1] = 15;  // regularization parameter C
	min[0] = 0.10; max[0] = 1.0; sections[0] = 10;   // kernel parameter gamma
	trainer.SetGrid(min, max, sections);
	//read labels
	std::vector<int> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}
	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	//convert data to shark format
	shark::Data<std::vector<SingularPointsPair<PropertiesSet>>> input_data = shark::createDataFromRange(sing_pts_pairs);
	shark::Data<unsigned int> input_labels = shark::createDataFromRange(unsigned_labels);
	auto labeled_data = shark::LabeledData<std::vector<SingularPointsPair<PropertiesSet>>, unsigned int>(input_data, input_labels);	
	std::vector<ElemWithIndexAndID<std::vector<SingularPointsPair<PropertiesSet>>>> index_elems;
	CreateSeqWithIndexAndIdFromSeq(sing_pts_pairs.begin(), sing_pts_pairs.end(), 0, index_elems);
	{
		std::ofstream t("m5.txt");
		std::ofstream t1("md_new.txt");
		for (size_t y = 0; y < labeled_data.numberOfElements(); ++y)
		{
			for (int ind = 0; ind < new_descr_mat_sparsed[y].size(); ++ind)
			{
				t << new_descr_mat_sparsed[y][ind] << " ";
			}
			t << "\n";
			for (int ind = 0; ind < new_descr_mat[y].size(); ++ind)
			{
				t1 << new_descr_mat[y][ind] << " ";
			}
			t1 << "\n";
			for (size_t x = 0; x < labeled_data.numberOfElements(); ++x)
			{
				t << y << " " << x << " " << (*trainer.m_kernel)(index_elems[y], index_elems[x]) << " ";
				t << "\n";
			}
		}
	}

	trainer.train(labeled_data, labeled_data.numberOfElements() / 2);
}

void SetManager::ProcessDescriptorsSVM()
{
	//read MD matrix
	molecule_descriptor::ReadMatrix(m_md_matrix, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrix());
	//read labels
	std::vector<int> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}


	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	DescriptorSVMTrainerManager svm_trainer;
	svm_trainer.SetData(m_md_matrix, unsigned_labels);
	svm_trainer.Train(m_mol_folder + m_mol_prefix);
	svm_trainer.Write(m_mol_folder + m_mol_prefix);
}



void SetManager::ProcessMGUASVMClassification(const int iter_num, const int descr_num)
{
	//read MD matrix
	const std::string md_mat_name = m_mol_folder + m_mol_prefix + 
		(m_use_triples ? 
		(m_use_levels ? Extensions::PairsMDMatrixLevels() : Extensions::PairsMDMatrix())
		: (m_use_levels ? Extensions::TriplesMDMatrixLevels() : Extensions::TriplesMDMatrix()));
	molecule_descriptor::ReadMatrix(m_md_matrix_double, md_mat_name);
	//read labels
	std::vector<int> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

	for (auto it = labels.begin(); it != labels.end(); ++it)
	{
		*it = *it == -1 ? 0 : *it;
	}

	std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
	auto svm_trainer = GetSVMTrainerForLinearClassification<shark::RealVector>();
	typedef decltype(svm_trainer) SVMType;
	const int kThreadsNum = 7;
	std::vector<SVMType> svm_vect(kThreadsNum);
	std::vector<double> param_min(1, -2.0), param_max(1, 1.0);
	std::vector<size_t> param_num(1, 4);

	for (size_t ind = 0; ind < svm_vect.size(); ++ind)
	{
		svm_vect[ind] = GetSVMTrainerForLinearClassification<shark::RealVector>();
		svm_vect[ind].SetGrid(param_min, param_max, param_num);
		svm_vect[ind].SetOutput(m_mol_folder + m_mol_prefix + "trainer_output1.txt");
	}

	m_mgua_trainer = MGUATrainer<unsigned int>();
	m_mgua_trainer.SetData(m_md_matrix_double, unsigned_labels, true);
	m_mgua_trainer.SetParameters(iter_num, descr_num);
	//m_mgua_trainer.TrainOneSet(svm_vect[0], m_mol_folder + m_mol_prefix);
	m_mgua_trainer.Train(svm_vect, m_mol_folder + m_mol_prefix,true);
}

void SetManager::ProcessSVMClassificationL0()
{
	
	const std::string md_mat_name = m_mol_folder + m_mol_prefix + 
		(! m_use_triples ? 
		(m_use_levels ? Extensions::PairsMDMatrixLevels() : Extensions::PairsMDMatrix())
		: (m_use_levels ? Extensions::TriplesMDMatrixLevels() : Extensions::TriplesMDMatrix()));

	if (!m_use_triples)			
	{
		molecule_descriptor::ReadMatrix(m_md_matrix_double, md_mat_name);
	}

	if (1)
	{
		for (int x = 0; x < m_md_matrix_double.cols; ++x)
		{
			double mean = 0;
			for (int y = 0; y < m_md_matrix_double.rows; ++y)
			{
				mean += m_md_matrix_double(y, x);
			}
			mean /= m_md_matrix_double.rows;
			double dev = 0;
			for (int y = 0; y < m_md_matrix_double.rows; ++y)
			{
				dev += (mean - m_md_matrix_double(y, x)) * (mean - m_md_matrix_double(y, x));
			}
			dev /= m_md_matrix_double.rows;
			for (int y = 0; y < m_md_matrix_double.rows; ++y)
			{
				m_md_matrix_double(y, x) = (m_md_matrix_double(y, x) - mean);
				if (dev > 0.0)
				{
					m_md_matrix_double(y, x) /= sqrt(dev);
				}
			}
		}
	}

	cv::Mat_<double> md_extended;
	std::vector<double> for_large_mat_extended;
	if (1)//!m_use_triples)
	{
		md_extended.create(m_md_matrix_double.rows, m_md_matrix_double.rows + m_md_matrix_double.cols);
	}
	else
	{
		for_large_mat_extended.resize(m_md_matrix_double.rows * (m_md_matrix_double.rows + m_md_matrix_double.cols));
		md_extended = cv::Mat_<double>(m_md_matrix_double.rows, m_md_matrix_double.rows + m_md_matrix_double.cols, for_large_mat_extended.data());
	}

	for (int lambda = 0; lambda <= 9; lambda += 3)//read MD matrix
	{
		for (int y = 0; y < m_md_matrix_double.rows; ++y)
		{
			for (int x = 0; x < m_md_matrix_double.cols; ++x)
			{
				md_extended(y, x) = m_md_matrix_double(y, x);
			}

			for (int x = 0; x < m_md_matrix_double.rows; ++x)
			{
				md_extended(y, x + m_md_matrix_double.cols) =  (x == y ? 1.0 / exp(lambda) : 0);
			}
		}
		cv::Mat_<double>& mat_to_use = md_extended;
				
		//convert data to shark format
		//read labels
		std::vector<int> labels;
		molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Labels(), labels);

		for (auto it = labels.begin(); it != labels.end(); ++it)
		{
			*it = *it == -1 ? 0 : *it;
		}

		std::vector<unsigned int> unsigned_labels(labels.begin(), labels.end());
		shark::LinearCSvmTrainer<shark::RealVector> svm_trainer(10000.0, true);
		std::vector<int> true_descriptor_nums(mat_to_use.cols);
		int non_zero_descr = 0;
		for (int x = 0; x < mat_to_use.cols; ++x)
		{
			int non_zero_elems = 0;
			for (int y = 0; y < mat_to_use.rows; ++y)
			{
				non_zero_elems += (mat_to_use(y, x) > 0);
			}
			if (non_zero_elems > 0)
			{
				true_descriptor_nums[non_zero_descr] = x;
				++non_zero_descr;
			}
		}
		true_descriptor_nums.erase(true_descriptor_nums.begin() + non_zero_descr, true_descriptor_nums.end());
		std::vector<double> weights(true_descriptor_nums.size(), 1.0), old_weights(weights);

		double diff = std::numeric_limits<double>::max();
	    int non_changed_iterations = 0;

		while (diff > 0.000005 && non_changed_iterations < 5 )
		{
			typedef shark::RealVector data_type;
			std::vector<data_type> temp(mat_to_use.rows, data_type(weights.size()));
			//rescale data with new weights
			for (int y = 0; y < mat_to_use.rows; ++y)
			{
				for (int x = 0; x < weights.size(); ++x)
				{
					temp[y][x] = mat_to_use(y, true_descriptor_nums[x]) * weights[x];
				}		
			}
			//train new weights
			auto balanced_data = PrepareDataForSVM(temp.begin(), temp.end(), unsigned_labels.begin(), unsigned_labels.end());
			shark::LinearClassifier<shark::RealVector> model;
			/*const double err = */svm_trainer.train(model, balanced_data);
			std::vector<double> new_weights(model.parameterVector().begin(), model.parameterVector().end());
			//calculate new weights characteristics
			CV_Assert(new_weights.size() == weights.size());
			diff = 0;
			int non_zero_weight = 0;
			int norm = 0;
			int changes = 0;
			for (size_t ind = 0; ind < weights.size(); ++ind)
			{
				diff += std::abs(old_weights[ind] - new_weights[ind]);
				norm += std::abs(old_weights[ind]);
				non_zero_weight += abs(new_weights[ind]) > 0.0;
				changes += (abs(new_weights[ind]) <= 0.0 && abs(old_weights[ind]) > 0.0);
			}
			non_changed_iterations = changes > 0 ? 0 : non_changed_iterations + 1;
			std::cout << changes << " " << non_zero_weight << " " << diff << " " << diff /norm /*<< err */<< "\n";
			//write new weights, delete zero weights
			non_zero_weight = 0;
			for (size_t ind = 0; ind < weights.size(); ++ind)
			{
				if (std::abs(new_weights[ind]) > 0.0)
				{
					weights[non_zero_weight] = weights[ind] * new_weights[ind];
					true_descriptor_nums[non_zero_weight] = true_descriptor_nums[ind];
					++non_zero_weight;
				}
			}
			weights.erase(weights.begin() + non_zero_weight, weights.end());
			true_descriptor_nums.erase(true_descriptor_nums.begin() + non_zero_weight, true_descriptor_nums.end());
			old_weights = new_weights;
		}
		//evaluate obtained decision
		int true_descriptors = 0;//non-zero weights which belong to true descriptors, not additional
		for (int ind = 0; ind < true_descriptor_nums.size(); ++ind)
		{
			true_descriptors += (true_descriptor_nums[ind] < m_md_matrix_double.cols);
		}
		std::cout << true_descriptors << "\n";
		//write non-zero data
		std::vector<shark::RealVector> temp(m_md_matrix_double.rows, shark::RealVector(true_descriptors));

		for (int y = 0; y < m_md_matrix_double.rows; ++y)
		{
			for (int ind = 0, x = 0; ind < weights.size(); ++ind)
			{
				if (true_descriptor_nums[ind] < m_md_matrix_double.cols)
				{
					temp[y][x] = m_md_matrix_double(y, true_descriptor_nums[ind]);
					++x;
				}
			}		
		}

		//evaluate with cross-validation
		shark::Data<shark::RealVector> input_data = shark::createDataFromRange(temp);
		shark::Data<unsigned int> input_labels = shark::createDataFromRange(unsigned_labels);
		auto labeled_data = shark::LabeledData<shark::RealVector, unsigned int>(input_data, input_labels);

		std::vector<double> param_min(1, 2.0), param_max(1, 2.0);
		std::vector<size_t> param_num(1, 1);
		auto svm_trainer_cross = GetSVMTrainerForLinearClassification<shark::RealVector>();

		svm_trainer_cross.SetGrid(param_min, param_max, param_num);

		double err = 0;
		const int kTrialsNum = 100;
		for (int trial_num = 0; trial_num < kTrialsNum; ++trial_num)
		{
			 err += svm_trainer_cross.train(labeled_data, 5);
		}
		err /= kTrialsNum;
		std::cout << err << "\n";

		std::ofstream f(m_mol_folder + m_mol_prefix + "_descr_set1.txt", std::ofstream::app);
		f << "======================================\n";
		f << m_command_line << std::endl;
			f << "descr = ";
			for (size_t descr_pos = 0; descr_pos < weights.size(); ++descr_pos)
			{
				f << true_descriptor_nums[descr_pos] << " ";
			}

			f << "precision = " << 1.0 - err;
			/*f << " recall = " << base_descr_set[ind].recall;
			f << " f_measure = " << base_descr_set[ind].FMeasure();*/
			f << "\n";
	}
}


void SetManager::ProcessMGUASVMRegression()
{
	//read MD matrix
	molecule_descriptor::ReadMatrix(m_md_matrix_double, m_mol_folder + m_mol_prefix + Extensions::PairsMDMatrix());
	//read labels
	std::vector<double> labels;
	molecule_descriptor::ReadVector(m_mol_folder + m_mol_prefix + Extensions::Properties(), labels);
	double mean = 0;
	std::vector<double> abs_labels(labels);
	for (auto it = labels.begin(), it_end = labels.end(); it != it_end; ++it)
	{
		abs_labels[it - labels.begin()] = abs(*it);
	}
	std::nth_element(abs_labels.begin(), abs_labels.begin() + abs_labels.size() / 2, abs_labels.end());
	mean = abs_labels[abs_labels.size() / 2];//labels.size();
	std::cout << "mean = " << mean << "\n";
	auto svm_trainer = GetSVMTrainerForRegression<shark::RealVector>(mean / 10.0);
	typedef decltype(svm_trainer) SVMType;
	const int kThreadsNum = 7;
	std::vector<SVMType> svm_vect(kThreadsNum);
	std::vector<double> param_min(2, -7.0), param_max(2, 3.0);
	std::vector<size_t> param_num(2, 13);
	param_max[1] = log(mean);
	param_min[1] = log(mean / 100.0);
	param_num[1] = 5;

	for (size_t ind = 0; ind < svm_vect.size(); ++ind)
	{
		svm_vect[ind] = GetSVMTrainerForRegression<shark::RealVector>(mean / 10.0);
		svm_vect[ind].SetGrid(param_min, param_max, param_num);
		svm_vect[ind].SetOutput(m_mol_folder + m_mol_prefix + "trainer_output1.txt");
	}

	MGUATrainer<shark::RealVector> mgua_trainer;
	/*cv::Mat_<double> md_mat(m_md_matrix.size());
	std::copy(m_md_matrix.begin(), m_md_matrix.end(), md_mat.begin());*/
	std::vector<shark::RealVector> labels_shark(labels.size(), shark::RealVector(1));
	for (size_t ind = 0; ind < labels_shark.size(); ++ind)
	{
		labels_shark[ind][0] = labels[ind];
	}
	mgua_trainer.SetData(m_md_matrix_double, labels_shark, false);
	mgua_trainer.SetParameters(20, 6);
	mgua_trainer.Train(svm_vect, m_mol_folder + m_mol_prefix, false);
}

std::string SetManager::MakeMoleculePrefix(const int ind)
{
	std::stringstream curr_file_prefix;
	curr_file_prefix << m_mol_folder << m_mol_prefix << "_" << ind;
	return curr_file_prefix.str();
}

void SetManager::CalculateDistThresholdsKMeans()
{
	const double null_elem = 0.0;
	const double dist_thresh = 0.1;
	const size_t max_iterations = 10000;
	std::vector<size_t> segmented;

	struct AbsMinus
	{
		double operator()(double arg1, double arg2) { 
			return abs(arg1 - arg2); }
	};

	std::function<double(double, double)> abs_minus_for_double = AbsMinus();

	SegmentKMeans(m_distances, m_distances_intervals, abs_minus_for_double, null_elem, 
		dist_thresh, max_iterations, 
		segmented, m_dist_thresholds);
	m_dist_thresholds.erase(m_dist_thresholds.begin());
	std::sort(m_dist_thresholds.begin(), m_dist_thresholds.end(), std::less<double>());
}



void SetManager::ReadParamsFromCommandLine(int argc, char** argv)
{
	ReadParamFromCommandLine(argc, argv, "-mol_folder", m_mol_folder);

	if (!m_mol_folder.empty())
	{
		if (m_mol_folder.back() != '\\')
		{
			m_mol_folder.push_back('\\');
		}
	}

	ReadParamFromCommandLine(argc, argv, "-mol_prefix", m_mol_prefix);
	ReadParamFromCommandLine(argc, argv, "-distances_intervals", m_distances_intervals);
	ReadParamFromCommandLineWithDefault(argc, argv, "-distances_intervals_levels", m_distance_interval_levels, 10);
	m_distance_quantile_for_thresh = m_distance_interval_levels - 1;
	ReadParamFromCommandLineWithDefault(argc, argv, "-mesh_levels_num", m_mesh_levels_num, 10);
	ReadParamFromCommandLineWithDefault(argc, argv, "-charges_threshes_num", m_charges_threshes_num, 2);
	ReadParamFromCommandLineWithDefault(argc, argv, "-lennard_jones_threshes_num", m_lennard_jones_threshes_num, 2);
	ReadParamFromCommandLineWithDefault(argc, argv, "-areas_threshes_num", m_areas_threshes_num, 0);
	ReadParamFromCommandLineWithDefault(argc, argv, "-use_levels", m_use_levels, true);
	ReadParamFromCommandLineWithDefault(argc, argv, "-use_triples", m_use_triples, false);

	ReadParamFromCommandLineWithDefault(argc, argv, "-levels_overlap", m_pairs_levels_overlap, 1);
	ReadParamFromCommandLineWithDefault(argc, argv, "-write_pairs", m_write_pairs, true);
	m_inited = true;
}

void SetManager::FindOutMoleculesNum()
{
	m_molecules_num = 0;
	SurfaceReader surface_reader;
	std::stringstream curr_surf_file;

	do
	{
		curr_surf_file.str("");
		curr_surf_file.clear();
		curr_surf_file << m_mol_folder << m_mol_prefix << "_" << m_molecules_num << Extensions::Surface();
		++m_molecules_num;		
	} 
	while (surface_reader.OpenFile(curr_surf_file.str()));

	m_molecules_num--;
}

}