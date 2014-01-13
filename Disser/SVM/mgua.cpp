#include "mgua.h"

namespace molecule_descriptor
{

void MGUATrainer::SetParameters(size_t max_iterations, size_t descr_to_select_num)
{
	m_max_iterations = max_iterations;
	m_descr_to_select_num = descr_to_select_num;

	//change descriptors set containers
	m_base_descr_set.resize(m_descr_to_select_num);
	m_curr_best_set.resize(m_descr_to_select_num);

	std::fill(m_base_descr_set.begin(), m_base_descr_set.end(), DescriptorsSetWithError(m_max_iterations, kBigError));
	std::fill(m_curr_best_set.begin(), m_curr_best_set.end(), DescriptorsSetWithError(m_max_iterations, kBigError));
}

void WriteSetToFile(const std::vector<DescriptorsSetWithError>& base_descr_set, const std::string& file_name)
{
	std::ofstream f(file_name, std::ofstream::app);
	f << "======================================\n";

	for (int ind = 0; ind < base_descr_set.size(); ++ind)
	{
		f << "descr = ";
		const auto& curr_descr_set = base_descr_set[ind].descriptors;

		for (int descr_pos = 0; descr_pos < curr_descr_set.size(); ++descr_pos)
		{
			f << curr_descr_set[descr_pos] << " ";
		}

		f << "err = " << base_descr_set[ind].error;
	}
}
}