#include "stdafx.h"
#include "mgua.h"

namespace molecule_descriptor
{



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

		f << "precision = " << base_descr_set[ind].precision;
		f << " recall = " << base_descr_set[ind].recall;
		f << " f_measure = " << base_descr_set[ind].FMeasure();
		f << "\n";
	}
}
}