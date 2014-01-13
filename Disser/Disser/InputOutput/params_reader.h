#pragma  once
#include <sstream>
#include <string>

namespace molecule_descriptor
{

template <typename T>
bool ReadParamFromCommandLine(int argc, char** argv, const std::string& param_key, T& param)
{
	for (int ind = 1; ind < argc; ++ind)
	{
		if (std::string(argv[ind]) == param_key)
		{
			if (ind == argc - 1)
			{
				return false;
			}

			std::stringstream string_param;
			string_param << argv[ind + 1];
			string_param >> param;
			return true;
		}
	}

	return false;
}

}