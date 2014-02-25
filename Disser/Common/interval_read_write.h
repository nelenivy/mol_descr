#pragma once
#include <vector>
#include <utility>
#include <tuple>
#include <array>
#include <iterator>
#include <string>
#include <fstream>
#include <strstream>

#include "opencv2/core/core.hpp"

namespace molecule_descriptor
{

template <typename T>
bool ReadVector(const std::string& file_name, std::vector<T>& vect_out);
template <class IteratorType>
bool WriteInterval(const std::string& file_name, IteratorType begin, IteratorType end);
template <typename T>
void WriteMatrix(const cv::Mat_<T>& matrix, const std::string& filename);
template <typename T>
void ReadMatrix(const cv::Mat_<T>& matrix, const std::string& filename);
//////////////////////////////////////////////////////////////////////////
template <typename T>
void WriteMatrix(const cv::Mat_<T>& matrix, const std::string& filename)
{
	std::ofstream out_file(filename);
	CV_Assert(out_file.is_open());

	for (int y = 0; y < matrix.rows; y++)
	{
		for (int x = 0; x < matrix.cols; x++)
		{
			out_file << matrix(y, x) << " ";
		}

		out_file << "\n";
	}
}

template <typename T>
void ReadMatrix(cv::Mat_<T>& matrix, const std::string& filename)
{
	std::ifstream in_file(filename);
	CV_Assert(in_file.is_open());
	//calculate rows number
	std::string row;
	int rows_num = 0;
	for (; std::getline(in_file, row); ++rows_num)
	{}
	//calculate columns number
	in_file.close();
	in_file.open(filename);
	std::getline(in_file, row);
	std::stringstream row_stream(row);
	int columns = 0;
	T buf;

	for (; !row_stream.eof(); ++columns)
	{
		row_stream >> buf;
	}
	//
	in_file.close();
	in_file.open(filename);
	matrix.create(rows_num, columns);
	for (int y = 0; y < matrix.rows; y++)
	{
		for (int x = 0; x < matrix.cols; x++)
		{
			in_file >> matrix(y, x);
		}
	}
}
//function is used in ReadElemsVector function
template <typename T>
struct ReadElemNonChecked;
template<typename T1, typename T2>
struct ReadElemNonChecked<std::pair<T1, T2>>;

template <typename T>
struct ReadElemNonChecked
{
	static inline void Do(std::ifstream& file_in, T& new_elem)
	{
		file_in >> new_elem;
		char buf;
		file_in.get(buf);
	}
};

template <typename T>
struct ReadElemNonChecked<cv::Point3_<T>>
{
	static inline void Do(std::ifstream& file_in, cv::Point3_<T>& new_elem)
	{
		char buf;
		file_in.get(buf);
		file_in >> new_elem.x;
		file_in.get(buf);
		file_in >> new_elem.y;
		file_in.get(buf);
		file_in >> new_elem.z;
		file_in.get(buf);
	}
};

template <typename T>
struct ReadElemNonChecked<cv::Point_<T>>
{
	static inline void Do(std::ifstream& file_in, cv::Point_<T>& new_elem)
	{
		char buf;
		file_in.get(buf);
		file_in >> new_elem.x;
		file_in.get(buf);
		file_in >> new_elem.y;
		file_in.get(buf);
	}
};

//specialization for pair
template<typename T1, typename T2>
struct ReadElemNonChecked<std::pair<T1, T2>>
{
	static inline void Do(std::ifstream& file_in, std::pair<T1, T2>& new_elem)
	{
		ReadElemNonChecked<T1>::Do(file_in, new_elem.first);
		ReadElemNonChecked<T2>::Do(file_in, new_elem.second);
	}
};
//specialization for tuple
//BUG WITH PARTIAL SPECIALIZATION USING VARIADIC TEMPLATES
//template <int kIndex, typename...Types>
//inline void ReadTupleIndexNonChecked(std::ifstream& file_in, std::tuple<Types...>& new_elem)
//{
//	const size_t kCurrIndex = sizeof...(Types) - kIndex;
//	typedef typename std::tuple_element<kCurrIndex, std::tuple<Types...>> CurrentElemType;
//	ReadElemNonChecked<CurrentElemType>::Do(file_in, std::get<kCurrIndex>(new_elem));
//
//	if (kIndex > 1)
//	{
//		ReadTupleIndexNonChecked<kIndex - 1, Types...>(file_in, new_elem);
//	}
//}

template <typename T1, typename T2, typename T3>
struct ReadElemNonChecked<std::tuple<T1, T2, T3>>
{
	static inline void Do(std::ifstream& file_in, std::tuple<T1, T2, T3>& new_elem)
	{
		//ReadTupleIndexNonChecked<3, T1, T2, T3>(file_in, new_elem);
		ReadElemNonChecked<T1>::Do(file_in, std::get<0>(new_elem));
		ReadElemNonChecked<T2>::Do(file_in, std::get<1>(new_elem));
		ReadElemNonChecked<T3>::Do(file_in, std::get<2>(new_elem));
	}
};
template <typename T, size_t kArrSize>
struct ReadElemNonChecked<std::array<T, kArrSize>>
{
	static inline void Do(std::ifstream& file_in, std::array<T, kArrSize>& new_elem)
	{
		for (size_t ind = 0; ind < kArrSize; ++ind)
		{
			ReadElemNonChecked<T>::Do(file_in, new_elem[ind]);
		}
	}
};
//////////////////////////////////////////////////////////////////////////
template <typename T>
inline bool ReadElem(std::ifstream& file_in, T& new_elem, std::string& string_buf)
{
	bool res = false;
	ReadElemNonChecked<T>::Do(file_in, new_elem);

	if (file_in.rdstate() == 0)
	{
		res = true;
		//char buf;
		getline(file_in, string_buf);
		//file_in >> string_buf;
	}

	return res;
}

template <typename T>
bool ReadVector(const std::string& file_name, std::vector<T>& vect_out)
{
	std::ifstream file_in(file_name);

	if (!file_in.is_open())
	{
		return false;
	}

	std::string string_buf;
	//read vertices
	T new_elem;	

	for(vect_out.clear(); ReadElem(file_in, new_elem, string_buf); vect_out.push_back(new_elem))
	{	}

	file_in.clear();
	file_in.seekg(0, std::ios_base::beg);
	return !(vect_out.empty());
}
//////////////////////////////////////////////////////////////////////////
//used in WriteInterval function
template <typename T>
struct WriteElemToFile
{
	inline void operator()(std::ofstream& file_out, const T& elem)
	{
		file_out << elem << " ";
	}
};

template <typename T1, typename T2>
struct WriteElemToFile<std::pair<T1, T2>>;

template <typename T1, typename T2, typename T3>
struct WriteElemToFile<std::tuple<T1, T2, T3>>;
//specialization for pair
template <typename T1, typename T2>
struct WriteElemToFile<std::pair<T1, T2>>
{
	inline void operator()(std::ofstream& file_out, const std::pair<T1, T2>& elem)
	{
		WriteElemToFile<T1>()(file_out, elem.first);
		WriteElemToFile<T2>()(file_out, elem.second);
	}
};
//specialization for tuple
//function for writing tuple elements
//BUG WITH VARIADIC TEMPLATES
//template <int kIndex, typename T1, typename T2, typename T3>
//inline void WriteTupleIndexToFile(std::ofstream& file_out, const std::tuple<T1, T2, T3>& elem)
//{
//	const size_t kCurrIndex = 3 - kIndex;
//	typedef typename std::tuple_element<kCurrIndex, std::tuple<T1, T2, T3>> CurrentElemType;
//	WriteElemToFile<CurrentElemType>()(file_out, std::get<kCurrIndex>(elem));
//
//	if (kIndex > 1)
//	{
//		WriteTupleIndexToFile<kIndex - 1, T1, T2, T3>(file_out, elem)
//	}
//}

template <typename T1, typename T2, typename T3>
struct WriteElemToFile<std::tuple<T1, T2, T3>>
{
	inline void operator()(std::ofstream& file_out, const std::tuple<T1, T2, T3>& elem)
	{
		/*typedef typename std::tuple_element<kCurrIndex, std::tuple<Types...>> CurrentElemType;*/
		WriteElemToFile<T1>()(file_out, std::get<0>(elem));
		WriteElemToFile<T2>()(file_out, std::get<1>(elem));
		WriteElemToFile<T3>()(file_out, std::get<2>(elem));
		/*WriteTupleIndexToFile<3>(file_out, elem);*/
	}
};

template <typename T, size_t kArrSize>
struct WriteElemToFile<std::array<T, kArrSize>>
{
	inline void operator()(std::ofstream& file_out, const std::array<T, kArrSize>& elem)
	{
		for (size_t ind = 0; ind < kArrSize; ++ind)
		{
			WriteElemToFile<T>()(file_out, elem[ind]);
		}
	}
};

template <class IteratorType>
bool WriteInterval(const std::string& file_name, IteratorType begin, IteratorType end)
{
	std::ofstream file_out(file_name);

	if (!file_out.is_open())
	{
		return false;
	}

	//read vertices
	typedef typename std::iterator_traits<IteratorType>::value_type value_type;
	for ( ; begin != end; ++begin)
	{
		WriteElemToFile<value_type>()(file_out, *begin);
		file_out << ";\n";
	}
	
	return true;
}

}