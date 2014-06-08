#pragma once 
#include "boost/archive/polymorphic_text_oarchive.hpp"
#include "../Common/singular_point.h"
#include "element_with_index_and_id.h"

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

		template <class Archive, typename PropT>
		void serialize(Archive& ar, molecule_descriptor::SingularPointsPair<PropT>& point, const unsigned int version)
		{
			ar & point.elem1;
			ar & point.elem2;
			ar & point.dist;
		}
	}
}
