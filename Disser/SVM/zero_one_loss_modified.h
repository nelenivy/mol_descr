//===========================================================================
/*!
 *  \brief Error measure for classication tasks, typically used for evaluation of results
 *
 *  \author T. Glasmachers
 *  \date 2010-2011
 *
 *
 *  <BR><HR>
 *  This file is part of Shark. This library is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 3, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library; if not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef SHARK_OBJECTIVEFUNCTIONS_LOSS_ZEROONELOSS_NON_AVERAGE_H
#define SHARK_OBJECTIVEFUNCTIONS_LOSS_ZEROONELOSS_NON_AVERAGE_H

#include "stdafx.h"
#include <shark/ObjectiveFunctions/Loss/AbstractLoss.h>

namespace shark {

///
/// \brief 0-1-loss for classification.
///
/// The ZeroOneLoss requires the existence of the comparison
/// operator == for its LabelType template parameter. The
/// loss function returns zero of the predictions exactly
/// matches the label, and one otherwise.
///
template<class LabelType = unsigned int, class OutputType = LabelType>
class ZeroOneLossNonAverage : public AbstractLoss<LabelType, LabelType>
{
public:
	typedef AbstractLoss<LabelType, LabelType> base_type;
	typedef typename base_type::BatchLabelType BatchLabelType;
	typedef typename base_type::BatchOutputType BatchOutputType;

	/// constructor
	ZeroOneLossNonAverage()
	{ }


	/// \brief From INameable: return the class name.
	std::string name() const
	{ return "ZeroOneLossNonAverage"; }

	using base_type::eval;
	virtual double eval(Data<LabelType> const& targets, Data<OutputType> const& predictions) const override{
		SIZE_CHECK(predictions.numberOfElements() == targets.numberOfElements());
		SIZE_CHECK(predictions.numberOfBatches() == targets.numberOfBatches());
		int numBatches = (int) targets.numberOfBatches();
		double error = 0;
		SHARK_PARALLEL_FOR(int i = 0; i < numBatches; ++i){
			double batchError= eval(targets.batch(i),predictions.batch(i));
			SHARK_CRITICAL_REGION{
				error+=batchError;
			}
		}
		return error;
	}
	///\brief Return zero if labels == predictions and one otherwise.
	double eval(BatchLabelType const& labels, BatchOutputType const& predictions) const{
		std::size_t numInputs = size(labels);
		SIZE_CHECK(numInputs == size(predictions));

		double error = 0;
		for(std::size_t i = 0; i != numInputs; ++i){
			error += (predictions(i) != labels(i))?1.0:0.0;
			//if (predictions(i) != labels(i))
				//std::cout << labels(i) << " ";
		}

		if (error > 0)
		{
			//std::cout << "\n";
		}

		return error;
	}
};


/// \brief 0-1-loss for classification.
template <>
class ZeroOneLossNonAverage<unsigned int, RealVector> : public AbstractLoss<unsigned int, RealVector>
{
public:
	typedef AbstractLoss<unsigned int, RealVector> base_type;
	typedef base_type::BatchLabelType BatchLabelType;
	typedef base_type::BatchOutputType BatchOutputType;

	/// constructor
	///
    /// \param threshold: in the case dim(predictions) == 1, predictions strictly larger than this parameter are regarded as belonging to the positive class
	ZeroOneLossNonAverage(double threshold = 0.0)
	{
		m_threshold = threshold;
	}

	/// \brief From INameable: return the class name.
	std::string name() const
	{ return "ZeroOneLossNonAverage"; }


	// annoyingness of C++ templates
	using base_type::eval;
	/// from AbstractCost
	///
	/// \param  targets      target values
	/// \param  predictions  predictions, typically made by a model
	virtual double eval(Data<unsigned int> const& targets, Data<RealVector> const& predictions) const{
		SIZE_CHECK(predictions.numberOfElements() == targets.numberOfElements());
		SIZE_CHECK(predictions.numberOfBatches() == targets.numberOfBatches());
		int numBatches = (int) targets.numberOfBatches();
		double error = 0;
		SHARK_PARALLEL_FOR(int i = 0; i < numBatches; ++i){
			double batchError= eval(targets.batch(i),predictions.batch(i));
			SHARK_CRITICAL_REGION{
				error+=batchError;
			}
		}
		return error;
	}
	/// Return zero if labels == arg max { predictions_i } and one otherwise,
	/// where the index i runs over the components of the predictions vector.
	/// A special version of dim(predictions) == 1 computes the predicted
	/// labels by thresholding at zero. Shark's label convention is used,
	/// saying that a positive value encodes class 0, a negative value
	/// encodes class 1.
	double eval(BatchLabelType const& labels, BatchOutputType const& predictions) const{
		std::size_t numInputs = size(labels);
		SIZE_CHECK(numInputs == (std::size_t)size(predictions));

		double error = 0;
		for(std::size_t i = 0; i != numInputs; ++i){
			const double res = evalSingle(labels(i),get(predictions,i));
			error+= res;
			//if (res > 0.0)
				//std::cout << i << " label " << labels(i) << " ";
		}
		if (error > 0.0)
		{
			//std::cout << "\n";
		}
		return error;
	}
private:
	template<class VectorType>
	double evalSingle(unsigned int label, VectorType const& predictions) const{
		std::size_t size = predictions.size();
		if (size == 1){
			// binary case, single real-valued predictions
			unsigned int t = (predictions(0) > m_threshold);
			if (t == label) return 0.0;
			else return 1.0;
		}
		else{
			// multi-class case, one prediction component per class
			RANGE_CHECK(label < size);
			double p = predictions(label);
			for (std::size_t i = 0; i<size; i++)
			{
				if (i == label) continue;
				if (predictions(i) >= p) return 1.0;
			}
			return 0.0;
		}
	}

	double m_threshold; ///< in the case dim(predictions) == 1, predictions strictly larger tha this parameter are regarded as belonging to the positive class
};


}
#endif
