#pragma once
#include "stdafx.h"
#include <shark/Algorithms/DirectSearch/GridSearch.h>

namespace molecule_descriptor
{

class GridSearchFromEnd : public shark::GridSearch
{
public:
	//! Please note that for the grid search optimizer it does
	//! not make sense to call step more than once, as the
	//! solution does not improve iteratively.
	GridSearchFromEnd() {
		m_out = &std::cout;
	}

	void SetOutput(std::ostream& out)
	{
		m_out = &out;
	}
	void step(const ObjectiveFunctionType& objectiveFunction) {
		const size_t dimensions = m_nodeValues.size();
		std::vector<size_t> index(dimensions, 0);
		m_best.value = 1e100;
		shark::RealVector point(dimensions);

		// loop through all grid points
		while (true)
		{
			// define the parameters
			for (size_t dimension = 0; dimension < dimensions; dimension++)
				point(dimension) = m_nodeValues[dimension][index[dimension]];

			// evaluate the model
			if (objectiveFunction.isFeasible(point))
			{
				double error = objectiveFunction.eval(point);

//#ifdef SHARK_CV_VERBOSE
				(*m_out) << point << "\t" << error << std::endl;
				std::cout << point << "\t" << error << std::endl;

//#endif	    
				if (error < m_best.value)
				{
					m_best.value = error;
					m_best.point = point;		// [TG] swap() solution is out, caused ugly memory bug, I changed this back
				}
			}

			{// next index
				int dimension = static_cast<int>(dimensions) - 1;
				for (; dimension >= 0; dimension--)
				{
					index[dimension]++;
					if (index[dimension] < m_nodeValues[dimension].size()) 
						break;
					index[dimension] = 0;
				}
				if (dimension == - 1) break;
			}

			for (int dimension = static_cast<int>(dimensions) - 1; dimension >= 0; dimension--)
			{
				std::cout << index[dimension] << " ";
			}
			std::cout << "\n";

		}
	}

private:
	std::ostream* m_out;
};

}