#pragma once

namespace molecule_descriptor
{
//mark curvature type
const unsigned char kConvexType = 1;
const unsigned char kConcaveType = 2;
const unsigned char kSaddleType = 3;

struct IntRange
{
	int min_val;
	int max_val;
	int step;

	int ElemsInRange() const
	{
		return (max_val - min_val + 1) / step;
	}
};

}