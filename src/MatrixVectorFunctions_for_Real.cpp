#include "MatrixVectorFunctions_for_Real.h"

double RealScalarProduct(
	vector < double > &vec1,
	vector < double > &vec2	
	)
{
	double result = 0;
	for (int i = 0; i < vec1.size(); i++)
	{
		result += vec1[i] * vec2[i];
	}
	return result;
}

void RealMultVector_by_Scalar(
	vector < double > &vec,
	double scalar,
	vector < double > &result	
	)
{
	for (int i = 0; i < vec.size(); i++)
	{
		result[i] = vec[i] * scalar;
	}
}