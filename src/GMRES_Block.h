#include "Main_includes.h"
#include "LLT.h"
#include "MatrixVectorFunctions_for_Complex.h"
#include "MatrixVectorFunctions_for_Real.h"
#include "IO.h"
void GMRES_for_complex(
	vector < int > &ig,
	vector < int > &jg,
	vector < double > &ggl,
	vector < double > &di,
	vector < int > &ijg,
	vector < int > &idi,
	vector < double > &right_part,
	int nb,
	vector < double > &result,
	double eps,
	int MaxIter
	);