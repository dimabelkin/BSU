#include "Main_includes.h"
#include "IO.h"
#include "COCG_for_Complex.h"
#include "GMRES_Block.h"
#include <iostream>

double CheckMatrixVectorMultiplication(
	
	)
{
	return 0;
}

void main()
{
	vector < double > Cdi, Cggl, Cright_part, Cresult;
	vector < int > Cidi, Cig, Cijg, Cjg;
	vector < double > VectorForCheckMatrixVectorSLAEMultiplication;
	int MaxIter, FullTaskSize, BlockSize;
	double eps;

	InputSparseComplexBlockSLAE("../resources/IN_data", Cdi, Cggl, Cig, Cjg, Cidi, Cijg, Cright_part, FullTaskSize, BlockSize, eps,
		MaxIter, VectorForCheckMatrixVectorSLAEMultiplication);

	//inputSparceBlockMatrix_from_txt("../resources/IN_data/testSLAE2.txt", Cdi, Cggl, Cig, Cjg,
	//	Cidi, Cijg, Cright_part, FullTaskSize, BlockSize, eps,
	//	MaxIter);

	Cresult.resize(FullTaskSize);

	COCG(Cig, Cjg, Cggl, Cdi, Cijg, Cidi, Cright_part, BlockSize, Cresult,eps, MaxIter);

	//GMRES_for_complex(Cig, Cjg, Cggl, Cdi,
	//	Cijg, Cidi, Cright_part,
	//	BlockSize, Cresult, eps, MaxIter);

	std::cout << "end" << std::endl;
}