#include "Main_includes.h"

void LLT_Factorization(
	vector < int > &ig,
	vector < int > &jg,
	vector < int > &ijg,
	vector < int > &idi,
	vector < double > &ggl,
	vector < double > &di,
	vector < int > &LLT_ig,
	vector < int > &LLT_jg,
	vector < int > &LLT_ijg,
	vector < int > &LLT_idi,
	vector < double > &LLT_ggl,
	vector < double > &LLT_di,
	int Nb
	);

void SLAE_Forward_Complex(
	vector < int > &LLT_ig,
	vector < int > &LLT_jg,
	vector < int > &LLT_ijg,
	vector < int > &LLT_idi,
	vector < double > &LLT_ggl,
	vector < double > &LLT_di,
	vector < double > &right_part,
	vector < double > &result,
	int Nb
	);

void SLAE_Backward_Complex(
	vector < int > &LLT_ig,
	vector < int > &LLT_jg,
	vector < int > &LLT_ijg,
	vector < int > &LLT_idi,
	vector < double > &LLT_ggl,
	vector < double > &LLT_di,
	vector < double > &right_part,
	vector < double > &result,
	int Nb
	);
void SqrtComplex(
	double *ab,
	double*xy
	);