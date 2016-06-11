#include "Main_includes.h"

void mult_block(
	double *x,
	double *y,
	double *a,
	int size
	);
void mult_MV(
	vector < int > &ig,
	vector < int > &jg,
	vector < double > &ggl,
	vector < double > &di,
	vector < int > &ijg,
	vector < int > &idi,
	vector < double > &x,
	vector < double > &y,
	int nb
	);
void mult_complex_numbers(
	double *x,
	double *y,
	double *res
	);
void div_complex_numbers(
	double *x,
	double *y,
	double *res
	);
void di_preconditioning(
	double * x,
	double *res
	);
void conjugate_SK_mult(
	vector < double > &x,
	vector < double > &y,
	double  *res,
	int nb
	);
void norm(
	vector < double > x,
	double &res,
	int nb
	);
void sub_vec(
	vector < double > &x,
	vector < double > &y,
	vector < double > &res
	);
void sum_vec(
	vector < double > &x,
	vector < double > &y,
	vector < double > &res
	);
/*void copy_vec(
	double *x,
	double *y,
	int size
	);*/
void mult_vector_by_scalar(
	vector < double > &x,
	double *skalar,
	vector < double > &res,
	int nb
	);
/*void mult_DI_matrix_vector(
	vector < double > &inverse_di,
	vector < double > &r,
	vector < double > &res,
	int nb
	)
	*/

void sum_complex_numbers(
	double *a,
	double *b,
	double *res
	);

void sub_complex_numbers(
	double *a,
	double *b,
	double *res
	);
void Complex_SK_mult(
	vector < double > &x,
	vector < double > &y,
	double *res,
	int nb
	);
void MultDiOnVect(
	vector < double > &di,
	vector < double > &vec,
	vector < double > &res,
	int Nb
	);
void MultLTMatrixOnVector(
	vector < int > &LLT_ig,
	vector < int > &LLT_jg,
	vector < int > &LLT_ijg,
	vector < int > &LLT_idi,
	vector < double > &LLT_ggl,
	vector < double > &LLT_di,
	vector < double > &vec,
	vector < double > &result,
	int Nb
	);
void MultVMatrixOnVector(
	vector < vector < double > > &V,
	vector < double > vec,
	vector < double > &res,
	int Nb,
	int m
	);
void NormalizeVec(
	vector < double > &r,
	double norma,
	int Nb
	);

void copy_vec(
	double *x,
	double *y,
	int size
	);
