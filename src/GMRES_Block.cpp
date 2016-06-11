#include "GMRES_Block.h"


void BuildVandHMatrix(
	vector < vector < double > > &VMatrix,
	vector < vector < double > >&HMatrixUpTr,
	vector < double > &HMatrixLowDiag,
	vector < double > &r,
	int Nb,
	int &m,
	int FactorisationFlag,
	vector < int > &LLT_ig,
	vector < int > &LLT_jg,
	vector < int > &LLT_ijg,
	vector < int > &LLT_idi,
	vector < double > &LLT_ggl,
	vector < double > &LLT_di,
	vector < int > &ig,
	vector < int > &jg,
	vector < int > &ijg,
	vector < int > &idi,
	vector < double > &ggl,
	vector < double > &di,
	vector < double > &InverseDI, 
	int iterNum

	)
{
	double norma = 0;
	norm(r, norma, Nb);
	VMatrix[0] = r;
	NormalizeVec(VMatrix[0], norma, Nb);
	vector < double > tmp, tmp2;
	tmp.resize(2 * Nb);
	vector < double > w;
	w.resize(2 * Nb);
	tmp2.resize(2 * Nb);
	//resize H matrix
	if (iterNum == 0)
	{
		HMatrixLowDiag.resize(m);
		HMatrixUpTr.resize(m);
		for (int i = 0; i < m; i++)
		{
			HMatrixUpTr[i].resize(2*(m - i));
		}
	}
	for (int i = 1; i <= m; i++)
	{
		for (int j = 0; j < 2 * Nb; j++)
		{
			VMatrix[i][j] = 0;
		}
	}
	for (int i = 0; i < m; i++)
	{
		HMatrixLowDiag[i] = 0;
		for (int j = 0; j < 2 * (m - i); j++)
		{
			HMatrixUpTr[i][j] = 0;
		}
	}
	for (int mu = 0; mu < m; mu++)
	{
		for (int i = 0; i < 2 * Nb; i++)
		{
			tmp[i] = 0;
			tmp2[i] = 0;
		}

		if (FactorisationFlag == 1)//diag factorisation
		{
			MultDiOnVect(InverseDI, VMatrix[mu], tmp, Nb);
			mult_MV(ig, jg, ggl, di, ijg, idi, tmp, tmp2, Nb);
			MultDiOnVect(InverseDI, tmp2, w, Nb);
		}
		else  //LLT
		{
			mult_MV(ig, jg, ggl, di, ijg,
				idi, VMatrix[mu], tmp, Nb);
			SLAE_Forward_Complex(LLT_ig, LLT_jg, LLT_ijg, 
				LLT_idi, LLT_ggl, LLT_di, tmp, tmp, Nb);
			SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg,
				LLT_idi, LLT_ggl, LLT_di, tmp, w, Nb);
		}
		//VMatrix.push_back(vector < double >(2 * Nb));
		for (int lambda = 0; lambda <= mu; lambda++)
		{
			double H[2];
			Complex_SK_mult(VMatrix[lambda], w, H, Nb);
			HMatrixUpTr[lambda][2 * (mu - lambda)] = H[0];
			HMatrixUpTr[lambda][2 * (mu - lambda) + 1] = H[1];
			mult_vector_by_scalar(VMatrix[lambda], H,tmp , Nb);//right
			sum_vec(VMatrix[mu + 1], tmp, VMatrix[mu + 1]);
			
		}
		sub_vec(w, VMatrix[mu+1], VMatrix[mu+1]);//right
		norm(VMatrix[mu+1], HMatrixLowDiag[mu], Nb);
		if (abs(HMatrixLowDiag[mu]) < 1e-10)
		{
			m = mu+1;
			return;
		}
		NormalizeVec(VMatrix[mu+1], HMatrixLowDiag[mu], Nb);
	}
}
void RotateHMatrix(
	vector < vector < double > > &HMatrixUpTr,
	vector < double > &HMatrixLowDiag,
	vector < double > &right_part,
	int MSize
	)
{
	
	vector < vector < double > > TmpHMatrixUpTr;
	vector < double > TmpHMatrixLowDiag;
	vector < double > tmpRight_part;
	vector < double > StartHMatrixLowDiag;
	StartHMatrixLowDiag = HMatrixLowDiag;
	double c[2];
	double c_conjugate[2];
	double s;
	double norma;
	double tmp[2];
	for (int i = 0; i < MSize; i++)
	{
		TmpHMatrixUpTr = HMatrixUpTr;
		TmpHMatrixLowDiag = HMatrixLowDiag;
		tmpRight_part = right_part;
		norma = 0;
		norm(HMatrixUpTr[i], norma, HMatrixUpTr[i].size()/2);
		norma = norma*norma;
		norma += StartHMatrixLowDiag[i] * StartHMatrixLowDiag[i];
		norma = sqrt(norma);
		s = StartHMatrixLowDiag[i] / norma;
		c[0] = HMatrixUpTr[i][0] / norma;
		c[1] = HMatrixUpTr[i][1] / norma;
		c_conjugate[0] = c[0];
		c_conjugate[1] = -c[1];
		tmp[0] = 0;
		tmp[1] = 0;
		mult_block(&HMatrixUpTr[i][0], tmp, c_conjugate, 2);
		TmpHMatrixUpTr[i][0] = tmp[0];
		TmpHMatrixUpTr[i][1] = tmp[1];
		TmpHMatrixUpTr[i][0] += s*HMatrixLowDiag[i];
		TmpHMatrixLowDiag[i] = (-s*HMatrixUpTr[i][0] + c[0]*HMatrixLowDiag[i]);//должен получаться 0
		//комплексная часть сокращается
		//проверка
		double CheckComplex = -s*HMatrixUpTr[i][1] + c[1] * HMatrixLowDiag[i];

		for (int j = 1; j < MSize - i; j++)
		{
			tmp[0] = 0;
			tmp[1] = 0;
			mult_block(&HMatrixUpTr[i][2 * j], tmp, c_conjugate, 2);
			TmpHMatrixUpTr[i][2 * j] = tmp[0] + s*HMatrixUpTr[i+1][2*(j-1)];
			TmpHMatrixUpTr[i][2 * j + 1] = tmp[1] + s*HMatrixUpTr[i + 1][2 * (j - 1)+1];

			tmp[0] = 0;
			tmp[1] = 0;
			mult_block(&HMatrixUpTr[i+1][2 * (j-1)], tmp, c, 2);
			TmpHMatrixUpTr[i+1][2 * (j-1)] = tmp[0] + s*HMatrixUpTr[i ][2 * (j)];
			TmpHMatrixUpTr[i+1][2 *(j -1) + 1] = tmp[1] + s*HMatrixUpTr[i ][2 * (j ) + 1];

		}
		tmpRight_part[2 * i] = 0;
		tmpRight_part[2 * i + 1] = 0;
		mult_block(c_conjugate, &tmpRight_part[2 * i], &right_part[2 * i], 2);
		tmpRight_part[2 * i] += s*right_part[2 * (i + 1)];
		tmpRight_part[2 * i + 1] += s*right_part[2 * (i + 1) + 1];

		tmpRight_part[2 * (i + 1)] = 0;
		tmpRight_part[2 * (i + 1)+1] = 0;
		mult_block(c, &tmpRight_part[2 * (i + 1)], &right_part[2 * (i + 1)], 2);
		tmpRight_part[2 * (i + 1)] += -s*right_part[2 * (i)];
		tmpRight_part[2 * (i + 1) + 1] += -s*right_part[2 * (i)+1];

		HMatrixUpTr = TmpHMatrixUpTr;
		HMatrixLowDiag = TmpHMatrixLowDiag;
		right_part = tmpRight_part;
	}
}

void SolveHMatrixSLAE(
	vector < vector < double > > &HMatrixRotated,
	vector < double > &right_part,
	vector < double > &result,
	int MSize
	)
{
	result = right_part;
	double tmp[2];
	for (int i = MSize - 1; i >= 0; i--)
	{
		
		for (int j = 1; j < MSize - i - 1; j++)
		{
			tmp[0] = 0;
			tmp[1] = 0;
			mult_block(&HMatrixRotated[i][2 * j], &result[2 * (i + j)], tmp, 2);
			sub_complex_numbers(&result[2 * (i + j)], tmp, &result[2 * (i + j)]);
		}
		div_complex_numbers(&result[2 * i], &HMatrixRotated[i][0], &result[2 * i]);
	}
}

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
	)
{
	vector < double > LLT_ggl, LLT_di;
	vector < int > LLT_ig, LLT_jg, LLT_ijg, LLT_idi;
	LLT_Factorization(ig, jg, ijg, idi, ggl, di, LLT_ig, 
		LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, nb);
	vector < vector < double > > VMatrix;
	vector < vector < double > > HMatrixUpTr;
	vector < double > HMatrixLowDiag;
	vector < double > InverseDiagSQRTFactorisation, DiagSQRTFactorisation;
	vector < double > x, r,  tmp, tmp2, d, z;
	int m = 50;
	int FactorizationType = 2;//тип факторизации 1- диагональная, 2 - LLT
	double complex_number[2];
	InverseDiagSQRTFactorisation.resize(2 * nb);
	DiagSQRTFactorisation.resize(2 * nb);
	tmp.resize(2 * nb);
	tmp2.resize(2 * nb);
	x.resize(2 * nb);
	r.resize(2 * nb);
	for (int i = 0; i < nb; i++)
	{
		int size = idi[i + 1] - idi[i];
		
		if (size == 2)
		{
			//di_preconditioning(&di[idi[i]], &InverseDiagSQRTFactorisation[2 * i]);
			DiagSQRTFactorisation[2 * i] = di[idi[i]]; //корень считается верно
			DiagSQRTFactorisation[2 * i+1] = di[idi[i] +1];
		}
		else
		{
			/*InverseDiagSQRTFactorisation[2 * i] = 1.0 / di[idi[i]];
			InverseDiagSQRTFactorisation[2 * i + 1] = 0;*/
			DiagSQRTFactorisation[2 * i] = di[idi[i]];
			DiagSQRTFactorisation[2 * i + 1] = 0;
		}
		SqrtComplex(&DiagSQRTFactorisation[2 * i], complex_number);
		DiagSQRTFactorisation[2 * i] = complex_number[0];
		DiagSQRTFactorisation[2 * i + 1] = complex_number[1];
		di_preconditioning(&DiagSQRTFactorisation[2*i], &InverseDiagSQRTFactorisation[2 * i]);//предобуславливание диагональное верно
		
	/*	InverseDiagSQRTFactorisation[2 * i] = complex_number[0];
		InverseDiagSQRTFactorisation[2 * i + 1] = complex_number[0];*/
	}
	double nev;
	result.resize(2 * nb);
	
	mult_MV(ig, jg, ggl, di, ijg, idi, result, tmp, nb);
	sub_vec(right_part, tmp, tmp);
	if (FactorizationType == 1)
	{
		MultDiOnVect(DiagSQRTFactorisation, result, x, nb);
		MultDiOnVect(InverseDiagSQRTFactorisation, tmp, r, nb);//умножение на диагональ верно
		
	}
	else
	{
		MultLTMatrixOnVector(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, result, x, nb);
		SLAE_Forward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, tmp, r, nb);		
	}
	double norm0;
	norm(r, nev, nb);
	norm0 = nev;
	int iterNumber = 0;
	VMatrix.resize(m + 1);
	for (int i = 0; i < m+1; i++)
	{
		VMatrix[i].resize(2 * nb);
	}
	printf("%le\r", nev / norm0);
	OutputItersResidual(nev / norm0, iterNumber, "../resources/OUT_data/nev1.txt");
	//OutputItersResidual(nevSecond, iter + 1, "../resources/OUT_data/nev2.txt");
	//до этого момента для диагонального предобуславливания все верно
	double t_start = omp_get_wtime();
	while (nev / norm0 >= eps && iterNumber < MaxIter)
	{
		for (int i = 0; i < 2 * nb; i++)
		{
			tmp[i] = 0;
			tmp2[i] = 0;
		}
		BuildVandHMatrix(VMatrix, HMatrixUpTr, HMatrixLowDiag, r, nb, m, FactorizationType, 
			LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, 
			ig, jg, ijg, idi, ggl, di, InverseDiagSQRTFactorisation, iterNumber);
		d.resize(2 * (m+1));
		z.resize(2 * m);
		for (int i = 0; i < 2 * (m + 1); i++)
		{
			d[i] = 0;
		}
		d[0] = nev;
		RotateHMatrix(HMatrixUpTr, HMatrixLowDiag, d, m);//right?
	//	d.resize(2 * m);
		SolveHMatrixSLAE(HMatrixUpTr, d, z, m);//right
		MultVMatrixOnVector(VMatrix, z, tmp, nb, m);//right
		sum_vec(x, tmp, x);
		if (FactorizationType == 1)
		{
			MultDiOnVect(InverseDiagSQRTFactorisation, x, tmp, nb);
			mult_MV(ig, jg, ggl, di, ijg, idi, tmp, tmp2, nb);
			sub_vec(right_part, tmp2, tmp);
			MultDiOnVect(InverseDiagSQRTFactorisation, tmp, r, nb);			
		}
		else
		{
			
			SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, x, tmp, nb);
			mult_MV(ig, jg, ggl, di, ijg, idi, tmp, tmp2, nb);
			sub_vec(right_part, tmp2, tmp2);
			SLAE_Forward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, tmp2, r, nb);
		}
		norm(r, nev, nb);
		printf("%le\r", nev / norm0);
		iterNumber++;
		OutputItersResidual(nev / norm0, iterNumber, "../resources/OUT_data/nev1.txt");
		//OutputItersResidual(nevSecond, iter + 1, "../resources/OUT_data/nev2.txt");
		
	
	}
	if (FactorizationType == 1)
	{
		MultDiOnVect(InverseDiagSQRTFactorisation, x, result, nb);
	}
	else
	{
		SLAE_Backward_Complex(LLT_ig, LLT_jg, LLT_ijg, LLT_idi, LLT_ggl, LLT_di, x, result, nb);
	}
	//calc true nev
	mult_MV(ig, jg, ggl, di, ijg, idi, result, tmp, nb);
	sub_vec(right_part, tmp, tmp);
	double true_nev;
	norm(tmp, true_nev, nb);
	true_nev /= norm0;
	double t_end = omp_get_wtime();
	FILE *fp;
	fp = fopen("Times.txt", "a");
	fprintf(fp, "%d\t%.3le\t%d\t%.3le\n", m, t_end - t_start, iterNumber, true_nev);
	fclose(fp);
	printf("\n\nTrue Nev = %le\n", true_nev );
	OutputItersResidual(true_nev, iterNumber, "../resources/OUT_data/nev1.txt");
	//OutputItersResidual(nevSecond, iter + 1, "../resources/OUT_data/nev2.txt");
}