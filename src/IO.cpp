#include "IO.h"

void InputSparseComplexBlockSLAE(
	char* dirname,
	vector < double > &di,
	vector < double > &ggl,
	vector < int > &ig,
	vector < int > &jg,
	vector < int > &idi,
	vector < int > &ijg,
	vector < double > &righr_part,
	int &FullTaskSize,
	int &BlockSize,
	double &eps,
	int &MaxIter,
	vector < double > &VectorForCheckMatrixVectorSLAEMultiplication
	)
{
	char filename[500];
	//kuslau
	sprintf(filename, "%s/kuslau", dirname);
	FILE *fp;
	fp = fopen(filename, "r");
	fscanf(fp, "%d%le%d", &FullTaskSize, &eps, &MaxIter);
	fclose(fp);

	BlockSize = FullTaskSize / 2;

	//ig
	ig.resize(BlockSize + 1);
	sprintf(filename, "%s/ig", dirname);
	fp = fopen(filename, "rb");
	for (int i = 0; i < BlockSize + 1; i++)
	{
		fread(&ig[i], sizeof(int), 1, fp);
		ig[i]--;
	}
	fclose(fp);

	//idi
	sprintf(filename, "%s/idi", dirname);
	idi.resize(BlockSize + 1);
	fp = fopen(filename, "rb");
	for (int i = 0; i < BlockSize + 1; i++)
	{
		fread(&idi[i], sizeof(int), 1, fp);
		idi[i]--;
	}
	fclose(fp);

	//jg
	sprintf(filename, "%s/jg", dirname);
	jg.resize(ig[BlockSize]);
	fp = fopen(filename, "rb");
	for (int i = 0; i < jg.size(); i++)
	{
		fread(&jg[i], sizeof(int), 1, fp);
		jg[i]--;
	}
	fclose(fp);
	//ijg
	sprintf(filename, "%s/ijg", dirname);
	ijg.resize(ig[BlockSize]+1);
	fp = fopen(filename, "rb");
	for (int i = 0; i < ijg.size(); i++)
	{
		fread(&ijg[i], sizeof(int), 1, fp);
		ijg[i]--;
	}
	fclose(fp);
	//ggl
	sprintf(filename, "%s/gg", dirname);
	ggl.resize(ijg[ig[BlockSize]]);
	fp = fopen(filename, "rb");
	for (int i = 0; i < ggl.size(); i++)
	{
		fread(&ggl[i], sizeof(double), 1, fp);	
	}
	fclose(fp);
	//di
	sprintf(filename, "%s/di", dirname);
	di.resize(idi[BlockSize]);
	fp = fopen(filename, "rb");
	for (int i = 0; i < di.size(); i++)
	{
		fread(&di[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	//right part
	sprintf(filename, "%s/pr", dirname);
	righr_part.resize(FullTaskSize);
	fp = fopen(filename, "rb");
	for (int i = 0; i < righr_part.size(); i++)
	{
		fread(&righr_part[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	//VectorForCheckMatrixVectorSLAEMultiplication
	sprintf(filename, "%s/y.txt", dirname);
	VectorForCheckMatrixVectorSLAEMultiplication.resize(FullTaskSize);
	fp = fopen(filename, "rb");
	for (int i = 0; i < VectorForCheckMatrixVectorSLAEMultiplication.size(); i++)
	{
		fread(&VectorForCheckMatrixVectorSLAEMultiplication[i], sizeof(double), 1, fp);
	}
	fclose(fp);
}

	void inputSparceBlockMatrix_from_txt(
		char* filename,
		vector < double > &di,
		vector < double > &ggl,
		vector < int > &ig,
		vector < int > &jg,
		vector < int > &idi,
		vector < int > &ijg,
		vector < double > &righr_part,
		int &FullTaskSize,
		int &BlockSize,
		double &eps,
		int &MaxIter
		)
	{
		
		
		FILE *fp;
		fp = fopen(filename, "r");
		fscanf(fp, "%d%le%d", &FullTaskSize, &eps, &MaxIter);
		BlockSize = FullTaskSize / 2;
		//ig
		ig.resize(BlockSize + 1);	
		for (int i = 0; i < BlockSize + 1; i++)
		{
			fscanf(fp, "%d", &ig[i]);
			ig[i]--;
		}
		//idi		
		idi.resize(BlockSize + 1);	
		for (int i = 0; i < BlockSize + 1; i++)
		{
			fscanf(fp, "%d", &idi[i]);
			idi[i]--;
		}
		//jg	
		jg.resize(ig[BlockSize]);	
		for (int i = 0; i < jg.size(); i++)
		{
			fscanf(fp, "%d", &jg[i]);
			jg[i]--;
		}	
		//ijg		
		ijg.resize(ig[BlockSize] + 1);		
		for (int i = 0; i < ijg.size(); i++)
		{
			fscanf(fp, "%d", &ijg[i]);
			ijg[i]--;
		}	
		//ggl		
		ggl.resize(ijg[ig[BlockSize]]);		
		for (int i = 0; i < ggl.size(); i++)
		{
			fscanf(fp, "%le", &ggl[i]);
		}		
		//di		
		di.resize(idi[BlockSize]);		
		for (int i = 0; i < di.size(); i++)
		{
			fscanf(fp, "%le",&di[i]);
		}		
		//right part		
		righr_part.resize(FullTaskSize);		
		for (int i = 0; i < righr_part.size(); i++)
		{
			fscanf(fp, "%le", &righr_part[i]);
		}
		fclose(fp);
	}

////////////////////////////////////////////////////////////////////////////////////////////////////
	void OutputItersResidual(
		double Residual, 
		int IterNumber,
		char *filename
		)
	{
		FILE *fp;
		if (IterNumber == 0)
		{
			fp = fopen(filename, "w");
		}
		else
		{
			fp = fopen(filename, "a");
		}
		fprintf(fp, "%d\t%.15le\n", IterNumber, Residual);
		fclose(fp);
	}