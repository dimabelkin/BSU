#include "IO.h"

const std::string directory = "../resources/IN_data/";

void InputSparseComplexBlockSLAE(
	vector < double > &di                                          ,
	vector < double > &ggl                                         ,
	vector < int    > &ig                                          ,
	vector < int    > &jg                                          ,
	vector < int    > &idi                                         ,
	vector < int    > &ijg                                         ,
	vector < double > &righr_part                                  ,
	int               &FullTaskSize                                ,
	int               &BlockSize                                   ,
	double            &eps                                         ,
	int               &MaxIter                                     ,
	vector < double > &VectorForCheckMatrixVectorSLAEMultiplication
	)
{
	printf("Begin reading files...\n");
	std::string filename;
	FILE *fp;

	printf("\tReading kuslau...");
	filename = directory + "kuslau";
	fopen_s(&fp, filename.c_str(), "r");
	fscanf_s(fp, "%d %le %d", &FullTaskSize, &eps, &MaxIter);
	fclose(fp);
	printf("\r\tkuslau readed...\t");

	BlockSize = FullTaskSize / 2;

	printf("\r\tReading ig...");
	filename  = directory + "ig";
	fopen_s(&fp, filename.c_str(), "rb");
	ig.resize(BlockSize + 1);
	for (int i = 0; i < BlockSize + 1; ++i)
	{
		fread(&ig[i], sizeof(int), 1, fp);
		ig[i]--;
	}
	fclose(fp);
	printf("\r\tig readed...\t");

	printf("\r\tReading idi...");
	filename = directory + "idi";
	fopen_s(&fp, filename.c_str(), "rb");
	idi.resize(BlockSize + 1);
	for (int i = 0; i < BlockSize + 1; ++i)
	{
		fread(&idi[i], sizeof(int), 1, fp);
		idi[i]--;
	}
	fclose(fp);
	printf("\r\tidi readed...\t");

	printf("\r\tReading jg...");
	filename = directory + "jg";
	fopen_s(&fp, filename.c_str(), "rb");
	jg.resize(ig[BlockSize]);
	for (int i = 0; i < jg.size(); ++i)
	{
		fread(&jg[i], sizeof(int), 1, fp);
		jg[i]--;
	}
	fclose(fp);
	printf("\r\tjg readed...\t");

	printf("\r\tReading ijg...");
	filename = directory + "ijg";
	fopen_s(&fp, filename.c_str(), "rb");
	ijg.resize(ig[BlockSize]+1);
	for (int i = 0; i < ijg.size(); ++i)
	{
		fread(&ijg[i], sizeof(int), 1, fp);
		ijg[i]--;
	}
	fclose(fp);
	printf("\r\tijg readed...\t");

	printf("\r\tReading gg...");
	filename = directory + "gg";
	fopen_s(&fp, filename.c_str(), "rb");
	ggl.resize(ijg[ig[BlockSize]]);
	for (int i = 0; i < ggl.size(); ++i)
	{
		fread(&ggl[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	printf("\r\tgg readed...\t");

	printf("\r\tReading di...");
	filename = directory + "di";
	fopen_s(&fp, filename.c_str(), "rb");
	di.resize(idi[BlockSize]);
	for (int i = 0; i < di.size(); ++i)
	{
		fread(&di[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	printf("\r\tdi readed...\t");

	printf("\r\tReading pr...");
	filename = directory + "pr";
	fopen_s(&fp, filename.c_str(), "rb");
	righr_part.resize(FullTaskSize);
	
	for (int i = 0; i < righr_part.size(); ++i)
	{
		fread(&righr_part[i], sizeof(double), 1, fp);
	}
	fclose(fp);
	printf("\r\tpr readed...\t");

	////VectorForCheckMatrixVectorSLAEMultiplication
	//sprintf(filename, "%s/y.txt", dirname);
	//VectorForCheckMatrixVectorSLAEMultiplication.resize(FullTaskSize);
	//fp = fopen(filename, "rb");
	//for (int i = 0; i < VectorForCheckMatrixVectorSLAEMultiplication.size(); i++)
	//{
	//	fread(&VectorForCheckMatrixVectorSLAEMultiplication[i], sizeof(double), 1, fp);
	//}
	//fclose(fp);

	printf("\rAll files have been read.\t\t\n\n");
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