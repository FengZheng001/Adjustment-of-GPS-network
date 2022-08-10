// CoGPS.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "MyGPS.h"
#include "public.h"

//////////////////////////////////////////////////////////////////////////
//    最小二乘平差算例
void  main1()
{
	CCoGPS coGPS;
	char *rfile="教材例题\\最小二乘平差\\r.txt";
	char *dfile ="教材例题\\最小二乘平差\\Data.txt";
	coGPS.resultfp=fopen(rfile,"w");
	if(coGPS.resultfp==NULL)
	{
		MyBreak("打开结果文件失败");
		return;
	}	
	fprintf(coGPS.resultfp," GPS网最小二乘平差算例:\n");
	
	coGPS.InputData(dfile);
	//coGPS.PrintData();
	
	coGPS.LeastSquares();
	fclose(coGPS.resultfp);
	
	printf("\n\n\n数据文件：%s",dfile);
	printf("\n结果文件：%s\n\n",rfile);
	
}

//////////////////////////////////////////////////////////////////////////
//  抗差估计算例
void main2()
{
	CCoGPS coGPS;
	char *rfile="教材例题\\robust\\r.txt";
	char *dfile ="教材例题\\robust\\Data.txt";
	
	coGPS.resultfp=fopen(rfile,"w");
	if(coGPS.resultfp==NULL)
	{
		MyBreak("打开结果文件失败!");
		return;
	}	
	
	coGPS.InputData(dfile); //输入向量，内含输入坐标
	// coGPS.PrintData();
	
	coGPS.Robust(IGG3,0.005);
	fclose(coGPS.resultfp);
	
	printf("\n\n\n数据文件：%s",dfile);
	printf("\n结果文件：%s\n\n",rfile);
	
}
//////////////////////////////////////////////////////////////////////////
//  粗差探测算例
void main3()
{
	CCoGPS coGPS;
	char *dfile="教材例题\\粗差探测\\Data.txt";
	char *rfile="教材例题\\粗差探测\\r.txt";
	
	coGPS.resultfp=fopen(rfile,"w");
	if(coGPS.resultfp==NULL)
	{
		MyBreak("打开结果文件失败");
		return;
	}	
	
	coGPS.InputData(dfile); //输入向量，内含输入坐标
	//coGPS.PrintData();
	
	coGPS.DataSnooping(0.001,20);
	fclose(coGPS.resultfp);
	
	printf("\n\n\n数据文件：%s",dfile);
	printf("\n结果文件：%s\n\n",rfile);
}

//////////////////////////////////////////////////////////////////////////
//    自由网平差算例
void main4()
{
	CCoGPS coGPS;	
	char *dfile="教材例题\\自由网平差\\Data.txt";
	char *rfile="教材例题\\自由网平差\\r.txt";
	
	coGPS.resultfp=fopen(rfile,"w");
	if(coGPS.resultfp==NULL)
	{
		MyBreak("打开结果文件失败");
		return;
	}	
	
	coGPS.InputData(dfile); //输入向量，内含输入坐标
	//coGPS.PrintData();
	
	coGPS.Free();
	fclose(coGPS.resultfp);			
	
	printf("\n\n\n数据文件：%s",dfile);
	printf("\n结果文件：%s\n\n",rfile);
	
}


//////////////////////////////////////////////////////////////////////////
//     拟稳平差算例
void main5()
{
	CCoGPS coGPS;	
	char *dfile="教材例题\\拟稳平差\\Data.txt";
	char *rfile="教材例题\\拟稳平差\\r.txt";
	coGPS.resultfp=fopen(rfile,"w");
	if(coGPS.resultfp==NULL)
	{
		MyBreak("打开结果文件失败");
		return;
	}	
	
	coGPS.InputData(dfile); //输入向量，内含输入坐标
	//coGPS.PrintData();
	
	coGPS.Quasi_Stable("教材例题\\拟稳平差\\拟稳点.txt");
	fclose(coGPS.resultfp);
	
	printf("\n\n\n数据文件：%s",dfile);
	printf("\n结果文件：%s\n\n",rfile);
	
}

//////////////////////////////////////////////////////////////////////////
//    已知点坐标加权平差算例
void  main6()
{
	CCoGPS coGPS;
	char *dfile="教材例题\\已知点加权平差\\Data.txt";
	char *rfile="教材例题\\已知点加权平差\\r.txt";
	coGPS.resultfp=fopen(rfile,"w");
	if(coGPS.resultfp==NULL)
	{
		MyBreak("打开结果文件失败");
		return;
	}
	
	coGPS.InputData(dfile);
	//coGPS.PrintData();
				
	char *XYZfile="教材例题\\已知点加权平差\\已知XYZ.txt";
	coGPS.InputKnownXYZ(XYZfile,false);
	
	coGPS.LeastSquares();
	fclose(coGPS.resultfp);
	printf("\n结果文件：%s\n\n",rfile);
	
}

//////////////////////////////////////////////////////////////////////////
//  已知点坐标加相关权平差算例
void main7()
{
	CCoGPS coGPS;	
	char *dfile="教材例题\\已知点加相关权平差\\Data.txt";
	char *rfile="教材例题\\已知点加相关权平差\\r.txt";
	char *XYZfile="教材例题\\已知点加相关权平差\\已知XYZ.txt";
	
	coGPS.resultfp=fopen(rfile,"w");
	if(coGPS.resultfp==NULL)
	{
		MyBreak("打开结果文件失败");
	}	
	
	coGPS.InputData(dfile); //输入向量，内含输入坐标
	//coGPS.PrintData();
	
	coGPS.InputKnownXYZ(XYZfile,true); //true表示相关权
	coGPS.LeastSquares();
	
	fclose(coGPS.resultfp);
	printf("\n结果文件：%s\n\n",rfile);
}


//////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{

	printf("\n  ==== GPS网平差计算 ====");
	printf("\n  1.最小二乘平差");
	printf("\n  2.最小二乘平差（参数加权）");
	printf("\n  3.最小二乘平差（参数加相关权）");
	
	printf("\n  4.自由网平差");
	printf("\n  5.拟稳平差\n");
	printf("\n  6.抗差估计");
	printf("\n  7.粗差探测\n");
	
	int sw=getch()-'0';
	
	switch(sw) {
	case 1: // 最小二乘平差
		{
			main1();
		}
		break;
	case 2:  //最小二乘平差（参数加权）
		{
			main6();
		}
		break;
	case 3: //最小二乘平差（参数加相关权）
		{
			main7();
		}
		break;
	case 4:
		{   // 自由网平差
			main4();
		}
		break;
	case 5:
		{   // 拟稳平差
			main5();
		}
		break;
	case 6: // 抗差估计
		{
			main2();
		}
		break;
	case 7:  // 粗差探测
		{  
			main3();			
		}
		break;
		
	default:
		break;
	}
	return 0;
}
