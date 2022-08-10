// MyGPS.h: interface for the CMyGPS class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYGPS_H__C7B8ED1B_4EB2_46BF_B706_63D12F6B803D__INCLUDED_)
#define AFX_MYGPS_H__C7B8ED1B_4EB2_46BF_B706_63D12F6B803D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <conio.h>
#include <stdlib.h>


#include "stdio.h"
#include "math.h"
#include "string.h"

class CCoGPS  
{
public:
    CCoGPS();
	virtual ~CCoGPS();
	FILE   *resultfp; //输出文件指针
	bool  m_IsPx;  // 是否参数加权平差	
	
	int    Anumber;      //同步区数，向量组总数
    int    Pnumber;      //总点数
    int    Vnumber;      //向量总数
    int    unPnumber;    //未知点数

    double *XYZ;     // 坐标数组	
    char   **Pname;  // 点名指针数组
	bool   *IsKnown;//已知点标志数组
	
    int    *dir0;        //各同步区首向量在观测值中的序号
    int    *dir1,*dir2;  //向量端点号

	double *L;  //全区观测向量
	double *ATPA,*ATPL;
	double MEDIAN;  //由中位数计算的单位权中误差
	double Sigma;   //由VTPV计算的单位权中误差
	double **P;     //权矩阵数组；
	double *QV;     //观测值残差权倒数
	double *V;      //残差
    double *PV;     //
	double *W;      //权因子	

    double *MxMyMz;  //已知点的中误差向量
	double *KnownXYZ;  //已知点坐标,参数加权平差时用
	int    *KnownPoints; //与已知点坐标相应的点号
    double *KnownPx;  //已知点权矩阵
	int    *AreaNumber;  //各向量所属向量组的编号	

public:
	bool  InputData(char *DataFile); //输入原始数据
	bool  InputKnownXYZ(char *file,bool corelated); //输入已知点坐标
    void  PrintData(); //打印原始数据
		
	void  PrintXYZ();  //打印平差坐标
    void  PrintLV(char *Title,double *LV);//打印观测值、残差

	void  LeastSquares(); //最小二乘平差计算
	void  Free();         //秩亏自由网平差计算
	void  Quasi_Stable(char *file);//拟稳平差计算
	void  Robust(int fname,double eps); // 抗差最小二乘平差	
	void  DataSnooping(double Arfa,int MaxLoop);  // 粗差探测计算

	void CCoGPS::Helmert();
		
private:
	void  CaATPA();	// 组成法方程
    void  CaATPAi(int n,int *dir1,int *dir2,double *Pi,double *L);  
	void  Known();//已知点处理函数
	void   CaLV(double V[]);//计算残差或自由项向量
	double CaVTPV(); //计算残差二次型
    double Ca_dX(bool IsPrint); //法方程解算,返回最大坐标改正数
	int  GetStationNumber(char *buf); //保存点名，返回点号
	double Ca_VTPV_XYZ();//已知点加相关权平差时坐标改正数二次型
	
	// 以下是粗差探测计算辅助函数
	void CaJTPJ(double JTPJ[], int Jint[],int s);
	double CaQg(int s,int Jint[],double gross[],
		         double Qg[],double JTPA[],double JTPV[] );
	void CaJTPA(double JTPA[],int s,int n0);
	void Ca_X_Qx(double gross[],double Qg[],double JTPA[],int s);

	//以下是抗差估计辅助函数
    void   Ca_UintV(); 
	void   caQV();
};





#endif // !defined(AFX_MYGPS_H__C7B8ED1B_4EB2_46BF_B706_63D12F6B803D__INCLUDED_)





















