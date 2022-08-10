// MyCpp.h: interface for the CMyCpp class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MYCPP_H__C63EB901_D617_4F46_96FE_5B8B28A71665__INCLUDED_)
#define AFX_MYCPP_H__C63EB901_D617_4F46_96FE_5B8B28A71665__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>


//////////////////////////////////////////////////////////////////////////
//       对称矩阵下标计算函数

#define ik (i>=k)? i*(i+1)/2+k : k*(k+1)/2+i
#define kj (k>=j)? k*(k+1)/2+j : j*(j+1)/2+k




//  权倒数计算函数
double Calculate_q(double *Q,double *F,int t);
//  权倒数计算函数
double Calculate_q(double Q[],double F[],int Fin[],int n);




void MyBreak(char* fmt, ...);

int ij(int i,int j);



bool inverse(double a[],int n);

/*
bool inverse3(double A[],int n);
bool inverse4(double A[],int n);
*/


// 权逆阵传播计算
void Calculate_BQBT(double B[],double Q[],int r,int n,double N[]);



//  将度分秒连写的(double型)角度化为弧度值
double dms_rad(double a);

//  将角度的弧度值化为度分秒连写的角度（double 型） 
double rad_dms(double a);


void PrintM(FILE *fp,double A[],int size, int t,char* fmt,
			char* title=NULL,bool IsLabel=true);

void  PrintM2(FILE* fp, double M[], int n,int t,char *fmt,
			  char* title=NULL,bool IsLabel=true);

void  PrintEquation(FILE* fp, double A[], double b[],
					int n, int m,char *fmt,char* title=NULL);


#define IGG1  0
#define IGG3  1
#define Huber 2

double Wi(int fname, double v, double k0, double k1);

double Median(double pp[],int n,bool IsAbs);//中位数计算


#endif // !defined(AFX_MYCPP_H__C63EB901_D617_4F46_96FE_5B8B28A71665__INCLUDED_)
