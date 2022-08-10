// MyCpp.cpp: implementation of the CMyCpp class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "stdio.h"

#include <stdarg.h>
#include <stdlib.h>
#include "math.h"
#include "public.h"


//////////////////////////////////////////////////////////////////////////
//   显示提示信息
void MyBreak(char* fmt, ...)
{
	char buffer[256];
	va_list argptr;
	va_start(argptr, fmt);
	vsprintf(buffer, fmt, argptr);
	va_end(argptr);
	
#ifdef VC_EXTRALEAN
	AfxMessageBox(buffer);
#else
	printf(buffer);
	getchar();
#endif // VC_EXTRALEAN
	
}



//////////////////////////////////////////////////////////////////////////
//       对称矩阵下标计算函数
int ij(int i,int j)
{
	return (i>=j)? i*(i+1)/2+j :j*(j+1)/2+i;
}


//////////////////////////////////////////////////////////////////////////
//  权倒数计算函数
double Calculate_q(double *Q,double *F,int t)
{
	double q=0.0;
	for(int i=0;i<t;i++)
		for(int j=0;j<t;j++) 
			q+=Q[ij(i,j)]*F[i]*F[j];
		
	return q;
		
}



//////////////////////////////////////////////////////////////////////////
//  权倒数计算函数(系数系列进行了压缩)
double Calculate_q(double Q[],double F[],int Fin[],int n)
{
	double q=0.0;
	for(int k=0;k<n;k++)
	{
		int i=Fin[k];
		for(int s=0;s<n;s++) 
		{
			int j=Fin[s];
			q+=Q[ij(i,j)]*F[k]*F[s];
		}
	}
	return q;
	
}




//////////////////////////////////////////////////////////////////////////
//  将度分秒连写的(double型)角度化为弧度值
double dms_rad(double a)
{
	//提取角度值的符号
	double sign=(a<0.0) ? -1.0 : 1.0;
	a=fabs(a);
	
	//提取角度值的整度
	int d=(int)((a+0.00001)/10000.0);
	a=a-d*10000.0;
	if(a<0.0){ d=d-1; a=a+10000;}
	
	//提取角度值的整分及秒值
	int m=(int)((a+0.00001)/100.0);
	a=a-m*100;
	if(a<0.0){ m=m-1; a=a+100.0;}
	
	a=sign*(d*3600.0+m*60.0+a)/206264.806247096363;
	
	return a;
}


//////////////////////////////////////////////////////////////////////////
//  将角度的弧度值化为度分秒连写的角度（double 型） 
double rad_dms(double a)
{
	a=a*206264.806247096363;
	
	double sign=(a<0.0) ? -1.0 : 1.0;
	a=fabs(a);
	
	int d=(int)(a/3600.0+0.0000001);
	a=a-d*3600.0;
	
	if(a<0.0){ d=d-1; a=a+3600.0; }
	
	int m=(int)(a/60.0+0.0001);
	a=a-m*60.0;
	if(a<0.0){ m=m-1; a=a+60.0; }
	
	a=d*10000.0+m*100.0+a;
	
	return a*sign;
}



//////////////////////////////////////////////////////////////////////////
//    权逆阵传播计算
void Calculate_BQBT(double B[],double Q[],int r,int n,double N[])
{
	for(int i=0;i<r;i++)
		for(int j=0;j<=i;j++)
		{
			double nij=0.0;
			for(int k=0;k<n;k++)
				for(int s=0;s<n;s++)
					nij+=B[i*n+k]*Q[ij(k,s)]*B[j*n+s];
				N[ij(i,j)]+=nij;
		}
}




//////////////////////////////////////////////////////////////////////////
//    向文件输出数组
void PrintM(FILE *fp,double A[],int size, int t,char* fmt,
			char* title,bool IsLabel)
{
	if(title)fprintf(fp,"\n %s: ",title);
	int j=0;
	for(int i=0;i<size;i++)
	{
		if(i%t==0)
		{
			j++;
			if(IsLabel)fprintf(fp,"\n%3d ",j);
			else fprintf(fp,"\n    ");
		}
		fprintf(fp,fmt,A[i]);
	}	
	fprintf(fp,"\n");
}



//////////////////////////////////////////////////////////////////////////
//  向文件输出对称矩阵（下三角元素）
void  PrintM2(FILE* fp, double M[], int n, int t,char *fmt,
			  char* title,bool IsLabel)
{
	if(title)fprintf(fp,"\n %s: ",title);
	
	int index=0;
	for(int i=0;i<n;i++)
	{
		if(IsLabel)fprintf(fp,"\n%3d ",i+1);
		else fprintf(fp,"\n    ");
		for(int j=0;j<=i;j++)
		{
			if(j>0 && j%t==0)fprintf(fp,"\n    ");
			fprintf(fp,fmt,M[index++]);
		}
	}
	fprintf(fp,"\n");	
}


//////////////////////////////////////////////////////////////////////////
//  向文件输出线性方程组
void  PrintEquation(FILE* fp, double A[], double b[],
					int n, int t, char *fmt, char* title)
{
	if(title)fprintf(fp,"\n %s: ",title);
	for(int i=0;i<n;i++)
	{
		fprintf(fp,"\n%3d ",i+1);
		for(int j=0;j<t;j++)
		{
			fprintf(fp,fmt,A[i*t+j]);
		}
		fprintf(fp,fmt,b[i]);
	}
}


//////////////////////////////////////////////////////////////////////////
//  对称正定矩阵求逆(仅存下三角元素)
bool inverse(double a[],int n)
{
    double *a0=new double[n];
    for(int k=0;k<n;k++)
	{ 
		double a00=a[0];
		if(a00+1.0==1.0)
		{
			delete []a0; 
			return false;
		}
		for(int i=1;i<n;i++)
		{
			double ai0 = a[i*(i+1)/2];
			
			if(i<=n-k-1)a0[i]= -ai0/a00;
			else        a0[i]=  ai0/a00;

			for(int j=1;j<=i;j++)
			{
				a[(i-1)*i/2+j-1]=a[i*(i+1)/2+j]+ai0*a0[j];
			} 
		}
		for(i=1;i<n;i++)
		{
			a[(n-1)*n/2+i-1]=a0[i];
		}
		a[n*(n+1)/2-1]=1.0/a00;
	} 
	delete []a0;
	return true;
}


#define IGG1  0
#define IGG3  1
#define Huber 2
double Wi(int fname, double v, double k0, double k1)
{
	double a;
	switch(fname)
	{
	case IGG1://IGG1函数
		v=fabs(v);
		if(v<=k0)return 1.0;
		if(v>k1)return 0.0;
		return k0/v;
	case IGG3://IGG3函数
		v=fabs(v);
		if(v<=k0)return 1.0;
		if(v>k1)return 0.0;
		a=(k1-v)/(k1-k0);
		return k0/v*a*a;
	case Huber: // Huber函数
		v=fabs(v);
		if(v<=k0)return 1.0;
		return k0/v;
	default:
		MyBreak("等价权函数名称错误！");
	return 1.0;
	}
}



//////////////////////////////////////////////////////////////////////////
// 中位数计算：将数组由大至小排序，返回值为中位数
double Median(double pp[],int n,bool IsAbs)
{
	double *p=new double [n];
	if(IsAbs)
	{
		for(int i=0;i<n;i++)p[i]=fabs(pp[i]);
	}
	else
	{
		for(int i=0;i<n;i++)p[i]=pp[i];
	}
	
	int k=n/2;
	while(k>0)
	{  
		for(int j=k;j<=n-1;j++)
		{  
			double t=p[j]; 
			int i=j-k;
			
			while( (i>=0) && (p[i]>t) )
			{
				p[i+k]=p[i];	
				i=i-k;  
			}
			p[i+k]=t;
		}
		k=k/2;
	}
	
	double mean = (n%2==1)? p[n/2] : (p[n/2]+p[n/2-1])/2.0;
	delete []p;
	
	return  mean;
	
}

