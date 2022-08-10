// MyGPS.cpp: implementation of the CMyGPS class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "MyGPS.h"
#include "public.h"


#define Old_Method  1
#define New_Method  0
#define Method  Old_Method


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CCoGPS::CCoGPS()
{
	m_IsPx = false; //不进行参数加权平差
}


CCoGPS::~CCoGPS()
{
  
}



/////////////////////////////////////////////////////////////
//        一个同步区的法方程,用于最小二乘平差
void CCoGPS::CaATPAi(int nk,  //观测值总数，不是向量总数
					 int *dir1,
					 int *dir2,
					 double *Pk,
					 double *Lk // 自由项
					)  
{
     for(int s1=0;s1<nk;s1++)
     {
         int i1=3*dir1[s1/3]+s1%3; //误差方程Vk中系数为-1的未知数编号
         int i2=3*dir2[s1/3]+s1%3; //误差方程Vk中系数为+1的未知数编号

         for(int s2=0; s2<nk; s2++)
         {
            int j1=3*dir1[s2/3]+s2%3;
            int j2=3*dir2[s2/3]+s2%3;

			double p12 = Pk[ij(s1,s2)];

            double l=Lk[s2];

            ATPL[i1]-=p12*l;
            ATPL[i2]+=p12*l;

            if(i1>=j1)ATPA[ij(i1,j1)]+=p12;
            if(i1>=j2)ATPA[ij(i1,j2)]-=p12;
            if(i2>=j1)ATPA[ij(i2,j1)]-=p12;
            if(i2>=j2)ATPA[ij(i2,j2)]+=p12;
         }
     }
	 
}




//////////////////////////////////////////////////////////////////////////
//    组成法方程式
void CCoGPS::CaATPA()
{
     int t=3*Pnumber;
     int tt = t*(t+1)/2;
     for(int i=0;i<t;i++) ATPL[i]=0.0;
     for(i=0;i<tt;i++) ATPA[i]=0.0;

     double *l=new double [3*Vnumber];//误差方程自由项
	 CaLV(l); //误差方程自由项计算

     for(i=0;i<Anumber;i++) //逐向量组组成法方程
     {
         int s=dir0[i];
         int n=dir0[i+1]-s;
         CaATPAi(3*n,dir1+s,dir2+s,P[i],l+3*s);
     }

     delete []l;
}






//////////////////////////////////////////////////////////////////////
//      输入原始数据（近似坐标、GPS向量）
bool CCoGPS::InputData(char *file)
{
	 FILE *fp=fopen(file,"r"); // 打开原始数据文件
	 if(fp==NULL)
	 {
		 MyBreak("原始数据文件打不开!");
		 return 0;
	 }

	 fscanf(fp,"%d",&Pnumber); //总点数
	 fscanf(fp,"%d",&Anumber); //同步区总数（向量组总数）
	 fscanf(fp,"%d",&Vnumber); //向量总数
	 
	 Pname = new char* [Pnumber]; //点名地址数据
	 IsKnown = new bool[Pnumber]; //已知/未知点标志数组
	 for(int i=0;i<Pnumber;i++)Pname[i]=NULL;
	 
     dir0=new int [Anumber+1];//各同步区首向量在总观测向量中的序号
     dir1=new int [Vnumber];  //向量起始点号
     dir2=new int [Vnumber];  //向量终点点号
	 AreaNumber=new int [Vnumber]; //各向量所属向量组的编号
	 
     int n=3*Vnumber;  //观测值总数
     L=new double [n]; //向量观测值
     V=new double [n]; //残差
 
	 P=new double* [Anumber]; //观测值权矩阵地址
     PV=new double [n]; //权矩阵于残差向量的乘积
	 
     int t=3*Pnumber;
     int tt=t*(t+1)/2;
	 XYZ =new double[t];   //坐标数组
     ATPA=new double [tt]; //法方程系数矩阵（或逆矩阵）
     ATPL=new double [t];  //法方程自由项向量

	 unPnumber=0; //未知点总数
	 for(int k=0;k<Pnumber;k++) //输入坐标
	 {
		 char name[20];
		 fscanf(fp,"%s",name); //点名
		 int i=GetStationNumber(name); //将点名换算成点号
		 if(i<0)
		 {
			 MyBreak("读取坐标是出错：点名多于总点数");  
			 return 0;
		 }
		 
		 fscanf(fp,"%lf%lf%lf",XYZ+3*i,XYZ+3*i+1,XYZ+3*i+2);//坐标值 
		 fscanf(fp,"%d",IsKnown+i); //已知点标志
		 if(!IsKnown[i])unPnumber++;
	 }

     dir0[0]=0;
     for( i=0;i<=Anumber-1;i++)   //逐同步区循环
     {
		 int d;
		 fscanf(fp,"%d",&d); //同步区编号		 
		 if(d!=i+1)
		 {
			 MyBreak("同步区编号不连续");
			 return 0;
		 }

         int jd;  //同步区内向量数
		 fscanf(fp,"%d",&jd);

		 dir0[i+1]=dir0[i]+jd;
         for(int j=dir0[i];j<dir0[i+1];j++)
         {			 
             char P1[20],P2[20];
			 fscanf(fp,"%s %s",P1,P2); //输入向量的端点点名
             dir1[j]=GetStationNumber(P1);
             dir2[j]=GetStationNumber(P2);
			 if(dir1[j]<0 || dir2[j]<0)
			 {
				 MyBreak("读取向量数据出错：点名多余总点数");
				 return 0;
			 }

             //输入向量观测值
			 fscanf(fp,"%lf %lf %lf",L+3*j,L+3*j+1,L+3*j+2);
			 AreaNumber[j]=i;
         }

		 //---------------------------------------
		 //     输入方差协方差矩阵
		 int ni=3*jd;
         double *Pi=new double [(ni+1)*ni/2];
		 P[i]=Pi;		 
		 
	 	 int s=0;
         for(j=0;j<ni;j++)
		 {
			 int d;
			 fscanf(fp,"%d",&d); //行号
		     for(int k=0;k<=j;k++)
			 {
				 fscanf(fp,"%lf",Pi+s); //向量的方差协方差
				 s++;
			 }
		 }
     }
     fclose(fp);

	 return 1;
}


//////////////////////////////////////////////////////////////////////
//       输出原始数据
void CCoGPS::PrintData()
{
	fprintf(resultfp,"\n  总点数: %d",Pnumber);
	fprintf(resultfp,"\n  向量组总数: %d",Anumber);
	fprintf(resultfp,"\n  向量总数: %d\n",Vnumber);

	//输出已知坐标和近似坐标	
	for(int k=0; k<Pnumber; k++) 
	{
		fprintf(resultfp,"\n %4s",Pname[k]); //输出点名
		fprintf(resultfp,"  %14.5lf",XYZ[3*k]);
		fprintf(resultfp,"  %14.5lf",XYZ[3*k+1]);
		fprintf(resultfp,"  %14.5lf",XYZ[3*k+2]);
		fprintf(resultfp,"  %d",IsKnown[k]);
	}		
	 
	//输出向量观测值及其方差阵
	for(int i=0; i<Anumber; i++)
	{
		 int jd=dir0[i+1]-dir0[i]; //组内向量总数
         fprintf(resultfp,"\n\n%4d  %d",i+1,jd);//向量组编号，组内向量数

         for(int j=dir0[i];j<dir0[i+1];j++)
         {
             int k1=dir1[j];
             int k2=dir2[j];

             double dx=L[3*j];
             double dy=L[3*j+1];
             double dz=L[3*j+2];
	         
			 fprintf(resultfp,"\n %8s %8s",Pname[k1],Pname[k2]);
			 fprintf(resultfp,"%17.5lf %17.5lf %17.5lf ",dx,dy,dz);
		 }

         int jj=0;
         for(j=0;j<3*jd;j++)
         {
			  fprintf(resultfp,"\n%3d",j+1);
              for(int k=0;k<=j;k++)
			  {
				   if(k%4==0 && k>0)fprintf(resultfp,"\n   "); 
				   fprintf(resultfp," %17.9e",P[i][jj++]);
			  }
		 }
	 }
}


//////////////////////////////////////////////////////////////////////////
//    未知数计算
double CCoGPS::Ca_dX(bool IsPrint)
{
     int t=3*Pnumber;
     if(!inverse(ATPA,t))
	 {
		 MyBreak("法方程不满秩!");
		 fclose(resultfp);
		 exit(0);
	 }
     
	 if(IsPrint)
	 {
		 fprintf(resultfp,"\n\n           ===== 坐标改正数 =====\n");		 
	 }
	 double max=0.0;
	 for(int i=0;i<t;i++)
	 {
          double xi=0.0;
   	      for(int j=0;j<t;j++)xi-=ATPA[ij(i,j)]*ATPL[j];
		  XYZ[i]+=xi;
		  if(fabs(xi)>fabs(max))max=xi;
		  if(IsPrint)
		  {
			  if(i%3==0) fprintf(resultfp,"\n%10s ",Pname[i/3]);
			  fprintf(resultfp," %8.4lf",xi);
		  }
     }
	 return max;
}


//////////////////////////////////////////////////////////////////////
//   计算误差方程自由项或最小二乘残差
void CCoGPS::CaLV(double V[])
{
     for(int i=0;i<Vnumber;i++)
     {
         int k1=dir1[i];
         int k2=dir2[i];
         
		 V[3*i]=XYZ[3*k2]-XYZ[3*k1]-L[3*i];         
         V[3*i+1]=XYZ[3*k2+1]-XYZ[3*k1+1]-L[3*i+1]; 
         V[3*i+2]=XYZ[3*k2+2]-XYZ[3*k1+2]-L[3*i+2]; 
     }
}

//////////////////////////////////////////////////////////////////////
//    计算 PV,VTPV
double CCoGPS::CaVTPV()
{
	 double VTPV=0.0;
	 for(int i=0;i<Anumber;i++) // Anumber：向量组总数
	 {
	    double *Pi=P[i];     //i区权矩阵起始地址
        int k1=3*dir0[i];    //i区首向量在观测值中的序号
        int k2=3*dir0[i+1];  //i+1区首向量在观测值中的序号

		double *PVi=PV+k1;
		double *Vi=V+k1;
		int n=k2-k1;

		double pvv=0.0;
		for(int j=0;j<n;j++)
		{
			double pvj=0.0;
			for(int k=0; k<n; k++) pvj+=Pi[ij(j,k)]*Vi[k];
			PVi[j]=pvj;
			pvv += pvj*Vi[j];
		}
		VTPV+=pvv;
	 }
	 return VTPV;
}



//////////////////////////////////////////////////////////////////////
//    输出残差V,也可以输出自由项向量，还可以输出权因子
void  CCoGPS::PrintLV(char *Title,double *V)
{
	fprintf(resultfp,"\n\n        	         ===== %s =====\n\n",Title);
	fprintf(resultfp,"    向量组  点 名       点 名"
		"     Vx         Vy         Vz");

	int last_ni = -1;
	for(int i=0;i<Vnumber;i++)
	{
		int ni=AreaNumber[i];
		if(ni!=last_ni)	fprintf(resultfp,"\n");
		last_ni=ni; 
		
		fprintf(resultfp,"\n%8d",ni+1);
		fprintf(resultfp,"%10s %10s",Pname[dir1[i]],Pname[dir2[i]]);
		fprintf(resultfp,"%10.4lf %10.4lf %10.4lf",V[3*i],V[3*i+1],V[3*i+2]);
	}  
}




//////////////////////////////////////////////////////////////////////
//  已知点处理（已知点固定，已知点加权）
void CCoGPS::Known()
{
	if(!m_IsPx) //已知坐标固定平差
	{
		for(int i=0;i<Pnumber;i++)
		{
			if(IsKnown[i])
			{
				ATPA[ij(3*i,3*i)]+=1.0e30;
				ATPA[ij(3*i+1,3*i+1)]+=1.0e30;
				ATPA[ij(3*i+2,3*i+2)]+=1.0e30;
			}
		}
		return;
	}

    //--------------------------------------
	//  已知坐标加权平差
	int kn=Pnumber-unPnumber;
	for(int i=0;i<3*kn;i++) //系数矩阵
	{
		int p1=3*KnownPoints[i/3]+i%3;
		for(int j=0;j<=i;j++)
		{
			int p2=3*KnownPoints[j/3]+j%3;
			ATPA[ij(p1,p2)]+=KnownPx[ij(i,j)];
		} 
	} 
		 
	//法方程自由项处理
	double *dX=new double[3*kn]; //已知坐标与近似坐标的坐标差
	for(i=0;i<kn;i++)
	{
		int k=KnownPoints[i];
		dX[3*i]  =XYZ[3*k]  -KnownXYZ[3*i];
		dX[3*i+1]=XYZ[3*k+1]-KnownXYZ[3*i+1];
		dX[3*i+2]=XYZ[3*k+2]-KnownXYZ[3*i+2];
	}

	for(i=0;i<3*kn;i++)
	{
		double dL=0.0;
		for(int j=0;j<3*kn;j++)
			dL+=KnownPx[ij(i,j)]*dX[j];
		int k=KnownPoints[i/3]; //点号
		ATPL[3*k+i%3]+=dL;
	} 
	delete []dX;
	return;
}



//////////////////////////////////////////////////////////////////////
//   自由网平差函数
void CCoGPS::Free()
{
	//将观测值的方差协方差矩阵求逆，化为权矩阵
	for(int i=0; i<Anumber; i++)
	{
		int ni = 3*(dir0[i+1]-dir0[i]); //向量组观测值总数
		inverse(P[i],ni);
	}
	
	CaATPA();     	 
	for( i=0; i<3*Pnumber; i++)
	{
		for(int j=0; j<=i;j++)
		{
			if(i%3==j%3) ATPA[ij(i,j)]+=1.0/Pnumber;
		}
	}
	
	double max=Ca_dX(true);     
	for( i=0; i<3*Pnumber; i++)
	{
		for(int j=0; j<=i;j++)
		{
			if(i%3==j%3) ATPA[ij(i,j)]-=1.0/Pnumber;
		}
	}
	
	CaLV(V);  //计算残差
	
	double VTPV=CaVTPV();
	Sigma=sqrt(VTPV/(3*Vnumber-3*(Pnumber-1)));

	PrintXYZ();
	PrintLV("V",V);
	fprintf(resultfp,"\n\n自由网平差:\n[pvv]=%e",VTPV);
	fprintf(resultfp,"\n单位权中误差μ=±%e",Sigma);	
	
}




//////////////////////////////////////////////////////////////////////
//  最小二乘平差函数
void CCoGPS::LeastSquares()
{
	//将观测值的方差协方差矩阵求逆，化为权矩阵
	for(int i=0; i<Anumber; i++)
	{
		int ni = 3*(dir0[i+1]-dir0[i]); //向量组观测值总数
		inverse(P[i],ni);
	}
	
	CaATPA();     	 
	Known();
	double max=Ca_dX(true);     
	
	CaLV(V);  //计算最小二乘残差
	double VTPV=CaVTPV();

	if(m_IsPx)
	{
		double pvv=Ca_VTPV_XYZ()+VTPV;
		Sigma=sqrt(pvv/(3*Vnumber-3*unPnumber));
	}
	else
	{
		Sigma=sqrt(VTPV/(3*Vnumber-3*unPnumber));
	}	
	PrintXYZ();
	PrintLV("V",V);

	fprintf(resultfp,"\n\n最小二乘平差:\n[pvv]=%e",VTPV);
	fprintf(resultfp,"\n单位权中误差μ=±%e",Sigma);	
			 
}

//////////////////////////////////////////////////////////////////////
//  最小二乘平差函数
void CCoGPS::Helmert()
{

	for(int kk=1;kk<=2;kk++)
	{
		CaATPA();     	 
		Known();
		double max=Ca_dX(1);     
		
		CaLV(V);  //计算最小二乘残差
		PrintLV("V",V);
		
		double VTPV=0.0;
		
		for(int i=0;i<Anumber;i++) // Anumber：向量组总数
		{
			double *Pi=P[i];     //i区权矩阵起始地址
			int k1=3*dir0[i];    //i区首向量在观测值中的序号
			int k2=3*dir0[i+1];  //i+1区首向量在观测值中的序号
			
			double *PVi=PV+k1;
			double *Vi=V+k1;
			int n=k2-k1;
			
			double pvv=0.0;
			int d0 = 3*dir0[i];
			for(int j=0;j<n;j++)
			{
				double pvj=0.0;
				for(int k=0; k<n; k++) pvj+=Pi[ij(j,k)]*Vi[k];
				
				PVi[j]=pvj;
				pvv += pvj*Vi[j];
			}
			
			
			fprintf(resultfp,"\n[%d]=%lf  \n",i+1,pvv/n);

			double m2=pvv/n;
//			if(m2<0.1)m2=0.1;
			
			for(j=0; j<(n+1)*n/2; j++)Pi[j]/=m2;
			
			//if(kk==2)inverse(Pi,n);
			VTPV+=pvv;
		}

			
	}
//	PrintData();
//	return;

	double VTPV=CaVTPV();
	if(m_IsPx)
	{
		double pvv=Ca_VTPV_XYZ()+VTPV;
		Sigma=sqrt(pvv/(3*Vnumber-3*unPnumber));
		
	}
	else
	{
		Sigma=sqrt(VTPV/(3*Vnumber-3*unPnumber));
	}	
	PrintXYZ();
	PrintLV("V",V);
	
	fprintf(resultfp,"\n     验后精度估计:\n");
	fprintf(resultfp,"\n     [pvv]:    %e",VTPV);
	fprintf(resultfp,"\n     m:        %e(验后单位权中误差)",Sigma);
}




//////////////////////////////////////////////////////////////////////
// 计算残差V的权倒数
void CCoGPS::caQV()
{
	//计算观测值平差值权倒数
	for(int i=0;i<3*Vnumber; i++)
    {
        int k1=3*dir1[i/3]+i%3;  //测站点未知数编号
        int k2=3*dir2[i/3]+i%3;  //照准点未知数编号
		
		double q=ATPA[ij(k1,k1)]+ATPA[ij(k2,k2)]-2.0*ATPA[ij(k1,k2)];

		QV[i] = QV[i]-q;
    }	
}



//////////////////////////////////////////////////////////////////////
//  计算中位数、:   (1)计算单位权残差; (2)中位数计算方差因子
void CCoGPS::Ca_UintV()
{
	for(int i=0; i<3*Vnumber; i++)
	{
		double mi=sqrt(QV[i]);
		if(fabs(QV[i])<1.0e-15) W[i] = 1.0;
		else W[i]=V[i]/mi;  
	}

    MEDIAN=Median(W,3*Vnumber,true)*1.4826;
	
	//fprintf(resultfp,"\nMEDIAN=%10e\n",MEDIAN);
}






//////////////////////////////////////////////////////////////////////////
//     输出坐标平差值及其中误差
void CCoGPS::PrintXYZ()
{
	double m=Sigma;
	double Maxmx=0.0;
	double Maxmy=0.0;
	double Maxmz=0.0;

	fprintf(resultfp,"\n\n\n STATION       X              Y              Z          RMS_X  RMS_Y  RMS_Z");
	fprintf(resultfp,"\n              (m)            (m)            (m)          (cm)   (cm)   (cm)\n");
	fprintf(resultfp," ----------------------------------------------------------------------------\n");

	for(int k=0;k<Pnumber;k++)
	{
		fprintf(resultfp,"\n%4d %4s",k+1,Pname[k]);

		double X=XYZ[3*k];
		double Y=XYZ[3*k+1];
		double Z=XYZ[3*k+2];				 

		fprintf(resultfp," %14.4lf %14.4lf %14.4lf",X,Y,Z);

		double mx,my,mz;

		mx=sqrt(ATPA[ij(3*k,3*k)])*m*100.0;
		my=sqrt(ATPA[ij(3*k+1,3*k+1)])*m*100.0;
		mz=sqrt(ATPA[ij(3*k+2,3*k+2)])*m*100.0;

		fprintf(resultfp," %6.2lf %6.2lf %6.2lf",mx,my,mz);

		if(mx>Maxmx)Maxmx=mx;
		if(my>Maxmy)Maxmy=my;
		if(mz>Maxmz)Maxmz=mz;

	}
	fprintf(resultfp,"\n\n  Max(RMS_X)=%5.2lf",Maxmx);
	fprintf(resultfp,"\n  Max(RMS_Y)=%5.2lf",Maxmy);
	fprintf(resultfp,"\n  Max(RMS_Z)=%5.2lf",Maxmz);
}




//////////////////////////////////////////////////////////////////////////
//    已知点加（相关）权平差时，坐标改正数的二次型
double CCoGPS::Ca_VTPV_XYZ()
{
	 int kn=Pnumber-unPnumber; //已知点数
	 double *Vxyz=new double [3*kn];

	 for(int i=0;i<kn;i++)
	 {
		 int k=KnownPoints[i];  //已知点的点号
		 double X=XYZ[3*k];
		 double Y=XYZ[3*k+1];
		 double Z=XYZ[3*k+2];				 

		 Vxyz[3*i]  =X-KnownXYZ[3*i];				 
		 Vxyz[3*i+1]=Y-KnownXYZ[3*i+1];				 
		 Vxyz[3*i+2]=Z-KnownXYZ[3*i+2];
	 }

	 double VTPV_XYZ=0.0;
	 for(i=0;i<3*kn;i++)
	 {
		 for(int j=0;j<3*kn;j++)
			 VTPV_XYZ+=Vxyz[i]*KnownPx[ij(i,j)]*Vxyz[j];
	 }

	 delete []Vxyz;
	 return VTPV_XYZ;
}




//////////////////////////////////////////////////////////////////////////
//	 点名存贮，返回点名对应的点号
int  CCoGPS::GetStationNumber(char *buf)
{
    for(int i=0; i<Pnumber; i++) 
	{
		if(Pname[i]!=NULL)
		{
			//与已经存入点名地址数组的点名比较
			if(strcmp(buf,Pname[i])==0)return i;
		}
		else
		{
			//将新点名存到内存，地址放到Pname数组中
			int len = strlen(buf);
			Pname[i] = new char[len+1];
			strcpy(Pname[i], buf);
			return i;
		}
	}
	return -1;
}


//////////////////////////////////////////////////////////////////////////
//   GPS网抗差估计
void CCoGPS::Robust(int fname,double eps)
{
	int n=3*Vnumber;    //观测值总数
	QV=new double [n];  //残差的权倒数
	W =new double [n];  //权因子
	
	//------------------------------------
	//  计算观测值的权矩阵
	int ii=0;  // 计数器
	for(int i=0; i<Anumber; i++)
	{
		int ni = 3*(dir0[i+1]-dir0[i]); //向量组观测值总数
		for(int j=0;j<ni;j++)
		{
			QV[ii]=P[i][ij(j,j)];
			ii++;
		}
		inverse(P[i],ni);
	}
	//--------------------------------------
	//  计算残差的权倒数
	CaATPA();   //组成法方程
	if(!inverse(ATPA,3*Pnumber))  //法方程系数矩阵求逆
	{
		MyBreak("法方程系数降秩降秩！");
		exit(0);
	}
	caQV();  //计算残差的权倒数
	
	//--------------------------------------
	//  将观测值的权备份到数组PP	
	double **PP= new double *[Anumber]; //等价权数组
	for(int k=0; k<Anumber; k++)
	{
		int nk=3*(dir0[k+1]-dir0[k]);
		PP[k]=new double[nk*(nk+1)/2];
		for(int i=0;i<nk*(nk+1)/2;i++) PP[k][i]=P[k][i];
	}
	
	//---------------------------------------
	// 抗差估计迭代
	for( i=0;i<n;i++) W[i]=1.0;//权因子赋初值
	int No=0; //迭代次数
	while(1)
	{
		No++;
		// 计算等价权
		for(int k=0; k<Anumber; k++) // k－向量组编号
		{
			int nk = 3*(dir0[k+1]-dir0[k]); //观测值个数
			double *Wk=W+3*dir0[k];         //权因子的数组地址
			
			int index = 0;
			for(int i=0;i<nk;i++)
			{
				for(int j=0;j<=i;j++)
				{
					if(Wk[i]<1.0 || Wk[j]<1.0)
						P[k][index] = PP[k][index]*sqrt(Wk[i]*Wk[j]);
					else P[k][index] = PP[k][index];
					index++;
				}
			}
		}

		CaATPA(); //组成法方程
		Known();  // 已知点处理
		double maxX=Ca_dX(0);     
		CaLV(V);  //计算残差

		//fprintf(resultfp,"\n抗差估计，第%d次迭代...\nmax(dX)=%e",No,maxX);
		// 		PrintLV("抗差估计残差",V);
		
		if(fabs(maxX)<eps && No>1)
		{
			break;
		}

		//-----------------------------------------------
		//权因子计算
		Ca_UintV();  //标准化残差计算(放在W里)
		for(int i=0;i<3*Vnumber;i++)
		{
			W[i]=Wi(fname,W[i]/MEDIAN,5.0 , 8.0);
		}
	}

	int dd=0; //降权观测值的个数
	for(i=0;i<3*Vnumber;i++)
	{
		if(W[i]<1.0)dd++;
	}
	
	if(dd==0)
	{
		fprintf(resultfp,"\n  所有观测值均未降权!\n");	
	}
	else
	{   //输出降权因子
		fprintf(resultfp,"\n  共有%d个观测值被降权!\n",dd);	
		fprintf(resultfp,"\n\n权因子: \n");
		for(int i=0;i<3*Vnumber;i++)
		{
			if(W[i]<1.0)
			{
				int ni=AreaNumber[i/3];
				fprintf(resultfp,"\n%8d",ni+1);
				fprintf(resultfp,"%10s %10s",
					Pname[dir1[i/3]],Pname[dir2[i/3]]);
				fprintf(resultfp,"   d%c %10.4lf",'x'+i%3,W[i]);
			}
		}  
	}
				
	double VTPV=CaVTPV();
	Sigma=sqrt(VTPV/(n-dd-3*unPnumber));
	PrintXYZ();
	PrintLV("V",V);

	fprintf(resultfp,"\n\n抗差估计:\n[pvv]=%e",VTPV);
	fprintf(resultfp,"\n单位权中误差μ=±%e",Sigma);	
	
	delete []PP;
	//PrintM(resultfp,W,3*Vnumber,3,"%6.3lf ","迭代权因子");	
}

//////////////////////////////////////////////////////////////////////////
//     拟稳平差
void CCoGPS::Quasi_Stable(char *file)
{
	//将观测值的方差协方差矩阵求逆，化为权矩阵
	for(int i=0; i<Anumber; i++)
	{
		int ni = 3*(dir0[i+1]-dir0[i]); //向量组观测值总数
		inverse(P[i],ni);
	}
	
	CaATPA();
	
	FILE *fp=fopen(file,"r");
	if(fp==NULL) 
	{
		MyBreak("拟稳点信息文件打不开！");
		exit(0);
	}
	
	bool *IsStable = new bool[Pnumber]; //是否为拟稳点
	for( i=0;i<Pnumber;i++)
	{
		IsStable[i]=false;
	}
	int StableNumber; //拟稳点总数
	fscanf(fp,"%d",	&StableNumber);
	for(i=0;i<StableNumber;i++)
	{
		char buf[50];
		fscanf(fp,"%s",buf);
		int k=GetStationNumber(buf);
		if(k<0)
		{
			MyBreak("拟稳点名文件中点名错误！"); 
			exit(0);
		}
		IsStable[k]=true;
	}
	fclose(fp);

	double f=1.0/StableNumber;
	for(i=0; i<3*Pnumber; i++)
	{
		if(!IsStable[i/3]) continue;//表示拟稳点，跳过去
		for(int j=0; j<=i;j++)
		{
			if(i%3==j%3 && IsStable[j/3])
				ATPA[ij(i,j)]+=f;
		}
	}
	
	delete []IsStable;
	Ca_dX(true);
	
	for( i=0; i<3*Pnumber; i++)
	{
		for(int j=0; j<=i;j++)
		{
			if(i%3==j%3) ATPA[ij(i,j)]-=f;
		}
	}
	
	CaLV(V);
	double VTPV=CaVTPV();
	Sigma = sqrt(VTPV/(3*Vnumber-3*(Pnumber-1)));
	
	PrintXYZ();
	fprintf(resultfp,"\n\n拟稳平差:\n[pvv]=%e",VTPV);
	fprintf(resultfp,"\n单位权中误差μ=±%e",Sigma);	
}

//////////////////////////////////////////////////////////////////////////
//  对称正定矩阵求逆(仅存下三角元素)
bool inverse2(double a[],int n)
{
    double *Y=new double[n];
    for(int k=0;k<n;k++)
	{ 
		double a0=a[0];
		if(a0+1.0==1.0)
		{
			delete []Y;  return false;
		}
		int m=n-k-1;
		for(int i=1;i<n;i++)
		{
			double g=a[i*(i+1)/2];
			Y[i]=g/a0;
			if(i<=m)Y[i]=-Y[i];
			for(int j=1;j<=i;j++)
			{
				a[(i-1)*i/2+j-1]=a[i*(i+1)/2+j]+g*Y[j];
			} 
		}
		a[n*(n+1)/2-1]=1.0/a0;
		for(i=1;i<n;i++)
			a[(n-1)*n/2+i-1]=Y[i];
	} 
	delete []Y;
	return true;
}


///////////////////////////////////////////////////////////////
// 计算粗差估值的权逆阵Qg,粗差估值gross,返回粗差二次型VTPVg
double CCoGPS::CaQg(
			int    s,      //粗差参数的个数
			int    Jint[], //加了粗差参数的观测值序号（集）
			double gross[],//粗差向量
			double Qg[],   //粗差权逆阵
			double JTPA[],  //PA的子矩阵
			double JTPV[] ) //V的子向量

{
     CaJTPJ( Qg, Jint, s);
     CaJTPA(JTPA,Jint[s-1],s);	 

	 // 粗差估值的权逆阵计算
     int t=3*Pnumber;   //参数总数
	 for(int j=0;j<s; j++) 
	 { 
		 for(int i=0; i<=j; i++)
		 {
			 for(int k=0; k<t; k++)
			 {
				 double a=JTPA[j*t+k];
				 if(a+1.0==1.0)continue;
				 for(int l=0;l<t;l++)
				 {
					 Qg[ij(j,i)]-=a*ATPA[ij(k,l)]*JTPA[i*t+l];
				 }
			 }
		 }
	 }  
	 if(!inverse2(Qg,s))return -2.0;
	 
	 // 计算粗差估值
	 JTPV[s-1]=PV[Jint[s-1]];
	 for(int i=0; i<=s-1; i++)
     {
		 gross[i]=0.0;
		 for(int j=0;j<=s-1;j++)
			 gross[i]-=Qg[ij(i,j)]*JTPV[j];
	 }

     // 计算粗差二次型
	 double VTPVg=0.0;
     for(i=0; i<=s-1; i++)
	 {
		 VTPVg-=JTPV[i]*gross[i];
	 }
     return VTPVg;
}

#include "Probability.h"

//////////////////////////////////////////////////////////////////////////
//    粗差探测函数      2008-08-29
void CCoGPS::DataSnooping(double Arfa,int Max)  // Max-最大搜索次数
{
     //  为数组申请内存
	 double *JTPA=new double [3*Pnumber*Max];
     double *JTPV=new double [Max];
     double *Qg= new double [Max*(Max+1)/2]; // 粗差权逆阵
     double *gross= new double [Max];  // 粗差估值
     int    *Jint= new int [Max];  // 可疑观测值序号

     CProbability MyProb; //概率计算类
	 
	 //-------------------------------------------------
	 //  最小二乘平差
	 for(int i=0; i<Anumber; i++)
	 {
		 int ni = 3*(dir0[i+1]-dir0[i]); //向量组观测值总数
		 inverse(P[i],ni);
	 }
	 
	 CaATPA();    
	 Known();	 
	 Ca_dX(false);     
	 CaLV(V);  //计算最小二乘残差
	 double VTPV=CaVTPV();

	 //---------------------------------------------------
	 //粗差探测
     fprintf(resultfp,"\n           逐个搜索法粗差探测:\n");

     int maxp;  //当前搜索中，统计量最大的观测值序号
     int Gn=0;  //当前搜索次数
	 double  new_Omg;      // 本次搜索最大粗差二次型;
     double  last_Omg=0.0; // 前次搜索最大粗差二次型
  	 int r=3*Vnumber-3*unPnumber; // 多余观测数		 
	 
	 bool *W=new bool[3*Vnumber];  //粗差标志数组
	 for( i=0; i<3*Vnumber; i++)W[i]=0;
			 
     while(1)
     {
        if(Gn+1>Max)
		{
			MyBreak("给定的搜索次数过小，请增大搜索次数，重新开始搜索。");
			fclose(resultfp);
			exit(0);
		}

  	    fprintf(resultfp,"\n第%d次搜索:",Gn+1);
		
        double Flimit=MyProb.re_F(1,r-(Gn+1),1.0-Arfa);

        new_Omg=0.0; 
        for(int i=0; i<3*Vnumber; i++)
		{
			 Jint[Gn]=i;

             if(W[i]) continue;	// 检查号观测值是否为已知粗差观测值

			 double Omg = CaQg(Gn+1,Jint,gross,Qg,JTPA,JTPV); //粗差二次型
             if(Omg<0.0)continue;

             if(Omg>new_Omg && Omg<VTPV)
			 {
				 new_Omg=Omg; 
				 maxp=i;  
			 }
		}

        double F=(new_Omg-last_Omg)/(VTPV-new_Omg)*(r-(Gn+1));
		   
		if(F>Flimit)fprintf(resultfp,"\n   F=%6.3lf(>%6.3lf)",F,Flimit);
		else	    fprintf(resultfp,"\n   F=%6.3lf(<%6.3lf)",F,Flimit);

		Jint[Gn]=maxp;

		double Omg=CaQg(Gn+1,Jint,gross,Qg,JTPA,JTPV);
 
		for(i=0;i<Gn+1;i++)
		{
			int j=Jint[i];
			fprintf(resultfp,"\n   No.%-2d",j);
			fprintf(resultfp,"  向量组:%-3d ",AreaNumber[j/3]+1); 
			fprintf(resultfp,"%5s ",Pname[dir1[j/3]]);
			fprintf(resultfp,"%5s ",Pname[dir2[j/3]]);
			fprintf(resultfp," (%c)",'X'+j%3);
			fprintf(resultfp,"%12.4lf",gross[i]);
		}

        if( F<=Flimit || fabs(gross[Gn])<0.01)break;
		Gn++;
        last_Omg=new_Omg;
		W[maxp]=1;
    }

    if(Gn==0)
    {
	      fprintf(resultfp,"\n\n粗差探测结果：未发现粗差");
    }
	else
	{
		 for(int j=1;j<=Gn;j++)
		 {
			 new_Omg=CaQg(j,Jint,gross,Qg,JTPA,JTPV);
		 }
		 fprintf(resultfp,"\n\n 逐个搜索法最后定位:(粗差改正)\n");
		 fprintf(resultfp," 粗差总数: %d\n",Gn);
		 
		 for(int i=0;i<Gn;i++)
		 {
			 j=Jint[i];
			fprintf(resultfp,"\n   No.%-2d",j);
			fprintf(resultfp,"  向量组:%-3d ",AreaNumber[j/3]+1); 
			fprintf(resultfp,"%5s ",Pname[dir1[j/3]]);
			fprintf(resultfp,"%5s ",Pname[dir2[j/3]]);
			fprintf(resultfp," (%c)",'X'+j%3);
			fprintf(resultfp,"%12.4lf",gross[i]);
		 }
		 
		 fprintf(resultfp,"\n\n剔除粗差后最小二乘平差结果");
		 VTPV=VTPV-new_Omg;
		 Sigma=sqrt(VTPV/(r-Gn));
		 Ca_X_Qx(gross,Qg,JTPA,Gn);
		 
		 PrintXYZ();
		 CaLV(V);//计算残差
		 PrintLV("V",V);	  
		 fprintf(resultfp,"\n\n粗差探测:\n[pvv]=%e",VTPV);
		 fprintf(resultfp,"\n单位权中误差μ=±%e",Sigma);	
		 
	}


	delete []JTPA;
	delete []JTPV;
	delete []Qg;
	delete []gross;
	delete []Jint;
	delete []W;
}


//////////////////////////////////////////////////////////////////////
//        计算JTPA的最后一行(n0-1)
//         JTPA是n0行t列的矩阵
void CCoGPS::CaJTPA(double JTPA[],int ks,int s) 
//s-可疑观测值总数 ;  ks-最后一个可疑观测值的编号
{
     int t=3*Pnumber;
     for(int i=0;i<t;i++)JTPA[(s-1)*t+i]=0.0;

	 int k =AreaNumber[ks/3];  // 第ks号观测值的组号
     int k1=dir0[k];           // k组首向量的序号
     int k2=dir0[k+1];         // k+1组区首向量的序号
	 int nk=3*(k2-k1);

	 int s2=ks-3*k1; // 组内观测值编号
	 int *Dir1=dir1+k1;
	 int *Dir2=dir2+k1;
	 
	 for(int s1=0;s1<nk;s1++)
     {
         int i1=3*Dir1[s1/3]+s1%3; 
         int i2=3*Dir2[s1/3]+s1%3; 
		 
		 double p12 = P[k][ij(s1,s2)];
		 		
		 JTPA[(s-1)*t+i1]-=p12;
         JTPA[(s-1)*t+i2]+=p12;
     }

	 /*   
     int num=AreaNumber[ks/3];  // 第ks号观测值的组号
     int k1=dir0[num];         //  num组首向量的序号
     int k2=dir0[num+1];       //  num+1区首向量的序号


     int j=ks-3*k1;
     for(int k=3*k1;k<3*k2;k++)
     {
         int i=k-3*k1;
         int ij=(j<=i) ? i*(i+1)/2+j : j*(j+1)/2+i;
         double pij = P[num][ij];

         int i1=3*dir1[k/3]+k%3;
         int i2=3*dir2[k/3]+k%3;
         JTPA[(s-1)*t+i1]-=pij;
         JTPA[(s-1)*t+i2]+=pij;
     }

*/

}


//////////////////////////////////////////////////////////////////////
//	计算JTPJ，s-粗差参数的个数
void CCoGPS::CaJTPJ(double JTPJ[], int Jint[],int s)
{
	for(int k=0; k<s; k++)
	{
		int k1=Jint[k]; //观测值编号
		int A1=AreaNumber[k1/3]; //观测值所属同步区编号
		double *Pk = P[A1];      //本同步区权矩阵数组
		int k0 = 3*dir0[A1];     //本同步区首观测值的编号    			

		for(int i=0;i<=k;i++)
		{ 
			int k2 = Jint[i]; //观测值编号
			if(AreaNumber[k2/3] == A1) //k1,k2属于同一同步区
			{
				JTPJ[ij(k,i)]=Pk[ij(k1-k0,k2-k0)];
			}
			else JTPJ[ij(k,i)]=0.0;
		}	  
	}
}


//////////////////////////////////////////////////////////////////////////
//   计算剔除s个粗差后坐标的平差值及权逆阵
void CCoGPS::Ca_X_Qx(double gross[],double Qg[],double JTPA[],int s)
// gross[] -粗差向量
// Qg[] - 粗差估值权逆阵
// J'PA 
// s      粗差个数
{
	int t=3*Pnumber;

	// M = inverse(N) * A'PJ
	double *M=new double[t*s];
	for(int i=0;i<t*s;i++) M[i]=0.0;
	for(i=0;i<t;i++)
	{
		for(int j=0;j<s;j++)
		{
			for(int k=0; k<t; k++)
				M[i*s+j]+=ATPA[ij(i,k)]*JTPA[j*t+k];
		}
	}

	for(i=0;i<t;i++)
	{
		double xi=0.0;
		for(int j=0;j<s;j++)
			xi+=M[i*s+j]*gross[j];
		XYZ[i]-=xi;
	}

	for(i=0;i<t; i++)
	{
		for(int j=0;j<=i; j++)
		{
			double nij=0.0;
			for(int k=0; k<s; k++)
				for(int h=0; h<s; h++)
					nij+=M[i*s+k]*Qg[ij(k,h)]*M[j*s+h];
			ATPA[ij(i,j)]+=nij;
		}
	}
	delete []M;
	return;
}




//////////////////////////////////////////////////////////////////////////
//   输入已知点的坐标、方差－协方差矩阵(已知点相关)、中误差(已知点独立)
bool CCoGPS::InputKnownXYZ(char *file,bool correlative)
{
	m_IsPx = true; //参数加权平差

    FILE *fp=fopen(file,"r");  
	if(fp==NULL )
	{
		MyBreak("已知坐标文件打不开");
		return 0;
	}
	
	int kn;
	fscanf(fp,"%d",&kn);   //已知点总数
	unPnumber=Pnumber-kn;
	
	int tt=(3*kn)*(3*kn+1)/2;   //权矩阵元素的个数
	KnownPoints=new int [kn];   //已知点点号
	KnownXYZ=new double [3*kn]; //已知点坐标
	KnownPx=new double [tt];    //已知点坐标的权矩阵
	for(int i=0;i<tt;i++) KnownPx[i]=0.0;
	for(i=0;i<Pnumber;i++)IsKnown[i]=false; //初始值：未知点

	if(correlative)//已知坐标是相关的
	{
		for(i=0;i<kn;i++)
		{
			char name[20];
			fscanf(fp,"%s%lf%lf%lf",name,
				KnownXYZ+3*i,KnownXYZ+3*i+1,KnownXYZ+3*i+2); 
			int k=GetStationNumber(name);
			if(k<0)
			{
				MyBreak("读取已知坐标文件出错：点名错误！");  
				return 0;
			}
			
			KnownPoints[i]=k; //点号
			IsKnown[k]=true;
		} 
		
		int jj=0;
		for(i=0;i<3*kn;i++)
		{
			int d;
			fscanf(fp,"%d",&d);   // 读取行号
			for(int j=0;j<=i;j++)
			{ 
				fscanf(fp,"%lf",KnownPx+jj); //方差协方差
				jj++;
			}
		} 
		
		fprintf(resultfp,"\n 已知点坐标:\n");
		for(i=0;i<kn;i++)
		{
			int k=KnownPoints[i];
			fprintf(resultfp,"\n%10s",Pname[k]);
			fprintf(resultfp," %15.5lf %14.5lf %14.5lf",
				KnownXYZ[3*i],KnownXYZ[3*i+1],KnownXYZ[3*i+2]);
		}
		PrintM2(resultfp,KnownPx,3*kn,4,"%18.10e ","已知点坐标的方差阵");

		bool rn=inverse(KnownPx,3*kn); //方差协方差矩阵求逆，的权矩阵
		if(!rn) 
		{
			MyBreak("输入已知坐标文件出错：坐标方差阵降秩");
			return 0;
		}
	}
	else //已知坐标是独立的
	{
		for(i=0;i<kn;i++)
		{
			char name[20];
			fscanf(fp,"%s%lf%lf%lf",name,
				KnownXYZ+3*i,KnownXYZ+3*i+1,KnownXYZ+3*i+2); 
			
			int k=GetStationNumber(name);
			if(k<0)
			{
				MyBreak("读取已知坐标文件出错：点名错误！");    
				return 0;
			}
			
			KnownPoints[i]=k;
			IsKnown[k]=true;
						
			double mx,my,mz;
			fscanf(fp,"%lf%lf%lf",&mx,&my,&mz);
			KnownPx[ij(3*i,3*i)] = 1.0/(mx*mx);
			KnownPx[ij(3*i+1,3*i+1)] = 1.0/(my*my);
			KnownPx[ij(3*i+2,3*i+2)] = 1.0/(mz*mz);
		} 
		fprintf(resultfp,"\n 已知点坐标及其中误差:\n");
		for(i=0;i<kn;i++)
		{
			int k=KnownPoints[i];
			fprintf(resultfp,"\n%10s",Pname[k]);
			fprintf(resultfp," %15.5lf %14.5lf %14.5lf",
				KnownXYZ[3*i],KnownXYZ[3*i+1],KnownXYZ[3*i+2]);

			double mx=KnownPx[ij(3*i,3*i)];
			double my=KnownPx[ij(3*i+1,3*i+1)];
			double mz=KnownPx[ij(3*i+2,3*i+2)];
			fprintf(resultfp," %6.4lf %6.4lf %6.4lf",
							sqrt(1.0/mx),sqrt(1.0/my),sqrt(1.0/mz));
		}
	}
	fclose(fp);	
	unPnumber=Pnumber-kn;
	return 1;
}

