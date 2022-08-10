// Chapter_7.cpp: implementation of the CChapter_7 class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "math.h"
#include "Probability.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//      伽玛函数的自然对数值
double CProbability::lngam(double x)
{
      double a,f,z,gam;
      int i;
      if(x<7.0)
      {
	     a=7.0-floor(x);
	     z=x-1;
	     f=1.0;
	     for(i=1;i<=a;i++)
		 { 
	       z=z+1.0;
	       f=f*z;
		 } 
	     x=z+1;
	     f=-log(f);
      }
      else f=0.0;
      z=1.0/x/x;
     
	  gam=f+(x-0.5)*log(x)-x+0.91893853320+
	       (((-0.00059523809524*z+0.00079365079365)
	           *z-0.0027777777778)*z+0.083333333333)/x;
      return gam;
}


//////////////////////////////////////////////////////////////////////////
//	  正态分布函数值:  p(-∞,u)
double CProbability::norm(double u)
{
      double y,er,p,q;
      y=fabs(u)/sqrt(2.0);
      if(u<-5.0)return 0.0;
      if(u>5.0)return 1.0;

      p=1.0+y*(0.0705230784+y*(0.0422820123+
	    y*(0.0092705272+y*(0.0001520143+
	    y*(0.0002765672+y*0.0000430638)))));
      p=-16.0*log(p);
      er=1.0-exp(p);
      q=0.5*er;
      p=(u<0.0)? 0.5-q: 0.5+q;
      return p;

}



//////////////////////////////////////////////////////////////////////////
// 正态分布的上侧分位数, p(-∞,u)=p ; 已知p, 返回u
double CProbability::re_norm(double p)
{
      double u,y;

      if(p==0.5) return 0.0;
	  int sign=1;
      if(p<0.5) return -re_norm(1.0-p);

//      if(p>0.999999999999||(1.0-p>0.999999999999))
//                       MyBreak("\n Error in \"pnorm.c\".");
      y=-log(4.0*p*(1.0-p));
      y=y*(1.570796+y*(0.3706988e-1
	      +y*(-0.836435359e-3+y*(-0.2250947e-3
	      +y*(0.68412183e-5+y*(0.58242385e-5
	      +y*(-0.1045275e-5+y*(0.8360937e-7
	      +y*(-0.3231081e-8+y*(0.3657763e-10
	      +y*0.6936234e-12))))))))));
      u=sqrt(y);

      return u;

}

/*******************************************
      chi 方分布函数值:区间(0,x)上的概率p
      已知x,n, 求p, d
      d: x处的概率密度值
      n: 自由度
 返回值: ch点的分布函数值q
*******************************************/
double CProbability::chi2(int n,double x,double *d)
{
     int i,iai;
     double lu,p,pi,pis,y,chs,u,pp,n2;

     pi=3.14159265358979312;

     pis=sqrt(pi);
     y=x/2.0;
     chs=sqrt(x);
     if(n%2)
     {
		 lu=log(sqrt(y))-y-log(pis);
	     pp=norm(chs);
	     p=2.0*(pp-0.5);
	     iai=1;
     }
     else
     {
		 lu=log(y)-y;
	     p=1.0-exp(-y);
	     iai=2;
     }
     u=exp(lu);
     if(iai!=n)
     {
		 n2=n-2;
         for(i=iai;i<=n2;i=i+2)
		 { 
			 p=p-2*u/i;
	         lu=lu+log(x/i);
	         u=exp(lu);
		 } 
     } 
     *d=u/x;
     return p;

}



//////////////////////////////////////////////////////////////////////////
//	chi方分布的上侧分位数:  p=p(0,x)
//  已知概率值p,自由度n，反求x 
//  返回值: 上侧分位数x
double CProbability::re_chi2(int n,double p)
{
     if((p>0.9999)) return -1.0e30;
     if(n==1)
     {
		 double x=re_norm((1.0-p)/2.0);
		 return x*x;
     }
	 
     if(n==2)return -2.0*log(1.0-p);

     double u=re_norm(p);
     double w=2.0/(9.0*n);
     double x0=1.0-w+u*sqrt(w);
     x0=n*x0*x0*x0;
     
	 int k=0;
     while(1)
     {
		 double f;
		 double pp=chi2(n,x0,&f);
		 if(f+1.0==1.0)return x0;
		 double xx=x0-(pp-p)/f;
		 if(fabs(x0-xx)<0.0001*fabs(xx)) return xx;
		 k++;
		 if(k>=30)return xx;
		 x0=xx;
     }

}


//////////////////////////////////////////////////////////////////////////
//      F分布函数值:区间(0,f)上的概率q
//      已知f,n1,n2, 求q,d
//      d: f处的概率密度值
double CProbability::Fdist(int n1,int n2, double f,double *d)
{
      int iai,ibi;
      double p,pi,x,u;

      pi=3.14159265358979312;
      x=n1*f/(n2+n1*f);

      if(n1%2 && n2%2)   /*  n1,n2均是奇数  */
      {
	      u=sqrt(x*(1.0-x))/pi;
	      p=1.0-2.0*atan(sqrt((1.0-x)/x))/pi;
	      iai=ibi=1;
      }

      if(n1%2 && !(n2%2) ) /*  n1是奇数,n2是偶数  */
      {
	       u=sqrt(x)*(1.0-x)/2.0;
	       p=sqrt(x);
	       iai=1;
	       ibi=2;
      }

      if(!(n1%2) && n2%2)  /*  n1是偶数,n2是奇数  */
      {
	        u=x*sqrt(1.0-x)/2.0;
	        p=1.0-sqrt(1.0-x);
	        iai=2;
	        ibi=1;
      }

      if(!(n1%2) && !(n2%2))  /*  n1,n2均是偶数  */
      {
	        u=x*(1.0-x);
	        p=x;
	        iai=ibi=2;
      }

      while((ibi!=n2) || (iai!=n1))
      {
	     if(ibi!=n2)
		  {
	         p=p+2.0*u/ibi;
	         u=u*(1.0+iai/(ibi+0.0))*(1.0-x);
	         ibi=ibi+2;
		  }
	     if(iai!=n1)
		  {
	       p=p-2.0*u/iai;
	       u=u*(1.0+ibi/(0.0+iai))*x;
	      iai=iai+2;
		  }
      }
      *d=u/f;
      return p;
}

//////////////////////////////////////////////////////////////////////////
//    F分布的反函数：p=F(0,x), 已知p,反求x
double CProbability::re_F(int n1,int n2,double p)
{
     double df12,df22,a,a1,b,b1,yq,e,f0,ff,pp,d,ga1,ga2,ga3,log2;

     df12=n1/2.;
     df22=n2/2.;
     a=2.0/(9.0*n1);
     a1=1.0-a;
     b=2.0/(9.0*n2);
     b1=1.0-b;
     yq=re_norm(p);
     e=b1*b1-b*yq*yq;

     log2=log(2.0);
     if(e<=0.8)
     {
	     ga1=lngam(df12+df22);
	     ga2=lngam(df12);
	     ga3=lngam(df22);
	     f0=(2.0/n2)*(log2+(df22-1.0)*log(n2+0.0)+ga1-ga2-ga3
				      -df22*log(n1+0.0)-log(1-p));
	     f0=exp(f0);
     }
     else
     {
	     f0=(a1*b1+yq*sqrt(a1*a1*b+a*e))/e;
	     f0=f0*f0*f0;
     }
     if(n1>=40 && n2>=40) return f0;

     while(1)
     {
	   double df;
	   pp=Fdist(n1,n2,f0,&d);
	   if(d+1.0==1.0)return f0;
	   df=(p-pp)/d;
	   while(fabs(f0)<2.0*fabs(df))df/=10.0;
	   ff=f0+df;
	   if(fabs(f0-ff)<0.00001)return ff;
	   f0=ff;
     }

}

//////////////////////////////////////////////////////////////////////////
//   t分布的反函数：p=F(x),已知p，反求x
double CProbability::re_t(int n,double p)
{
	double a,b;

	if(p>0.5)
	{ 
		a=0.0;
		b=5.0;
	}
	else
	{
		a=-5.0;
		b=0.0;
	}


	while(1)
	{
		double pa=t(n,a);
     	double pb=t(n,b);
		if(pa<=p && p<=pb)
		{
			if(b-a<0.0001)return a;
			double ab=(a+b)/2.0;
			double pm=t(n,ab);

			if(p>pm)a=ab;
			else b=ab;
		}
		else if(p<pa)a=a-5.0;
		else if(p>pb)b=b+5.0;

	}
	return 0.0;
}

//////////////////////////////////////////////////////////////////////////
//     t分布的分布密度值，计算累计函数值时，要调用这个函数
double CProbability::ff(int n, double x)
{
	static int nn=-1;
	static double gm;
	double pai=3.14159265358979312;
	if(n!=nn)
	{
		nn=n;
		gm=exp(lngam(0.5*(n+1.0))-lngam(0.5*(n)));
		gm=gm/sqrt(n*pai);

   	}

	double z=pow(1.0+x*x/n,-0.5*(n+1));

	return gm*z;

}

//////////////////////////////////////////////////////////////////////////
//        t分布的分布函数值（负无穷到x的积分值）  //nn 自由度
double CProbability::t(int nn, double x) 
{
	if(x+1.0==1.0)return 0.5;

	static int nnn=-1; //上次的nn
	static double xx;  //上次的x
	static double ss;  //上次的结果

	double a,b;
	if(nn==nnn)
	{
		if(x>xx){ a=xx; b=x; }
		else    { b=xx; a=x; }
	}
	else 
	{
		a=0.0;
		b=fabs(x);
	}
	
	double s2;

	int n=1;
	double h=b-a;
	double t1=h*(ff(nn,a)+ff(nn,b))/2.0;
	double s1=t1;

	double eps=0.00000001;
	double ep=eps+1.0;
	while(ep>=eps)
	{
		double p=0.0;
		for(int k=0;k<=n-1;k++)
		{
			double x=a+(k+0.5)*h;
			p=p+ff(nn,x);
		}
		double t2=(t1+h*p)/2.0;
		s2=(4.0*t2-t1)/3.0;
		ep=fabs(s2-s1);
		t1=t2;
		s1=s2;
		n=+n+n;
		h=h/2.0;
	}

	if(nn==nnn)
	{
		if(x>xx)s2=ss+s2;
		else    s2=ss-s2; 
	}
	else 
	{
		if(x>0.0)s2=s2+0.5;
		else s2=0.5-s2;
	}

	ss=s2;
	nnn=nn;
	xx=x;

	return s2;
}



