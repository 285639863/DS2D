#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include <math.h>

#include "alloc.h"

float dotP (float **x1, float **x2, int n1, int n2)
{

        float result;
        result = 0.0f;

        for (int i2=0; i2<n2; ++i2)
            for (int i1=0; i1<n1; ++i1)
                result += x1[i2][i1]*x2[i2][i1];

        return result;
}


float *mul1(float *x1, float *x2, int n)
{
      float *pt;
      pt = (float*)malloc(sizeof(float)*n);
      for (int i=0; i<n; ++i)
          pt[i]=x1[i]*x2[i];
      
      return pt ;
}

float **mul2(float **x1, float **x2, int n1, int n2)
{
      float **pt;
      pt = alloc2float(n1, n2);
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i2][i1]=x1[i2][i1]*x2[i2][i1];

      return pt ;
}


float ***mul3(float ***x1, float ***x2, int n1, int n2, int n3)
{
      float ***pt;
      pt = alloc3float(n1, n2, n3);
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i3][i2][i1]=x1[i3][i2][i1]*x2[i3][i2][i1];

      return pt ;
}



float *mul1c(float *x1, float c, int n)
{
      float *pt;
      pt = (float*)malloc(sizeof(float)*n);
      for (int i=0; i<n; ++i)
          pt[i]=x1[i]*c;
      
      return pt ;
}

float **mul2c(float **x1, float c, int n1, int n2)
{
      float **pt;
      pt = alloc2float(n1, n2);
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i2][i1]=x1[i2][i1]*c;

      return pt ;
}


float ***mul3c(float ***x1, float c, int n1, int n2, int n3)
{
      float ***pt;
      pt = alloc3float(n1, n2, n3);
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i3][i2][i1]=x1[i3][i2][i1]*c;

      return pt ;
}



float *add1(float *x1, float *x2, int n)
{
      float *pt;
      pt = (float*)malloc(sizeof(float)*n);
      for (int i=0; i<n; ++i)
          pt[i]=x1[i]+x2[i];
      
      return pt ;
}

float **add2(float **x1, float **x2, int n1, int n2)
{
      float **pt;
      pt = alloc2float(n1, n2);
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i2][i1]=x1[i2][i1]+x2[i2][i1];

      return pt ;
}


float ***add3(float ***x1, float ***x2, int n1, int n2, int n3)
{
      float ***pt;
      pt = alloc3float(n1, n2, n3);
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i3][i2][i1]=x1[i3][i2][i1]+x2[i3][i2][i1];

      return pt ;
}

float *sub1(float *x1, float *x2, int n)
{
      float *pt;
      pt = (float*)malloc(sizeof(float)*n);
      for (int i=0; i<n; ++i)
          pt[i]=x1[i]-x2[i];
      
      return pt ;
}

float **sub2(float **x1, float **x2, int n1, int n2)
{
      float **pt;
      pt = alloc2float(n1, n2);
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i2][i1]=x1[i2][i1]-x2[i2][i1];

      return pt ;
}


float ***sub3(float ***x1, float ***x2, int n1, int n2, int n3)
{
      float ***pt;
      pt = alloc3float(n1, n2, n3);
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i3][i2][i1]=x1[i3][i2][i1]-x2[i3][i2][i1];

      return pt ;
}



float *div1c(float *x1, float c, int n)
{
      float *pt;
      pt = (float*)malloc(sizeof(float)*n);
      for (int i=0; i<n; ++i)
          pt[i]=c/x1[i];
      
      return pt ;
}

float **div2c(float **x1, float c, int n1, int n2)
{
      float **pt;
      pt = alloc2float(n1, n2);
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i2][i1]=c/x1[i2][i1];

      return pt ;
}


float ***div3c(float ***x1, float c, int n1, int n2, int n3)
{
      float ***pt;
      pt = alloc3float(n1, n2, n3);
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i3][i2][i1]=c/x1[i3][i2][i1];

      return pt ;
}


float max1(float *x, int n)
{
      float m ;
      m = -9999999999;
      for(int i=0; i<n; ++i){
         if(x[i]>m)m=x[i];
      }
      
      return m;
}



float max2(float **x, int n1, int n2)
{
      float m ;
      m = -9999999999;
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1){
               if(x[i2][i1]>m)m=x[i2][i1];
      }
      
      return m;
}


float max3(float ***x, int n1, int n2, int n3)
{
      float m ;
      m = -9999999999;
      for (int i3=0; i3<n3; ++i3)
         for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1){
               if(x[i3][i2][i1]>m)m=x[i3][i2][i1];
      }
      
      return m;
}



float min1(float *x, int n)
{
      float m ;
      m = 9999999999;
      for(int i=0; i<n; ++i){
         if(x[i]<m)m=x[i];
      }
      
      return m;
}



float min2(float **x, int n1, int n2)
{
      float m ;
      m = 9999999999;
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1){
               if(x[i2][i1]<m)m=x[i2][i1];
      }
      
      return m;
}


float min3(float ***x, int n1, int n2, int n3)
{
      float m ;
      m = 9999999999;
      for (int i3=0; i3<n3; ++i3)
         for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1){
               if(x[i3][i2][i1]<m)m=x[i3][i2][i1];
      }
      
      return m;
}



float *sqrt1(float *x1,  int n)
{
      float *pt;
      pt = (float*)malloc(sizeof(float)*n);
      for (int i=0; i<n; ++i)
          pt[i]=sqrt(x1[i]);
      
      return pt ;
}

float **sqrt2(float **x1,  int n1, int n2)
{
      float **pt;
      pt = alloc2float(n1, n2);
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i2][i1]=sqrt(x1[i2][i1]);

      return pt ;
}


float ***sqrt3(float ***x1,  int n1, int n2, int n3)
{
      float ***pt;
      pt = alloc3float(n1, n2, n3);
      for (int i3=0; i3<n3; ++i3)
        for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i3][i2][i1]=sqrt(x1[i3][i2][i1]);

      return pt ;
}


float **abs2(float **x1,  int n1, int n2)
{
      float **pt;
      pt = alloc2float(n1, n2);
      for (int i2=0; i2<n2; ++i2)
           for (int i1=0; i1<n1; ++i1)
               pt[i2][i1]=abs(x1[i2][i1]);

      return pt ;
}
