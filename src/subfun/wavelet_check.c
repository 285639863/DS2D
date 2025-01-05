#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include <math.h>

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif

#include "alloc.c"
#include "rw.c"

void    gauss0 (int nt,  float f, float dt,float *s);
void    gauss1 (int nt,  float f, float dt,float *s);
void    gauss2 (int nt,  float f, float dt,float *s);
void    ricker1 (int nt,  float f, float dt,float *s);



int main (int argc, char *argv[])
{
    
    int     nt = 2000;
    int     f0 = 25;
    float   dt = 0.0005;
    
    float   *sou = NULL;

    sou = alloc1float(nt);  

    gauss2  (nt, f0, dt, sou);
    ricker1 (nt, f0, dt, sou);

    write_1d_float_wb(sou, nt, "wavelet_ricker.dat");    

    exit (0);

}

void    gauss0 (int nt,  float f, float dt,float *s)
{
        float pi = 3.14159265358979f;
        for(int i=0; i<nt; i++){
           float tt = i*dt-100*dt;
           float sp = pi*f*tt;
           s[i] = -2000.*pi*pi*f*f*tt*exp(-sp*sp)*(3.-2.*sp*sp);
        }
}

void    gauss1 (int nt,  float f, float dt,float *s)
{
        float pi = 3.14159265358979f;
        float t0 = 80/1000.0;
        int   kt = (int)(t0/dt); 
        for(int i=0; i<nt; i++){
           float tt = i*dt-kt*dt;
           float sp = pi*f*tt;
           s[i] = -2000.*pi*pi*f*f*tt*exp(-sp*sp)*(3.-2.*sp*sp);
        }
}

void    gauss2 (int nt,  float f, float dt,float *s)
{
        float pi = 3.14159265358979f;
        float t0 = 1/f;
        int   kt = (int)(t0/dt); 
        for(int i=0; i<nt; i++){
           float tt = i*dt-kt*dt;
           float sp = pi*f*tt;
           s[i] = -2000.*pi*pi*f*f*tt*exp(-sp*sp)*(3.-2.*sp*sp);
        }
}

void    ricker1 (int nt,  float f, float dt,float *s)
{
        float pi = 3.14159265358979f;
        float timeshift = 1/f;
        float pi2 = sqrt(pi)/2.0f;
        float b = sqrt(6.0f)/(pi*f);
        float constant = 2.0f*sqrt(6.0f)/b;
        float smax = 0.0f;
        float t1, t2, amp, u,amax;

        amax = 0 ;
        for (int i=0; i<nt; ++i) {
          t1 = dt*(float)i; 
          t2 = t1 - timeshift;
          u = constant*t2;
          s[i] = -1.0f*((u*u)/4.0f-0.5f)*pi2*(exp(-u*u/4.0f));
          amax = MAX(fabs(s[i]),amax);
        }

        for (int i=0; i<nt; ++i) s[i] /= amax ;
        //s = div(s,max(abs(s)));

}

