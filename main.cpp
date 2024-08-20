//#define _CRT_SECURE_NO_WARNINGS
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include <math.h>
#include <sys/stat.h>
#define PI 3.14159265358979
#define FILE_NAME_MAX_LENGTH 256
    
#ifndef MAX 
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#include "laplace.h"
#include "./subfun/alloc.c"

//#include "../subfun/rw.c"
#include "./subfun/segy.h"
//#include "./subfun/read_write_su.c"
//#include "matrixc.c"

//OpenMP
#include<omp.h>
#include "mpi.h"

#include "FD2DGPU.cuh"

#include <algorithm>

// extern float *alloc1float(size_t n1);
// extern int **alloc2int(size_t n1, size_t n2);
// extern void free1float(float *p);
// extern void free2int(int **p);
// extern void free2 (void **p);
// extern void **alloc2 (size_t n1, size_t n2, size_t size);


void normalize(float* matrix, int length, int tracenum) {
	float max = 0;
	for (int i = 0; i < length * tracenum; i++)
		max = (fabs(matrix[i]) > max) ? fabs(matrix[i]) : max;
	//normalized	
	for (int i = 0; i < length * tracenum; i++) {
		matrix[i] = matrix[i] / max;
	}
}


int    index_shot(char *fn, int *nt, float *dt, int *ns, int **table,int coordinate_scale=1)
{
       bhed   Bh;
       segy1  Th;
       Y_3200 Y3200;
       FILE   *fp ;
       int    cx_min_s, cx_max_s;
       int    ntr, pos;

       fp = fopen(fn,"rb");
       if(fp == NULL) {
          printf("Sorry,cann't open seismic file!\n");
          return 1;
       }
 

       fseek(fp, 0, SEEK_SET);
       fread(&Th, 240, 1, fp); //240 can be displaced by sizeof(segy1)
       *nt = (int)Th.ns;
       *dt = Th.dt/1000000.0;

	// std::cout<<*nt<<std::endl;
	// std::cout<<*dt<<std::endl;

       int TL = (*nt)*sizeof(float) ;
       
       fseek(fp, 0, SEEK_SET);
     
       *ns     = 0;
       pos     = 0;
       int sx0 = -999999;
       int sy0 = -999999;
       for( ; ; ){
          fread(&Th, 240, 1, fp);
          if(feof(fp)){
             int ns0 = *ns ;
             table[ns0-1][0] = ns0;
             table[ns0-1][1] = ntr;
             table[ns0-1][2] = sx0/coordinate_scale;
             table[ns0-1][3] = 0; //sy
             table[ns0-1][4] = cx_min_s/coordinate_scale;
             table[ns0-1][5] = 0;
             table[ns0-1][6] = cx_max_s/coordinate_scale;
             table[ns0-1][7] = 0;
             table[ns0-1][8] = pos - ntr;
             break;
          }

          int sx = Th.sx;
          int sy = 0;//2d,sy=0
          int gx = Th.gx;
          int gy = 0;
     
	// if(sx==14795640||sx==14790640){ 
	// 	std::cout<<"sx = "<<sx<<"gx = "<<gx<<std::endl;
	// 	sleep(0.1);
	// }    

          int xmin = MIN(sx, gx);
          int xmax = MAX(sx, gx);

          if(sx != sx0){
              if(pos > 0){
                  int ns0 = *ns ;
                  table[ns0-1][0] = ns0;
                  table[ns0-1][1] = ntr;
                  table[ns0-1][2] = sx0/coordinate_scale;
                  table[ns0-1][3] = 0; //sy
                  table[ns0-1][4] = cx_min_s/coordinate_scale;
                  table[ns0-1][5] = 0;
                  table[ns0-1][6] = cx_max_s/coordinate_scale;
                  table[ns0-1][7] = 0;
                  table[ns0-1][8] = pos - ntr;
              }
              (*ns) ++;
              if((*ns)%50==0)printf(" %dth shot has been indexed!\n", (*ns));
        //       printf("ns:%d ntr:%d\n",*ns,ntr);
		// sleep(3);
		// if(*ns==5)return 0;
              ntr = 1;

              sx0 = sx;
              sy0 = sy;
             
              cx_min_s = 99999999 ;
              cx_max_s = -999999 ;
          }else{
              ntr ++ ; 
//		std::cout<<"ntr = "<<ntr<<std::endl;
//		sleep(1); 
          }

          pos ++ ;
          if(xmin < cx_min_s) cx_min_s = xmin;
          if(xmax > cx_max_s) cx_max_s = xmax;      
  
          fseek(fp, TL, SEEK_CUR);

	
      } 
      fclose(fp); 
   
      return 0 ;
}

void index_shot(int ns,int **table,int dis_shot,int disx)
{
	for(int i=0;i<ns;i++)
	{
		table[i][0] = i+1;
		table[i][1] = (disx/2)+dis_shot*i;
		table[i][2] = dis_shot*i;
		table[i][3] = disx+dis_shot*i;
	}

}


void index_shot_marchenko(int ns,int **table,int dis_shot,int disx)
{
	for(int i=0;i<ns;i++)
	{
		table[i][0] = i+1;
		table[i][1] = dis_shot*i;
		table[i][2] = 0;
		table[i][3] = disx;
	}

}


void index_shot_update(int min_shot,int max_shot,int **table,float dis_shot,float disx)
{
	for(int i=0;i<=max_shot-min_shot;i++)
	{
		table[i][0] = i+1;
		table[i][1] = dis_shot*(i+min_shot);
		table[i][2] = 0;
		table[i][3] = disx;
	}

}





void index_shot_disx(int ns,int **table,int dis_shot,int disx,int distancex)
{
	for(int i=0;i<ns;i++)
	{
		table[i][0] = i+1;
		table[i][1] = dis_shot*i;
		if(table[i][1]<=(disx/2))
		{
			table[i][2] = 0;
			table[i][3] = disx;
		}
		else if(table[i][1]>(disx/2)&&table[i][1]<distancex-(disx/2))
		{
			table[i][2] = table[i][1] - (disx/2);
			table[i][3] = table[i][1] + (disx/2);
		}
		else        //table[i][1]>=distancex-(disx/2)
		{ 
			table[i][2] = distancex - disx;
			table[i][3] = distancex; 
		}
	}
}


void   read_shot_gc_su(char *fn, long long int pos, int ntr, int nt,int *gc)
{
	int i;
    FILE   *fp ;
     float *dat = new float[ntr*nt];

    fp = fopen(fn,"rb");
    if(fp == NULL) {
            printf("Sorry,cann't open input seismic file!\n");
            exit(0);
    }
          
    int TL = 240 + nt*sizeof(float);
           
    fseek(fp, (long long int)TL*pos, SEEK_SET);

	segy1 Th;

    for(i=0; i<ntr; i++){
           fread(&Th, 240, 1, fp);
           gc[i]=Th.gx;               
           fread(&dat[i*nt], sizeof(float), nt, fp);
    }
    delete[] dat;
    fclose(fp);
}



void   read_shot_gather_su(char *fn, long long int pos, int ntr, int nt, float *dat,int *gc,int coordinate_scale=1)
{
	int i;
    FILE   *fp ;
     
    fp = fopen(fn,"rb");
    if(fp == NULL) {
            printf("Sorry,cann't open input seismic file!\n");
            exit(0);
    }
          
    int TL = 240 + nt*sizeof(float);
           
    fseek(fp, (long long int)TL*pos, SEEK_SET);

	segy1 Th;

    for(i=0; i<ntr; i++){
           fread(&Th, 240, 1, fp);
           gc[i]=Th.gx/coordinate_scale;               
           fread(&dat[i*nt], sizeof(float), nt, fp);
    }

    fclose(fp);
}

void   read_shot_gather_su2(char *fn, long long int pos, int ntr, int nt, float *dat,int *gc,int coordinate_scale=1)
{
	int i;
    FILE   *fp ;
    
	std::cout<<"==================enter differ read shot gather su======================="<<std::endl;
 
    fp = fopen(fn,"rb");
    if(fp == NULL) {
            printf("Sorry,cann't open input seismic file!\n");
            exit(0);
    }
          
    int TL = 240 + nt*sizeof(float);
           
    fseek(fp, (long long int)TL*pos, SEEK_SET);

	segy1 Th;

    for(i=0; i<ntr; i++){
           fread(&Th, 240, 1, fp);
           gc[i]=Th.gx/coordinate_scale;         
//	   std::cout<<Th.gx<<std::endl;
//        	sleep(0.5); 
           fread(&dat[i*nt], sizeof(float), nt, fp);
    }

    fclose(fp);
}



void    setPml(float **v, int nx, int nz, float dz, float dx, int pml,float **damp_pml)
{
        int   nzpml = nz+2*pml;
        int   nxpml = nx+2*pml;
        float *damp ;
        float vmin  = 999999 ;

        for(int i=0; i<nxpml; i++)
          for(int j=0; j<nzpml; j++){
            vmin = MIN(v[i][j], vmin);
            //if(vmin < 500)printf(" WARNING! Vmin =%f,i=%d,j=%d\n", vmin,i,j);
            }
        if(vmin < 500)printf(" WARNING! Vmin =%f\n", vmin);
        vmin=5000;

        damp = alloc1float(pml);

        float a = dx*(float)(pml-1);
        //float b = 1.5f*vmin*log(1000.0f)/a;
        float b = 1.5f*vmin*log(10000000.0f)/a;

        float xa;

        for (int ix=0; ix<pml; ++ix) {
          xa = (float)ix/(pml-1.0f);
          damp[ix] = b*xa*xa;
        }

        for (int ix=0; ix<pml; ++ix) {
          for (int iz=0; iz<nzpml; ++iz) {
            damp_pml[pml-ix-1    ][iz] = damp[ix];
            damp_pml[nxpml-pml+ix][iz] = damp[ix];
          }
        }

        a = dz*(float)(pml-1);
        //b = 1.5f*vmin*log(1000.0f)/a;
        b = 1.5f*vmin*log(10000000.0f)/a;

        float za;

        for (int iz=0; iz<pml; ++iz) {
          za = (float)iz/(pml-1.0f);
          damp[iz] = b*za*za;
        }

        for (int iz=0; iz<pml; ++iz) {
          for (int ix=pml-iz; ix<nxpml-(pml-iz); ++ix) {
            damp_pml[ix][pml-iz-1    ] = damp[iz];
            damp_pml[ix][nzpml-pml+iz] = damp[iz];
          }
        }

}

void    setPml1d(float *v, int nx, int nz, float dz, float dx, int pml,float *damp_pml)
{
        int   nzpml = nz+2*pml;
        int   nxpml = nx+2*pml;
        float *damp ;
        float vmin  = 999999 ;

        for(int i=0; i<nxpml; i++)
          for(int j=0; j<nzpml; j++){
            vmin = MIN(v[i*nzpml+j], vmin);
            //if(vmin < 500)printf(" WARNING! Vmin =%f,i=%d,j=%d\n", vmin,i,j);
            }
        if(vmin < 500)printf(" WARNING! Vmin =%f\n", vmin);
        // vmin=5000;

        damp = alloc1float(pml);

        float a = dx*(float)(pml-1);
        //float b = 1.5f*vmin*log(1000.0f)/a;
        float b = 1.5f*vmin*log(10000000.0f)/a;

        float xa;

        for (int ix=0; ix<pml; ++ix) {
          xa = (float)ix/(pml-1.0f);
          damp[ix] = b*xa*xa;
        }

        for (int ix=0; ix<pml; ++ix) {
          for (int iz=0; iz<nzpml; ++iz) {
            damp_pml[(pml-ix-1)*nzpml+iz] = damp[ix];
            damp_pml[(nxpml-pml+ix)*nzpml+iz] = damp[ix];
          }
        }

        a = dz*(float)(pml-1);
        //b = 1.5f*vmin*log(1000.0f)/a;
        b = 1.5f*vmin*log(10000000.0f)/a;

        float za;

        for (int iz=0; iz<pml; ++iz) {
          za = (float)iz/(pml-1.0f);
          damp[iz] = b*za*za;
        }

        for (int iz=0; iz<pml; ++iz) {
          for (int ix=pml-iz; ix<nxpml-(pml-iz); ++ix) {
            damp_pml[ix*nzpml+pml-iz-1] = damp[iz];
            damp_pml[ix*nzpml+nzpml-pml+iz] = damp[iz];
          }
        }

}

void    setVel(float **v,int nx, int nz, int pml)
{
        int nzpml = nz+2*pml;
        int nxpml = nx+2*pml;

        for (int iz=0; iz<pml; ++iz)
            for (int ix=0; ix<nxpml; ++ix) {
                v[ix][iz] = v[ix][pml];
                v[ix][nzpml-pml+iz] = v[ix][nzpml-pml-1];
            }

        for (int ix=0; ix<pml; ++ix)
            for (int iz=0; iz<nzpml; ++iz) {
                v[ix][iz] = v[pml][iz];
                v[nxpml-pml+ix][iz] = v[nxpml-pml-1][iz];
            }

}

void    setVel1d(float *v,int nx, int nz, int pml)
{
        int nzpml = nz+2*pml;
        int nxpml = nx+2*pml;

        for (int iz=0; iz<pml; ++iz)
            for (int ix=0; ix<nxpml; ++ix) {
                v[ix*nzpml+iz] = v[ix*nzpml+pml];
                v[ix*nzpml+nzpml-pml+iz] = v[ix*nzpml+nzpml-pml-1];
            }

        for (int ix=0; ix<pml; ++ix)
            for (int iz=0; iz<nzpml; ++iz) {
                v[ix*nzpml+iz] = v[pml*nzpml+iz];
                v[(nxpml-pml+ix)*nzpml+iz] = v[(nxpml-pml-1)*nzpml+iz];
            }

}



void    ricker1 (int nt,  float f, float dt,float *s)
{
        float pi = 3.14159265358979f;
        float t0 = 1/f;
        int   kt = (int)(t0/dt);

        for(int i=0; i<nt; i++){
           float tt = i*dt-kt*dt;
           float sp = pi*f*tt;
           //s[i] = 1000.*exp(-sp*sp)*(1.-2.*sp*sp);
           s[i] = 1.0*exp(-sp*sp)*(1.-2.*sp*sp); 
        }
}

void  negative_staggered_grid_extrapolation_2d_pml(float **p,float **px,float **pz,float **vx,float **vz,float **damp,float **rho,float **k,float dt,int nzpml,int nxpml,float dz,float dx,int sord)
{
      float  C[6];
      float  C1, C2;
      int    nop ;


      if( sord == 4 ){
          C[0] = 0;
          C[1] = 9.0/8.0;
          C[2] = -1.0/24.0;
          C[3] = 0;
          C[4] = 0;
          C[5] = 0;
      }

      else if( sord == 6 ){
	  C[0] = 0;
          C[1] = 75.0/64.0;
          C[2] = -25.0/384.0;
          C[3] = 3.0/640.0;
          C[4] = 0;
          C[5] = 0;
      }

      else if( sord == 8 ){
          C[0] = 0;
          C[1] = 1225.0/1024.0;
          C[2] = -245.0/3072.0;
          C[3] = 49.0/5120.0;
          C[4] = -5.0/7168.0;
          C[5] = 0;
      }

      else if( sord == 10 ){
          C[0] = 0;
          C[1] = 711.0/587.0;
          C[2] = -158.0/1761.0;
          C[3] = 271.0/19577.0;
          C[4] = -25.0/14159.0;
          C[5] = 18.0/151669.0;
      }


      else{
          printf(" space approximation order not in 4,6,8,10 please check!\n ");
//          goto end;
      }

      nop = sord/2;

	int ix,iz;
	float dampx1,dampx2,dampz1,dampz2;
	float tmp_p,tmp_vx,tmp_vz;
	
#pragma omp parallel for private(ix,iz,dampx1,dampx2,dampz1,dampz2,tmp_vx,tmp_vz)
      for(ix=nop-1; ix<nxpml-nop; ++ix)
          for(iz=nop-1; iz<nzpml-nop; ++iz){
	      dampx1 = 1 - dt*damp[ix][iz] / 2;
	      dampx2 = 1 + dt*damp[ix][iz] / 2;
	      dampz1 = 1 - dt*damp[ix][iz] / 2;
	      dampz2 = 1 + dt*damp[ix][iz] / 2;

	      tmp_vx = 0;
	      tmp_vz = 0;

              for (int i=1; i<=nop;i++){
                    tmp_vx += C[i]*(vx[ix+i][iz]-vx[ix-i+1][iz]);
                    tmp_vz += C[i]*(vz[ix][iz+i]-vz[ix][iz-i+1]);
              }

	      px[ix][iz] = (dampx1*px[ix][iz]+k[ix][iz]*(dt/dx)*tmp_vx)/dampx2;
	      pz[ix][iz] = (dampz1*pz[ix][iz]+k[ix][iz]*(dt/dz)*tmp_vz)/dampz2;

//	      px[ix][iz] = (dampx1*px[ix][iz]-k[ix][iz]*(dt/dx)*tmp_vx)/dampx2;
//	      pz[ix][iz] = (dampz1*pz[ix][iz]-k[ix][iz]*(dt/dz)*tmp_vz)/dampz2;
              p[ix][iz]  = px[ix][iz]+pz[ix][iz];

      }


#pragma omp parallel for private(ix,iz,dampx1,dampx2,tmp_p)
      for (ix=nop; ix<nxpml-nop+1; ++ix)
          for (iz=nop; iz<nzpml-nop; ++iz) {

	      dampx1 = 1 - dt*damp[ix][iz] / 2;
	      dampx2 = 1 + dt*damp[ix][iz] / 2;


	      tmp_p = 0;	   

              for (int i=1; i<=nop;i++){
                    tmp_p += C[i]*(p[ix+i-1][iz]-p[ix-i][iz]);
              }

	      vx[ix][iz] = (dampx1*vx[ix][iz]+(1.0/rho[ix][iz])*(dt/dx)*tmp_p)/dampx2;

//	      vx[ix][iz] = (dampx1*vx[ix][iz]-(1.0/rho[ix][iz])*(dt/dx)*tmp_p)/dampx2;

      }




#pragma omp parallel for private(ix,iz,dampz1,dampz2,tmp_p)
      for (ix=nop; ix<nxpml-nop; ++ix)
          for (iz=nop; iz<nzpml-nop+1; ++iz) {

	      dampz1 = 1 - dt*damp[ix][iz] / 2;
	      dampz2 = 1 + dt*damp[ix][iz] / 2;


	      tmp_p = 0;	   

              for (int i=1; i<=nop;i++){
                    tmp_p += C[i]*(p[ix][iz+i-1]-p[ix][iz-i]);
              }

	      vz[ix][iz] = (dampz1*vz[ix][iz]+(1.0/rho[ix][iz])*(dt/dz)*tmp_p)/dampz2;

//	      vz[ix][iz] = (dampz1*vz[ix][iz]-(1.0/rho[ix][iz])*(dt/dz)*tmp_p)/dampz2;
      }

   end: printf("");

}

void  staggered_grid_extrapolation_2d_pml(float **p,float **px,float **pz,float **vx,float **vz,float **damp,float **rho,float **k,float dt,int nzpml,int nxpml,float dz,float dx,int sord)
{
      float  C[6];
      float  C1, C2;
      int    nop ;


      if( sord == 4 ){
          C[0] = 0;
          C[1] = 9.0/8.0;
          C[2] = -1.0/24.0;
          C[3] = 0;
          C[4] = 0;
          C[5] = 0;
      }

      else if( sord == 6 ){
	  C[0] = 0;
          C[1] = 75.0/64.0;
          C[2] = -25.0/384.0;
          C[3] = 3.0/640.0;
          C[4] = 0;
          C[5] = 0;
      }

      else if( sord == 8 ){
          C[0] = 0;
          C[1] = 1225.0/1024.0;
          C[2] = -245.0/3072.0;
          C[3] = 49.0/5120.0;
          C[4] = -5.0/7168.0;
          C[5] = 0;
      }

      else if( sord == 10 ){
          C[0] = 0;
          C[1] = 711.0/587.0;
          C[2] = -158.0/1761.0;
          C[3] = 271.0/19577.0;
          C[4] = -25.0/14159.0;
          C[5] = 18.0/151669.0;
      }


      else{
          printf(" space approximation order not in 4,6,8,10 please check!\n ");
          goto end;
      }

      nop = sord/2;

	int ix,iz;
	float dampx1,dampx2,dampz1,dampz2;
	float tmp_p,tmp_vx,tmp_vz;
	
#pragma omp parallel for private(ix,iz,dampx1,dampx2,dampz1,dampz2,tmp_vx,tmp_vz)
      for(ix=nop-1; ix<nxpml-nop; ++ix)
          for(iz=nop-1; iz<nzpml-nop; ++iz){
	      dampx1 = 1 - dt*damp[ix][iz] / 2;
	      dampx2 = 1 + dt*damp[ix][iz] / 2;
	      dampz1 = 1 - dt*damp[ix][iz] / 2;
	      dampz2 = 1 + dt*damp[ix][iz] / 2;

	      tmp_vx = 0;
	      tmp_vz = 0;

              for (int i=1; i<=nop;i++){
                    tmp_vx += C[i]*(vx[ix+i][iz]-vx[ix-i+1][iz]);
                    tmp_vz += C[i]*(vz[ix][iz+i]-vz[ix][iz-i+1]);
              }

//	      px[ix][iz] = (dampx1*px[ix][iz]+k[ix][iz]*(dt/dx)*tmp_vx)/dampx2;
//	      pz[ix][iz] = (dampz1*pz[ix][iz]+k[ix][iz]*(dt/dz)*tmp_vz)/dampz2;

	      px[ix][iz] = (dampx1*px[ix][iz]-k[ix][iz]*(dt/dx)*tmp_vx)/dampx2;
	      pz[ix][iz] = (dampz1*pz[ix][iz]-k[ix][iz]*(dt/dz)*tmp_vz)/dampz2;
              p[ix][iz]  = px[ix][iz]+pz[ix][iz];

      }


#pragma omp parallel for private(ix,iz,dampx1,dampx2,tmp_p)
      for (ix=nop; ix<nxpml-nop+1; ++ix)
          for (iz=nop; iz<nzpml-nop; ++iz) {

	      dampx1 = 1 - dt*damp[ix][iz] / 2;
	      dampx2 = 1 + dt*damp[ix][iz] / 2;


	      tmp_p = 0;	   

              for (int i=1; i<=nop;i++){
                    tmp_p += C[i]*(p[ix+i-1][iz]-p[ix-i][iz]);
              }

//	      vx[ix][iz] = (dampx1*vx[ix][iz]+(1.0/rho[ix][iz])*(dt/dx)*tmp_p)/dampx2;

	      vx[ix][iz] = (dampx1*vx[ix][iz]-(1.0/rho[ix][iz])*(dt/dx)*tmp_p)/dampx2;

      }




#pragma omp parallel for private(ix,iz,dampz1,dampz2,tmp_p)
      for (ix=nop; ix<nxpml-nop; ++ix)
          for (iz=nop; iz<nzpml-nop+1; ++iz) {

	      dampz1 = 1 - dt*damp[ix][iz] / 2;
	      dampz2 = 1 + dt*damp[ix][iz] / 2;


	      tmp_p = 0;	   

              for (int i=1; i<=nop;i++){
                    tmp_p += C[i]*(p[ix][iz+i-1]-p[ix][iz-i]);
              }

//	      vz[ix][iz] = (dampz1*vz[ix][iz]+(1.0/rho[ix][iz])*(dt/dz)*tmp_p)/dampz2;

	      vz[ix][iz] = (dampz1*vz[ix][iz]-(1.0/rho[ix][iz])*(dt/dz)*tmp_p)/dampz2;
      }

   end: printf("");

}


void ricker_wave(float *w, int Tn, float dt, float FM)
{
	int t;
	float  Nk = PI*PI*FM*FM*dt*dt;
	int t0 = ceil(1.0 / (FM*dt));          
	for (t = 0; t<Tn; t++)
	{
		w[t] = (1.0 - 2.0*Nk*(t - t0)*(t - t0))*exp(-Nk*(t - t0)*(t - t0));
	}
}

// 		lsrtm   use   model   par
typedef struct{
    char fn1[1024];
    char fn2[1024];	
    char imagedir[1024];	    
    float *velp;
    float *vels;
    float *sou;
    float *record_z;
    float *record_x;
    int *gc;
    float dx,dz,dt;
    int minshot,maxshot,nx,nz,ns,nxpml,nzpml,allnx,allnz,scale,pml,nt,nop,ntr_pre;
    bool light,rbc,cpu_mem;
    int iointerval;
    int all_left;
} modelpar;


void init_modelparameters(modelpar *model,float dx,float dz,float dt,int minshot,int maxshot, int nx, int nz,\
            int ns,int nxpml,int nzpml,int allnx,int allnz,int scale,int pml,int nt,int nop,int ntr_pre,bool light,bool rbc,bool cpu_mem,int iointerval,int all_left){
    model->dx = dx; 
    model->dz = dz;
    model->dt = dt;
    model->minshot = minshot;
    model->maxshot = maxshot;
    model->nx = nx;
    model->nz = nz;    
    model->ns = ns;
    model->nxpml = nxpml;
    model->nzpml = nzpml;
    model->allnx = allnx;
    model->allnz = allnz;
    model->scale = scale;
    model->pml = pml;
    model->nt = nt;
    model->nop = nop;
    model->ntr_pre = ntr_pre;
    model->light  = light;
    model->rbc = rbc;
    model->cpu_mem = cpu_mem;
    model->iointerval = iointerval;
    model->all_left = all_left;
}
//

// elsrtm one iteration ///
void lsrtm_all_iteration(int myid,int np,int sy,int gy,MPI_Status status,int **table,modelpar *model,float *image_pp,float *image_ps,float *pp_grad,float *ps_grad,float *illumination,int maxiter,int record_left_in_v){

//	MPI_Status status;
//    	int myid,np;
	// MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	// MPI_Comm_size(MPI_COMM_WORLD,&np);


    float *b_pp_grad = new float[model->allnx*model->allnz];
    float *b_ps_grad = new float[model->allnx*model->allnz];
    float *b_pp_grad2 = new float[model->allnx*model->allnz];
    float *b_ps_grad2 = new float[model->allnx*model->allnz];

    float *sum_pp_grad = new float[model->allnx*model->allnz];
    float *sum_ps_grad = new float[model->allnx*model->allnz];
    if(myid!=0){
        cudaSetDevice((myid-1)%8);			
        // FD2DGPU_ELASTIC image2delastic(model->sou,model->dx,model->dz,model->dt,model->nxpml,model->nzpml,model->allnx,model->allnz,model->scale,model->pml,model->nt,model->nop,model->ntr_pre);	
        // image2delastic.GPUbufferVPVS(model->velp,model->vels);
    }

    int *members = new int[np-1];
    for(int i=0;i<np-1;i++){
        members[i]=i+1;
    }

    MPI_Group group_world,group_new;
    MPI_Comm groupcomm;
    MPI_Comm_group(MPI_COMM_WORLD,&group_world);
    MPI_Group_incl(group_world, np-1, members, &group_new);
    MPI_Comm_create(MPI_COMM_WORLD, group_new, &groupcomm);


    // int id_ingroup,np_ingroup;
	// MPI_Comm_rank(groupcomm,&id_ingroup);
	// MPI_Comm_size(groupcomm,&np_ingroup);

    // if(myid!=0){
    //     std::cout<<"myid = "<<myid<<std::endl;    //<<"\t id_ingroup = "<<id_ingroup<<"\t np_ingroup = "<<np_ingroup<<std::endl;
    // }

	int ip;
	int send[9],recv[9];
	int nsend,ntask;
	ntask = model->ns;   

    double misfit =0;
    double max_misfit =0;
    double rms_misfit =0;

    double fenzi =0;
    double fenmu =0;        

    double beta_fenzi =0;
    double beta_fenmu =0;
    double beta =0;

    float *pp_cg = new float[model->allnx*model->allnz]{};
    float *ps_cg = new float[model->allnx*model->allnz]{};

    float *pp_cg_old = new float[model->allnx*model->allnz]{};
    float *ps_cg_old = new float[model->allnx*model->allnz]{};

    float array_misfit[2][maxiter];

    float *s_k_p = new float[model->allnx*model->allnz];
    float *s_k_s = new float[model->allnx*model->allnz];
    float *r_k_p = new float[model->allnx*model->allnz];
    float *r_k_s = new float[model->allnx*model->allnz];
    float *s_k0_p = new float[model->allnx*model->allnz];
    float *s_k0_s = new float[model->allnx*model->allnz];
    float *r_k0_p = new float[model->allnx*model->allnz];
    float *r_k0_s = new float[model->allnx*model->allnz];    
    float *s_kmod_p = new float[model->allnx*model->allnz];
    float *s_kmod_s = new float[model->allnx*model->allnz];
    float *r_kmod_p = new float[model->allnx*model->allnz];
    float *r_kmod_s = new float[model->allnx*model->allnz]; 

    float *delta_mp = new float[model->allnx*model->allnz];
    float *delta_ms = new float[model->allnx*model->allnz];

    float beta1 = 0.8;
    float beta2 = 0.999;
    // float beta2 = 0.9;    
    float alpha = 3e-5;

    // float epsilon = 3e7;     
    float epsilon = 5e8;       
    float gama = 0.99;

    // float beta1 = 0.7;
    // float beta2 = 0.999;
    // float alpha = 5e-5;
    // float epsilon = 1e7;
    // float gama = 0.99;

    // float beta1 = 0.7;
    // float beta2 = 0.999;
    // float alpha = 1e-5;
    // float epsilon = 3e7;
    // float gama = 0.99;

    int k;
    float alpha_k;

    
	if(myid==0)
	{
    for(int iter=0;iter<maxiter;iter++)
    {
        printf("======================================================\n\n");      
        printf(" iter           : %d\n",iter);
        printf(" maxiter        : %d\n",maxiter);
        printf(" process        : %f%\n\n",(float)iter/maxiter*100);
        printf("======================================================\n\n");      
		nsend = 0;
		for(int i=0;i<ntask+np-1;i++)
		{
			MPI_Recv(recv,9,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			ip = status.MPI_SOURCE;
			if(i<model->ns)
			{
                         send[0] = table[i][0];
                         send[1] = table[i][1];
                         send[2] = table[i][2];
                         send[3] = table[i][3];
                         send[4] = table[i][4];
                         send[5] = table[i][5];
                         send[6] = table[i][6];
                         send[7] = table[i][7];
                         send[8] = table[i][8];
			}
			else
			{
				// printf("shotnum = %d\n",i);				
				send[0] = 0;
			}
			
			MPI_Send(send,9,MPI_INT,ip,99,MPI_COMM_WORLD);
			nsend = nsend+1;
			// if(i<ntask)printf("Calculating Gradient. Send No.=%d. Shot to Processor %d\n",send[0],ip);
//			fflush(stdout);
		}
      
	MPI_Barrier(MPI_COMM_WORLD);
    }
	}
	else    
	{

    // cudaSetDevice((myid-1)%8);			
    FD2DGPU_ELASTIC image2delastic(model->sou,model->dx,model->dz,model->dt,model->nxpml,model->nzpml,model->allnx,model->allnz,model->scale,model->pml,model->nt,model->nop,model->ntr_pre);	
    image2delastic.GPUbufferVPVS(model->velp,model->vels);
    for(int iter=0;iter<maxiter;iter++)
    {
        misfit =0;

        memset(pp_grad,0,sizeof(float)*model->allnx*model->allnz);        
        memset(ps_grad,0,sizeof(float)*model->allnx*model->allnz); 
        memset(illumination,0,sizeof(float)*model->allnx*model->allnz);                     
     
        memset(b_pp_grad,0,sizeof(float)*model->allnx*model->allnz);        
        memset(b_ps_grad,0,sizeof(float)*model->allnx*model->allnz); 
        memset(b_pp_grad2,0,sizeof(float)*model->allnx*model->allnz);        
        memset(b_ps_grad2,0,sizeof(float)*model->allnx*model->allnz); 

        image2delastic.GPUbufferM(image_pp,image_ps); 
		MPI_Send(send,9,MPI_INT,0,0,MPI_COMM_WORLD);
		for(;;)
		{
			MPI_Recv(recv,9,MPI_INT,0,99,MPI_COMM_WORLD,&status);

            int sno = recv[0];
            int ntr = recv[1];
            int scx = recv[2];
            int scy = recv[3];
            int cx_min_s = recv[4];
            int cy_min_s = recv[5];
            int cx_max_s = recv[6];
            int cy_max_s = recv[7];
            int pos = recv[8];

			if(sno == 0)
			{
				// printf("myid=%d,calculating gradient finished,waiting...\n",myid);	
				break;
			}

            if( sno<model->minshot||sno>model->maxshot ) {
                // printf(" the %dth shot is out of imaging range\n", sno);
                // fflush(stdout);
			    MPI_Send(send,9,MPI_INT,0,myid,MPI_COMM_WORLD); 
                continue;               
            }

			printf("Calculating Gradient,  %d Shot in Processor %d\n",sno,myid);

		    int   ngx_left  = cx_min_s;
			int   ngx_right = cx_max_s;
			int   nx_s = (cx_max_s - cx_min_s)/model->dx + 1 ;
            // printf("all_left=%d,ngx_left=%d,ngx_right=%d\n nx_s=%d,pos=%d,ntr=%d\n",model->all_left,ngx_left,ngx_right,nx_s,pos,ntr);
	        int delta_left =int((ngx_left-model->all_left)/model->dx);

            if(ntr==image2delastic.ntr&&ntr==model->ntr_pre){
                read_shot_gather_su(model->fn1, pos, ntr, model->nt, model->record_z, model->gc);
                read_shot_gather_su(model->fn2, pos, ntr, model->nt, model->record_x, model->gc);             
                image2delastic.record_copytoGPU(model->record_z,model->record_x,model->gc);
            }
            else{
                printf("Trace numbers differ from preload ntr_pre in this shot. Need reallocate.\n");
		        image2delastic.ntr = ntr;
                delete[] model->record_z;
                delete[] model->record_x;
                delete[] model->gc;                
                model->record_z = new float[model->nt*ntr]{};
                model->record_x = new float[model->nt*ntr]{};                
                model->gc = new int[ntr]{};
                read_shot_gather_su2(model->fn1, pos, ntr, model->nt, model->record_z, model->gc);  
                read_shot_gather_su2(model->fn2, pos, ntr, model->nt, model->record_x, model->gc); 				
                image2delastic.record_copytoGPU(model->record_z,model->record_x,model->gc);
		        std::cout<<"reallocate completed"<<std::endl;
            }
            if(nx_s==model->nx){
			    image2delastic.bufferVpVsHtoD(delta_left+record_left_in_v);
            }
            else{
                printf("Nx differ from preload ntr_pre in this shot. Need reallocate.\n");
                //TODO reallocate space and memory copy
			    image2delastic.bufferVpVsHtoD(delta_left+record_left_in_v);
            }

            image2delastic.delta_left = delta_left;
            image2delastic.record_left_in_v = record_left_in_v;                        
	        image2delastic.isx = model->pml+(int)((scx-ngx_left)/model->dx);
            image2delastic.isz = model->pml + sy;
            image2delastic.igz = model->pml + gy;
            image2delastic.cmin = delta_left+record_left_in_v;          //most left position(grid) in velocity/image model 
            image2delastic.ngx_left = ngx_left;                         //most left positon   
			int Dim = model->nxpml*model->nzpml;
            // printf("scx = %d,isx = %d\n",scx,image2delastic.isx);
			// printf("source num=%d\n",sno);
		    // printf("gc1=%d,gc2=%d\n",model->gc[0],model->gc[ntr-1]);

			ELSRTM_ROT_EXP_IMAGE_SINGLESHOT(image2delastic,sno,Dim,myid,model->record_z,model->record_x,model->gc,model->rbc,model->cpu_mem,model->iointerval,iter,model->imagedir,image_pp,image_ps,pp_grad,ps_grad,misfit);

			MPI_Send(send,9,MPI_INT,0,myid,MPI_COMM_WORLD);
		}
        image2delastic.imagebuffer_resettozero(pp_grad,ps_grad,illumination,b_pp_grad,b_ps_grad,b_pp_grad2,b_ps_grad2);             //TODO:bug exist!
		printf("myid=%d,calculating gradient finished,waiting...\n",myid);	        
        MPI_Barrier(groupcomm);
        if(myid==1){
        MPI_Reduce(MPI_IN_PLACE, &pp_grad[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);
        MPI_Reduce(MPI_IN_PLACE, &ps_grad[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);	
        MPI_Reduce(MPI_IN_PLACE, &illumination[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);
        MPI_Reduce(MPI_IN_PLACE, &misfit, 1, MPI_DOUBLE, MPI_SUM, 0, groupcomm);        


        MPI_Reduce(MPI_IN_PLACE, &b_pp_grad[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);
        MPI_Reduce(MPI_IN_PLACE, &b_ps_grad[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);	 

        MPI_Reduce(MPI_IN_PLACE, &b_pp_grad2[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);
        MPI_Reduce(MPI_IN_PLACE, &b_ps_grad2[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);	                
        }
        else{
        MPI_Reduce(&pp_grad[0], &pp_grad[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);
        MPI_Reduce(&ps_grad[0], &ps_grad[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);	
        MPI_Reduce(&illumination[0], &illumination[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);        
        MPI_Reduce(&misfit, &misfit, 1, MPI_DOUBLE, MPI_SUM, 0, groupcomm);           

        MPI_Reduce(&b_pp_grad[0], &b_pp_grad[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);
        MPI_Reduce(&b_ps_grad[0], &b_ps_grad[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);	   

        MPI_Reduce(&b_pp_grad2[0], &b_pp_grad2[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);
        MPI_Reduce(&b_ps_grad2[0], &b_ps_grad2[0], model->allnx*model->allnz, MPI_FLOAT, MPI_SUM, 0, groupcomm);	                
        }
		if(myid==1){

            if(iter==0){
                misfit = 0.0f;
            }

            array_misfit[0][iter] = 1.0f - fabs((1.0f/model->ns)*misfit);
            array_misfit[1][iter] = 1.0f - fabs((1.0f/model->ns)*misfit);       

            std::cout<<"iter = "<<iter<<"  misfit = "<<misfit<<"  rmsmisfit = "<<array_misfit[0][iter]<<std::endl;            

//  Bandpass filter to bgrad & bgrad2
            // int x1 = 100;
            // int x2 = 60;
            // int z1 = 50;            // z1 = x1 * (nz/nx)
            // int z2 = 30;            // z2 = x2 * (nz/nx)         
            // bandpass_filter(model->allnx,model->allnz,x1,x2,z1,z2,b_pp_grad2,b_ps_grad2);
            // bandpass_filter(model->allnx,model->allnz,x1,x2,z1,z2,b_pp_grad,b_ps_grad);      

// new band pass using gaussian windows
            // int x1 = 0;       // lowcutoff_x    distance between origin and low cut off x
            // int x2 = 0;        // highcutoff x   distance between origin and high cut off x
            // int z1 = 0;            // z1 = x1 * (nz/nx)
            // int z2 = 80;            // z2 = x2 * (nz/nx)         
            // bandpass_filter(model->allnx,model->allnz,x1,x2,z1,z2,b_pp_grad2,b_ps_grad2);
            // bandpass_filter(model->allnx,model->allnz,x1,x2,z1,z2,b_pp_grad,b_ps_grad);                      
//  Bandpass filter to bgrad & bgrad2


            if(iter==0){
                alpha_k = alpha;
            }
            else{
                alpha_k = (alpha/2)*pow(gama,iter);
            }

            k = iter + 1;

            for(int i=0;i<model->allnx;i++)
                for(int j=0;j<model->allnz;j++)
                {

                    sum_pp_grad[i*model->allnz+j] = pp_grad[i*model->allnz+j] + b_pp_grad[i*model->allnz+j] + b_pp_grad2[i*model->allnz+j];
                    sum_ps_grad[i*model->allnz+j] = ps_grad[i*model->allnz+j] + b_ps_grad[i*model->allnz+j] + b_ps_grad2[i*model->allnz+j];   

                    s_k_p[i*model->allnz+j] = beta1*s_k0_p[i*model->allnz+j] + (1-beta1)*sum_pp_grad[i*model->allnz+j];
                    s_k_s[i*model->allnz+j] = beta1*s_k0_s[i*model->allnz+j] + (1-beta1)*sum_ps_grad[i*model->allnz+j];

                    r_k_p[i*model->allnz+j] = beta2*r_k0_p[i*model->allnz+j] + (1-beta2)*(sum_pp_grad[i*model->allnz+j]*sum_pp_grad[i*model->allnz+j]);
                    r_k_s[i*model->allnz+j] = beta2*r_k0_s[i*model->allnz+j] + (1-beta2)*(sum_ps_grad[i*model->allnz+j]*sum_ps_grad[i*model->allnz+j]);

                    s_kmod_p[i*model->allnz+j] = s_k_p[i*model->allnz+j]/(1 - pow(beta1,k));                    
                    s_kmod_s[i*model->allnz+j] = s_k_s[i*model->allnz+j]/(1 - pow(beta1,k));                    
                    r_kmod_p[i*model->allnz+j] = r_k_p[i*model->allnz+j]/(1 - pow(beta2,k));   
                    r_kmod_s[i*model->allnz+j] = r_k_s[i*model->allnz+j]/(1 - pow(beta2,k));  

                    delta_mp[i*model->allnz+j] = alpha_k*(s_kmod_p[i*model->allnz+j]/(sqrt(r_kmod_p[i*model->allnz+j]) + epsilon));
                    delta_ms[i*model->allnz+j] = alpha_k*(s_kmod_s[i*model->allnz+j]/(sqrt(r_kmod_s[i*model->allnz+j]) + epsilon));

                    image_pp[i*model->allnz+j] = image_pp[i*model->allnz+j] - delta_mp[i*model->allnz+j];
                    image_ps[i*model->allnz+j] = image_ps[i*model->allnz+j] - delta_ms[i*model->allnz+j];                         

                }

            memcpy(s_k0_p,s_k_p,model->allnx*model->allnz*sizeof(float));
            memcpy(s_k0_s,s_k_s,model->allnx*model->allnz*sizeof(float));
            memcpy(r_k0_p,r_k_p,model->allnx*model->allnz*sizeof(float));
            memcpy(r_k0_s,r_k_s,model->allnx*model->allnz*sizeof(float));



            char imagepath_pp[1024];
            FILE *fpp = NULL;
     
            sprintf(imagepath_pp, "%s/sum_bgrad2_pp_%d_%d_%d.dat",model->imagedir,iter,model->allnx,model->allnz);
            fpp = fopen(imagepath_pp,"wb");
            fwrite(b_pp_grad2,model->allnx*model->allnz*sizeof(float),1,fpp);
            fclose(fpp);  

            sprintf(imagepath_pp, "%s/sum_bgrad2_ps_%d_%d_%d.dat",model->imagedir,iter,model->allnx,model->allnz);
            fpp = fopen(imagepath_pp,"wb");
            fwrite(b_ps_grad2,model->allnx*model->allnz*sizeof(float),1,fpp);
            fclose(fpp);   

            sprintf(imagepath_pp, "%s/all_grad_pp_%d_%d_%d.dat",model->imagedir,iter,model->allnx,model->allnz);
            fpp = fopen(imagepath_pp,"wb");
            fwrite(sum_pp_grad,model->allnx*model->allnz*sizeof(float),1,fpp);
            fclose(fpp);  

            sprintf(imagepath_pp, "%s/all_grad_ps_%d_%d_%d.dat",model->imagedir,iter,model->allnx,model->allnz);
            fpp = fopen(imagepath_pp,"wb");
            fwrite(sum_ps_grad,model->allnx*model->allnz*sizeof(float),1,fpp);
            fclose(fpp);              

            sprintf(imagepath_pp, "%s/sum_image_pp_%d_%d_%d.dat",model->imagedir,iter,model->allnx,model->allnz);
            fpp = fopen(imagepath_pp,"wb");
            fwrite(image_pp,model->allnx*model->allnz*sizeof(float),1,fpp);
            fclose(fpp); 

            sprintf(imagepath_pp, "%s/sum_image_ps_%d_%d_%d.dat",model->imagedir,iter,model->allnx,model->allnz);
            fpp = fopen(imagepath_pp,"wb");
            fwrite(image_ps,model->allnx*model->allnz*sizeof(float),1,fpp);
            fclose(fpp);    


//*******************************  Using as RTM *******************************************
            // exit(0);
//*******************************  Using as RTM *******************************************
        }

	    MPI_Bcast(image_pp, model->allnx*model->allnz, MPI_FLOAT, 0, groupcomm);
	    MPI_Bcast(image_ps, model->allnx*model->allnz, MPI_FLOAT, 0, groupcomm);
	    MPI_Barrier(MPI_COMM_WORLD);
	}
	// MPI_Barrier(MPI_COMM_WORLD);
        if(myid==1){
            char misfitdir[1024];
            sprintf(misfitdir, "%s/misfit.txt",model->imagedir);            
            FILE *fpout = fopen(misfitdir,"w");
            for(int i=0;i<maxiter;i++){
                fprintf(fpout,"i = %d,misfit = %e,rms_misfit = %f\n",i,array_misfit[0][i],array_misfit[1][i]);
            }
            fclose(fpout);
        }

    }

    delete[] s_k_p;
    delete[] s_k_s;
    delete[] r_k_p;
    delete[] r_k_s;
    delete[] s_k0_p;
    delete[] s_k0_s;
    delete[] r_k0_p;
    delete[] r_k0_s;   
    delete[] s_kmod_p;
    delete[] s_kmod_s;
    delete[] r_kmod_p;
    delete[] r_kmod_s;
    delete[] delta_mp;
    delete[] delta_ms;

    delete[] pp_cg;
    delete[] ps_cg;
    delete[] pp_cg_old;
    delete[] ps_cg_old;    

    delete[]b_pp_grad;
    delete[]b_ps_grad;
    delete[]b_pp_grad2;
    delete[]b_ps_grad2;    

    delete[]sum_pp_grad;
    delete[]sum_ps_grad;    
}
//  elsrtm one iteration ///


// CG


void sjvecaddf(float *z, int n, float a, float *x, float b, float *y){
    //! z[] = a*x[] + b*y[]
    int ii;
    for (ii = 0; ii < n; ++ii){
        z[ii] = a * x[ii] + b * y[ii];
    }
}

float sjcgbeta(int n, float *cg, float *g1, float *g0, int iter){

    if (iter == 0) {
        return 0.0f;
    } else {
        int ii;
        float a, b, c;
        a = b = c = 0;
        for (ii = 0; ii < n; ++ii) {
            a += g1[ii] * (g1[ii] - g0[ii]);
            b += cg[ii] * (g1[ii] - g0[ii]);
            c += g1[ii] * g1[ii];
        }
        printf("a=%f,b=%f,c=%f\n",a,b,c);
        printf("fabsf(b)=%f\n",fabsf(b));

        float beta_HS = (fabsf(b) > 0.0f) ? (a / b) : 0.0f;
        float beta_DY = (fabsf(b) > 0.0f) ? (c / b) : 0.0f;

        printf("beta_HS=%f,beta_DY=%f\n",beta_HS,beta_DY);        
        return MAX(0.0f, MIN(beta_HS, beta_DY));
    }
}

int sjfindabsmaxf(float *a, int n){
    //! Find the index of a vector' absolute maxmium
    int ii, index = 0;
    for (ii = 0; ii < n; ++ii){
        if (fabs(a[ii]) > fabs(a[index])) index = ii;
    }
    return index;
}

float sjcglength(int n, float *s, float *x, float err, int iter){

    int index_s, index_x;
    if (iter == 0) {
      return 0.0001f;
//       return 10.0f;        
    } else {
        index_s = sjfindabsmaxf(s, n);
        index_x = sjfindabsmaxf(x, n);
        return err * fabsf(x[index_x]) / fabsf(s[index_s]);
    }
}

void sjcgsolver(float *m, int n, float *cg, float *g1, float *g0, int iter){

    float alpha, beta;

    //! Calculate CG direction
    if(iter==0) {
        beta = 0.0f;
        sjvecaddf(cg, n, -1.0f, g1, beta, cg);
    } else {
        beta = sjcgbeta(n, cg, g1, g0, iter);              // calculate beta.      use g0 g1 cg.
        sjvecaddf(cg, n, -1.0f, g1, beta, cg);          //calculate conjugate gradient.     cg = -g1 + beta*cg;
    }
    memcpy(g0, g1, n * sizeof(float));
    //! Calculate CG step size
    alpha = sjcglength(n, cg, m, 0.1f, iter);    // old is 0.1
    //! Update model
    printf("iter=%d,alpha=%f,beta=%lf\n",iter,alpha,beta);

    // if(iter==0){
    //     sjvecaddf(m,n,1.0f,m,1.0f,g1);
    // }
    // else{
         sjvecaddf(m, n, 1.0f, m, alpha, cg);    //m = m + alpha*cg;
    // }
//     sjvecaddf(m, n, 1.0f, m, 100, g1);    //m = m + alpha*cg;    
}

// CG





void gradientmethod(float *m, int n, float *g1, float alpha, int iter){

    if(iter==0){
        sjvecaddf(m,n,1.0f,m,1.0f,g1);
    }
    else{
       sjvecaddf(m, n, 1.0f, m, -alpha, g1);    //m = m + alpha*cg;
    }
   
}




////  main function of reverse time migration /////
int main(int argc,char *argv[]){

	char parfn[1024];
    int i,j,isx,isz;
	strcpy(parfn,argv[1]);
	modelpar model;
	char fvelp[1024];
	char fvels[1024];    
	// char fn1[1024];
	// char fn2[1024];    
	float dx,dz;
	int pml,allnx,allnz,nx,nz,dis_shot,scale,order;
	float t;
	int mode;
    int nop;
    int minshot,maxshot;
    int sy,gy;
	float f0;
	int maxiter;
	bool light,rbc,cpu_mem;
	int light_temp,rbc_temp,cpumem_temp,iointerval;
    scale =1 ;     
	FILE *fp = NULL;
	fp = fopen(parfn,"r");
	fscanf(fp,"%f %f %d %d %d %d %d %d %d",&dx,&dz,&pml,&allnx,&allnz,&minshot,&maxshot);       //TODO:  pml exist bug, only 100 multiply work
	fscanf(fp,"%s",fvelp);
	fscanf(fp,"%s",fvels);    
	fscanf(fp,"%s",model.fn1);
	fscanf(fp,"%s",model.fn2);
	fscanf(fp,"%s",model.imagedir);     
	fscanf(fp,"%f %d",&f0,&maxiter);
	fscanf(fp,"%d %d",&sy,&gy);    
	fscanf(fp,"%d %d %d %d",&light_temp,&rbc_temp,&cpumem_temp,&iointerval);  	   	
	fclose(fp);

	int myid,np;
	MPI_Status status;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&np);

	light = light_temp;
	rbc = rbc_temp;
	cpu_mem = cpumem_temp;
	
	if(myid==0){
        if(access(model.imagedir, F_OK) == 0){
            std::cout<<"image output directory exist."<<std::endl;
        }
        else{
            mkdir(model.imagedir,0755);
            std::cout<<"mkdir image directory."<<std::endl;
        }      
		printf("minshot=%d,maxshot=%d\n",minshot,maxshot);
	}

	float *image_pp = new float[allnx*allnz]{};
	float *image_ps = new float[allnx*allnz]{};

	float *image_pp_m = new float[allnx*allnz]{};
	float *image_ps_m = new float[allnx*allnz]{};

	float *pp_grad = new float[allnx*allnz]{};
	float *ps_grad = new float[allnx*allnz]{};

    float *illumination;
    illumination = new float[allnx*allnz]{};

	float idz = 1.0f/dz;
	float idx = 1.0f/dx;

	int nt,ns;
	float dt;
	float x0 = 0.0f;
    float cmin,cmax,cmleft,cmright;
	cmin = x0;
    cmax = x0 + (allnx-1)*dx; 

    int **table= NULL ;   //
    table  = alloc2int( 9, 100000 );
    int all_left = 999999999;

    if(myid ==0){
        if(index_shot(model.fn1, &nt, &dt, &ns, table)){
             printf("Can not read the shot file!\n");
             return 0;
        }

        for(i=0;i<ns;i++){
        //	std::cout<<table[i][4]<<std::endl;
            all_left = MIN(all_left,table[i][4]);
        }        
        //edges of Imaging
        cmleft  = 999999;
        cmright = -999999;
        for (i=0;i<ns;i++){
            if(cmleft>table[i][4])cmleft=table[i][4];
            if(cmright<table[i][6])cmright=table[i][6];
            if(cmin>table[i][4])printf(" Warning! The %dth Shot's minimum coordinate is on the left of velocity model\n",i);
            if(cmax<table[i][6])printf(" Warning! The %dth Shot's maximum coordinate is on the right of velocity model\n",i);
        }
    }
    if(myid==0){
          printf("==========Parameters of input seismic file============\n");
          printf(" Shot number           : %d\n",ns);
          printf(" Sampling point number : %d\n",nt);
          printf(" Sampling interval     : %f s\n",dt);
          printf("======================================================\n\n");
        //   for(i=0;i<ns;i++)printf(" table[%d][8] : %d\n",i,table[i][8]);
    }
    int ntr_pre;
    if(myid==0){
        ntr_pre = table[0][1];
        nx = (int)(table[0][6]/dx+0.5) - (int)(table[0][4]/dx+0.5) + 1;
        nz = allnz;
        printf("nx=%d,nz=%d\n",nx,nz);     
    }
    MPI_Bcast(&all_left, 1, MPI_INT, 0, MPI_COMM_WORLD);    
    MPI_Bcast(&ns, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nt, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dt, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ntr_pre, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nx, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nz, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);    
	MPI_Barrier(MPI_COMM_WORLD);  
	int nzpml = nz + 2*pml;
	int nxpml = nx + 2*pml;

	model.velp = NULL;	
	model.vels = NULL;	
	model.velp = alloc1float(allnz*allnx);
	model.vels = alloc1float(allnz*allnx);

	FILE *fp1 = NULL;
    FILE *fp2 = NULL;
	fp1 = fopen(fvelp,"rb");
	fp2 = fopen(fvels,"rb");

	for(i=0; i<allnx; i++)
	{
		fread(&model.velp[i*allnz],sizeof(float),allnz,fp1);
		fread(&model.vels[i*allnz],sizeof(float),allnz,fp2);        
	}
	fclose(fp1);
	fclose(fp2);

	model.sou = NULL;
	model.sou = alloc1float(nt);

//	float f0 = 25;
	ricker1(nt,f0,dt,model.sou);

    // cudaMallocHost((void**)&(model.record_z),nt*ntr_pre*sizeof(float));
    // cudaMallocHost((void**)&(model.record_x),nt*ntr_pre*sizeof(float));
    // cudaMemset(model.record_z,0,nt*ntr_pre*sizeof(float));
    // cudaMemset(model.record_x,0,nt*ntr_pre*sizeof(float));

    model.record_z = new float[nt*ntr_pre];
    model.record_x = new float[nt*ntr_pre];
	memset(model.record_z,0,sizeof(float)*nt*ntr_pre);
	memset(model.record_x,0,sizeof(float)*nt*ntr_pre);

    model.gc = new int[ntr_pre];

    nop = 4;

    // FILE* fimagepp = NULL;
    // fimagepp = fopen("/home/rondo/ertm-cuda-del/velmodel/refl_no_prismatic_vp_601_201.bin","rb");
    // fread(image_pp,sizeof(float),allnx*allnz,fimagepp);
    // fclose(fimagepp);
    // fimagepp = fopen("/home/rondo/ertm-cuda-del/velmodel/refl_no_prismatic_vs_601_201.bin","rb");
    // fread(image_ps,sizeof(float),allnx*allnz,fimagepp);
    // fclose(fimagepp);

    int coordinate_scale = 1;
    int record_left_in_v = 0;    
	MPI_Barrier(MPI_COMM_WORLD);
	init_modelparameters(&model,dx,dz,dt,minshot,maxshot,nx,nz,ns,nxpml,nzpml,allnx,allnz,scale,pml,nt,nop,ntr_pre,light,rbc,cpu_mem,iointerval,all_left);
    lsrtm_all_iteration(myid,np,sy,gy,status,table,&model,image_pp,image_ps,pp_grad,ps_grad,illumination,maxiter,record_left_in_v);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	free1float(model.velp);
	free1float(model.vels);    
	free1float(model.sou);
	free2int(table);

    delete[] model.record_z;
    delete[] model.record_x; 
	// cudaFreeHost(model.record_z);
	// cudaFreeHost(model.record_x);
	   
    delete[] model.gc;
    delete[] image_pp;
    delete[] image_ps;  
    delete[] image_pp_m;
    delete[] image_ps_m;          
    delete[] pp_grad;
    delete[] ps_grad;
	return 0;

}











