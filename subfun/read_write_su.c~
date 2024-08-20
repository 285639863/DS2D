/**********************************************************
*     Author:Duan Xinbiao,SGRI.   2011.8.9-13             *
**********************************************************/
#include <stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include <math.h>
#include "segy.h"

bhed   Bh;
segy1  Th;

void write_shot_gather_su_3d(char *path, int sno, int *sx, int *sy, int **gc_x, int **gc_y, int nx, int ny,  int nt , float dt, float ****seis)
{
    int i,j,k;
    short  ns;
    short  hdt;
    char   path0[1024], tmps[1024], filenm[1024];
    FILE   *fpo;
    
    ns = nt ;
    hdt= dt*1000000;
/*
    strcpy(path0, path);
    sprintf(tmps, "/directshot_%d.su", sno);
    strcat(path0, tmps);
    strcpy(filenm, path0);
    fpo = fopen(filenm,"wb");
*/
    fpo = fopen(path,"wb");


    if(fpo == NULL) {
        printf("Sorry,cann't open output data file!\n");
        exit(0);
    }
for(j=0;j<sno;j++)
  for(i=0; i<ny; i++)
    for(k=0; k<nx; k++)
    {
       // trace header
       Th.tracl =i;         // trace sequence number within line //
//     Th.tracr =i;  //5-8     // trace sequence number within reel //
       Th.fldr  =j+1;           //9-12    // field record number //
       Th.tracf =i*nx+k;  //13-16   // trace number within field record //     
       Th.ep    =j+1;  //17-20   // energy source point number //
       Th.cdp   =k;  //21-24
//     Th.offset=0;          //37-40
//     Th.gelev =gc[i][0];           //41-44
//       Th.selev =sz;           //45-48
       Th.sx    = sx[j];           //73-76
       Th.sy    = sy[j];           //77-80
       Th.gx    = gc_x[j][k];
//       Th.gx    =gc[i][1] ;  //81-84
       Th.gy    = gc_y[j][i];           //85-88
       Th.offset=Th.gx-Th.sx;
       Th.ns    =ns;          //
       Th.dt    =hdt;         //    
       Th.trid  =1;           //    // trace identification code:
       Th.scalel=1 ;
       Th.scalco=1 ;
   
       fwrite(&Th,240,1,fpo);
       fwrite(&seis[j][k][i][0],sizeof(float),ns,fpo);
    }

    fclose(fpo);
}











void write_shot_gather_su(char *path, int sno, int *sx, int sz, int **gc, int ntr,  int nt , float dt, float ***seis)
{
    int i,j;
    short  ns;
    short  hdt;
    char   path0[1024], tmps[1024], filenm[1024];
    FILE   *fpo;
    
    ns = nt ;
    hdt= dt*1000000;
/*
    strcpy(path0, path);
    sprintf(tmps, "/directshot_%d.su", sno);
    strcat(path0, tmps);
    strcpy(filenm, path0);
    fpo = fopen(filenm,"wb");
*/
    fpo = fopen(path,"wb");


    if(fpo == NULL) {
        printf("Sorry,cann't open output data file!\n");
        exit(0);
    }
for(j=0;j<sno;j++)
    for( i=0; i<ntr; i++)
    {
       // trace header
       Th.tracl =i;         // trace sequence number within line //
       Th.tracr =i;  //5-8     // trace sequence number within reel //
       Th.fldr  =j+1;           //9-12    // field record number //
       Th.tracf =i;  //13-16   // trace number within field record //     
       Th.ep    =i;  //17-20   // energy source point number //
       Th.cdp   =0;  //21-24
//     Th.offset=0;          //37-40
//     Th.gelev =gc[i][0];           //41-44
       Th.selev =sz;           //45-48
       Th.sx    =sx[j];           //73-76
       Th.sy    =0;           //77-80
       Th.gx = gc[j][i];
//       Th.gx    =gc[i][1] ;  //81-84
       Th.gy    =0;           //85-88
       Th.offset=Th.gx-Th.sx;
       Th.ns    =ns;          //
       Th.dt    =hdt;         //    
       Th.trid  =1;           //    // trace identification code:
       Th.scalel=1 ;
       Th.scalco=1 ;
   
       fwrite(&Th,240,1,fpo);
       fwrite(&seis[j][i][0],sizeof(float),ns,fpo);
    }

    fclose(fpo);
}

void write_single_shot_su(char *path, int sno, int sx, int sz, int **gc, int ntr,  int nt , float dt, float **seis)
{
	int i,j;
    short  ns;
    short  hdt;
    char   path0[256], tmps[256], filenm[256];
    FILE   *fpo;
    
    ns = nt ;
    hdt= dt*1000000;

    strcpy(path0, path);
    sprintf(tmps, "/shot_%d.su", sno);
    strcat(path0, tmps);
    strcpy(filenm, path0);
    
    fpo = fopen(filenm,"wb");
    if(fpo == NULL) {
        printf("Sorry,cann't open output data file!\n");
        exit(0);
    }

    for( i=0; i<ntr; i++)
    {
       // trace header
       Th.tracl =i;         /* trace sequence number within line */
       Th.tracr =i;  //5-8     /* trace sequence number within reel */
       Th.fldr  =sno;           //9-12    /* field record number */
       Th.tracf =i;  //13-16   /* trace number within field record */     
       Th.ep    =i;  //17-20   /* energy source point number */
       Th.cdp   =0;  //21-24
       //Th.offset=0;          //37-40
       Th.gelev =gc[i][0];           //41-44
       Th.selev =sz;           //45-48
       Th.sx    =sx;           //73-76
       Th.sy    =0;           //77-80
       Th.gx    =gc[i][1] ;  //81-84
       Th.gy    =0;           //85-88
       Th.offset=Th.gx-Th.sx;
       Th.ns    =ns;          //
       Th.dt    =hdt;         //    
       Th.trid  =1;           //    /* trace identification code:
       Th.scalel=1 ;
       Th.scalco=1 ;
   
       fwrite(&Th,240,1,fpo);
       fwrite(&seis[i][0],sizeof(float),ns,fpo);
    }

    fclose(fpo);
}





void write_single_shot_su_extract(char *path, int sno, int sx, int sz, int **gc, int ntr,  int nt , float dt, float **seis, int nextr)
{
	int i,j;
    short  ns;
    short  hdt;
    char   path0[256], tmps[256], filenm[256];
    FILE   *fpo;
    
    ns = nt / nextr;
    dt = dt * nextr;

    hdt= dt*1000000;

    strcpy(path0, path);
    sprintf(tmps, "/shot_%d.su", sno);
    strcat(path0, tmps);
    strcpy(filenm, path0);
    
    fpo = fopen(filenm,"wb");
    if(fpo == NULL) {
        printf("Sorry,cann't open output data file!\n");
        exit(0);
    }

    for(i=0; i<ntr; i++)
    {
       // trace header
       Th.tracl =i;         /* trace sequence number within line */
       Th.tracr =i;  //5-8     /* trace sequence number within reel */
       Th.fldr  =sno;           //9-12    /* field record number */
       Th.tracf =i;  //13-16   /* trace number within field record */     
       Th.ep    =i;  //17-20   /* energy source point number */
       Th.cdp   =0;  //21-24
       //Th.offset=0;          //37-40
       Th.gelev =gc[i][0];           //41-44
       Th.selev =sz;           //45-48
       Th.sx    =sx;           //73-76
       Th.sy    =0;           //77-80
       Th.gx    =gc[i][1] ;  //81-84
       Th.gy    =0;           //85-88
       Th.offset=Th.gx-Th.sx;
       Th.ns    =ns;          //
       Th.dt    =hdt;         //    
       Th.trid  =1;           //    /* trace identification code:
       Th.scalel=1 ;
       Th.scalco=1 ;
   
       fwrite(&Th,240,1,fpo);
       for(j=0 ;j<ns ;j++)seis[i][j]=seis[i][nextr*j];
       fwrite(&seis[i][0],sizeof(float),ns,fpo);
    }

    fclose(fpo);

}


void   read_shot_gather_su(char *fn, int pos, int ntr, int nt, float **dat,int **gc)
{
	int i;
    FILE   *fp ;
     
    fp = fopen(fn,"rb");
    if(fp == NULL) {
            printf("Sorry,cann't open input seismic file!\n");
            exit(0);
    }
          
    int TL = 240 + nt*sizeof(float);
           
    fseek(fp, TL*pos, SEEK_SET);

    for(i=0; i<ntr; i++){
           fread(&Th, 240, 1, fp);
           gc[i][0]=Th.gelev;
           gc[i][1]=Th.gx;
           gc[i][2]=Th.gy;
               
           fread(&dat[i][0], sizeof(float), nt, fp);
    }

    fclose(fp);
}

void   read_shot_gather_coordinate_su(char *fn, int pos, int ntr, int nt, int **gc)
{
	int i;
    FILE   *fp ;
     
    fp = fopen(fn,"rb");
    if(fp == NULL) {
            printf("Sorry,cann't open input seismic file!\n");
            exit(0);
    }
          
    int TL = 240 + nt*sizeof(float);
           
    fseek(fp, TL*pos, SEEK_SET);

    for( i=0; i<ntr; i++){
           fread(&Th, 240, 1, fp);
           gc[i][0]=Th.gelev;
           gc[i][1]=Th.gx ;
           gc[i][2]=Th.gy;
               
           //fread(&dat[i][0], sizeof(float), nt, fp);
           fseek(fp, nt*sizeof(float), SEEK_CUR);
    }

    fclose(fp);
}

void   write_image_single_shot_su(char *path, int sno, int nx, int nz, float dx, float dz, int cmin,float **image)
{

	int i;
    short  ns;
    short  hdt, dt;
    char   path0[256], tmps[256], filenm[256];
    FILE   *fpo;

    ns = nz;
    dt = dz;

    hdt= dt*1000;

    strcpy(path0, path);
    sprintf(tmps, "/mig_of_shot_%d.su", sno);
    strcat(path0, tmps);
    strcpy(filenm, path0);
    
    fpo = fopen(filenm,"wb");
    if(fpo == NULL) {
        printf("Sorry,cann't open output data file!\n");
        exit(0);
    }


    for( i=0; i<nx; i++)
    {
       // trace header
       Th.tracl =i;         /* trace sequence number within line */
       Th.tracr =i;  //5-8     /* trace sequence number within reel */
       Th.fldr  =sno;           //9-12    /* field record number */
       Th.tracf =i;  //13-16   /* trace number within field record */     
       Th.ep    =sno;  //17-20   /* energy source point number */
       Th.cdp   =i;  //21-24
       //Th.offset=0;          //37-40
       Th.gelev =0;           //41-44
       Th.selev =0;           //45-48
       Th.sx    =i*dx+cmin ;  //73-76
       Th.sy    =0;           //77-80
       Th.gx    =i*dx+cmin ;  //81-84
       Th.gy    =0;           //85-88
       Th.offset=0;
       Th.ns    =ns;          //
       Th.dt    =hdt;         //    
       Th.trid  =1;           //    /* trace identification code:
       Th.scalel=1 ;
       Th.scalco=1 ;

       fwrite(&Th,240,1,fpo);
       fwrite(&image[i][0],sizeof(float),ns,fpo);
    }

    fclose(fpo);
}

void  write_stack_image_all_shot_su(char *path, int nx, int nz, float dx, float dz, int cmin,float **image)
{

	int i;
    short  ns;
    short  hdt, dt;
    char   path0[256], tmps[256], filenm[256];
    FILE   *fpo;

    ns = nz;
    dt = dz;

    hdt= dt*1000;

    strcpy(path0, path);
    sprintf(tmps, "/stack_image_of_all_shot.su");
    strcat(path0, tmps);
    strcpy(filenm, path0);
    
    fpo = fopen(filenm,"wb");
    if(fpo == NULL) {
        printf("Sorry,cann't open output data file!\n");
        exit(0);
    }


    for( i=0; i<nx; i++)
    {
       // trace header
       Th.tracl =i;         /* trace sequence number within line */
       Th.tracr =i;  //5-8     /* trace sequence number within reel */
       Th.fldr  =1;           //9-12    /* field record number */
       Th.tracf =i;  //13-16   /* trace number within field record */     
       Th.ep    =1;  //17-20   /* energy source point number */
       Th.cdp   =i;  //21-24
       //Th.offset=0;          //37-40
       Th.gelev =0;           //41-44
       Th.selev =0;           //45-48
       Th.sx    =i*dx+cmin ;  //73-76
       Th.sy    =0;           //77-80
       Th.gx    =i*dx+cmin ;  //81-84
       Th.gy    =0;           //85-88
       Th.offset=0;
       Th.ns    =ns;          //
       Th.dt    =hdt;         //    
       Th.trid  =1;           //    /* trace identification code:
       Th.scalel=1 ;
       Th.scalco=1 ;

       fwrite(&Th,240,1,fpo);
       fwrite(&image[i][0],sizeof(float),ns,fpo);
    }

    fclose(fpo);
}


void   write_illumination_single_shot_su(char *path, int sno, int nx, int nz, float dx, float dz, int cmin,float **illum)
{

	int i;
    short  ns;
    short  hdt, dt;
    char   path0[256], tmps[256], filenm[256];
    FILE   *fpo;

    ns = nz;
    dt = dz;

    hdt= dt*1000;

    strcpy(path0, path);
    sprintf(tmps, "/illumination_of_shot_%d.su", sno);
    strcat(path0, tmps);
    strcpy(filenm, path0);
    
    fpo = fopen(filenm,"wb");
    if(fpo == NULL) {
        printf("Sorry,cann't open output data file!\n");
        exit(0);
    }


    for( i=0; i<nx; i++)
    {
       // trace header
       Th.tracl =i;         /* trace sequence number within line */
       Th.tracr =i;  //5-8     /* trace sequence number within reel */
       Th.fldr  =sno;           //9-12    /* field record number */
       Th.tracf =i;  //13-16   /* trace number within field record */     
       Th.ep    =sno;  //17-20   /* energy source point number */
       Th.cdp   =i;  //21-24
       //Th.offset=0;          //37-40
       Th.gelev =0;           //41-44
       Th.selev =0;           //45-48
       Th.sx    =i*dx+cmin ;  //73-76
       Th.sy    =0;           //77-80
       Th.gx    =i*dx+cmin ;  //81-84
       Th.gy    =0;           //85-88
       Th.offset=0;
       Th.ns    =ns;          //
       Th.dt    =hdt;         //    
       Th.trid  =1;           //    /* trace identification code:
       Th.scalel=1 ;
       Th.scalco=1 ;

       fwrite(&Th,240,1,fpo);
       fwrite(&illum[i][0],sizeof(float),ns,fpo);
    }

    fclose(fpo);
}



void   write_grad_single_shot_su(char *path, int sno,  int nx, int nz, float dx, float dz, int cmin,float **illum)
{

	int i;
    short  ns;
    short  hdt, dt;
    char   path0[256], tmps[256], filenm[256];
    FILE   *fpo;

    ns = nz;
    dt = dz;

    hdt= dt*1000;

    strcpy(path0, path);
    sprintf(tmps, "/grad_of_shot_%d.su", sno);
    strcat(path0, tmps);
    strcpy(filenm, path0);
    
    fpo = fopen(filenm,"wb");
    if(fpo == NULL) {
        printf("Sorry,cann't open output data file!\n");
        exit(0);
    }


    for( i=0; i<nx; i++)
    {
       // trace header
       Th.tracl =i;         /* trace sequence number within line */
       Th.tracr =i;  //5-8     /* trace sequence number within reel */
       Th.fldr  =sno;           //9-12    /* field record number */
       Th.tracf =i;  //13-16   /* trace number within field record */     
       Th.ep    =sno;  //17-20   /* energy source point number */
       Th.cdp   =i;  //21-24
       //Th.offset=0;          //37-40
       Th.gelev =0;           //41-44
       Th.selev =0;           //45-48
       Th.sx    =i*dx+cmin ;  //73-76
       Th.sy    =0;           //77-80
       Th.gx    =i*dx+cmin ;  //81-84
       Th.gy    =0;           //85-88
       Th.offset=0;
       Th.ns    =ns;          //
       Th.dt    =hdt;         //    
       Th.trid  =1;           //    /* trace identification code:
       Th.scalel=1 ;
       Th.scalco=1 ;

       fwrite(&Th,240,1,fpo);
       fwrite(&illum[i][0],sizeof(float),ns,fpo);
    }

    fclose(fpo);
}


void   read_syn_shot_gather_su(char *path, int is, int pos, int ntr, int nt, float **dat,int **gc)
{
	int i;
    char   path0[256], tmps[256], fn[256];
    FILE   *fp ;

    strcpy(path0, path);
    sprintf(tmps, "/shot_%d.su",is);
    strcat(path0, tmps);
    strcpy(fn, path0);
    
    fp = fopen(fn,"rb");
    if(fp == NULL) {
            printf("Sorry,cann't open input seismic file!\n");
            exit(0);
    }
          
    int TL = 240 + nt*sizeof(float);
           
    fseek(fp, TL*pos, SEEK_SET);

    for( i=0; i<ntr; i++){
           fread(&Th, 240, 1, fp);
           gc[i][0]=Th.gelev;
           gc[i][1]=Th.gx ;
           gc[i][2]=Th.gy;
               
           fread(&dat[i][0], sizeof(float), nt, fp);
    }

    fclose(fp);
}
