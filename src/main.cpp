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

#include "./subfun/alloc.c"

#include "./subfun/segy.h"

#include<omp.h>
#include "mpi.h"

#include "FD2DGPU.cuh"

#include <algorithm>


void normalize(float* matrix, int length, int tracenum) {
	float max = 0;
	for (int i = 0; i < length * tracenum; i++)
		max = (fabs(matrix[i]) > max) ? fabs(matrix[i]) : max;
	//normalized	
	for (int i = 0; i < length * tracenum; i++) {
		matrix[i] = matrix[i] / max;
	}
}


void normalize_all_rho(float* matrix1,float* matrix2, float* matrix3, int length, int tracenum) {
	float max = 0;
	for (int i = 0; i < length * tracenum; i++){
		max = (fabs(matrix1[i]) > max) ? fabs(matrix1[i]) : max;
    }
    for (int i = 0; i < length * tracenum; i++){
		max = (fabs(matrix2[i]) > max) ? fabs(matrix2[i]) : max;
    }
    for (int i = 0; i < length * tracenum; i++){
		max = (fabs(matrix3[i]) > max) ? fabs(matrix3[i]) : max;
    }    
	//normalized	
	for (int i = 0; i < length * tracenum; i++) {
		matrix1[i] = matrix1[i] / max;
		matrix2[i] = matrix2[i] / max;  
		matrix3[i] = matrix3[i] / max;                
	}
}


void normalize_all(float* matrix1,float* matrix2, int length, int tracenum) {
	float max = 0;
	for (int i = 0; i < length * tracenum; i++){
		max = (fabs(matrix1[i]) > max) ? fabs(matrix1[i]) : max;
    }
    for (int i = 0; i < length * tracenum; i++){
		max = (fabs(matrix2[i]) > max) ? fabs(matrix2[i]) : max;
    }
	//normalized	
	for (int i = 0; i < length * tracenum; i++) {
		matrix1[i] = matrix1[i] / max;
		matrix2[i] = matrix2[i] / max;              
	}
}







int    index_shot(char *fn, int *nt, float *dt, int *ns, int **table,int coordinate_scale=1)    //read seismic record
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

              ntr = 1;

              sx0 = sx;
              sy0 = sy;
             
              cx_min_s = 99999999 ;
              cx_max_s = -999999 ;
          }else{
              ntr ++ ; 
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

// enter interation loop
void lsrtm_all_iteration(int myid,int np,int sy,int gy,MPI_Status status,int **table,modelpar *model,float *image_pp,float *image_ps,float *pp_grad,float *ps_grad,float *illumination,int maxiter,int record_left_in_v){


    float *b_pp_grad = new float[model->allnx*model->allnz];
    float *b_ps_grad = new float[model->allnx*model->allnz];
    float *b_pp_grad2 = new float[model->allnx*model->allnz];
    float *b_ps_grad2 = new float[model->allnx*model->allnz];

    float *sum_pp_grad = new float[model->allnx*model->allnz];
    float *sum_ps_grad = new float[model->allnx*model->allnz];
    if(myid!=0){
        cudaSetDevice((myid-1)%8);			
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
    float beta2 = 0.9;

    float alpha_1 = 2e-3;
    float alpha = 2e-3;

    float epsilon = 1e-2;   

    float gama = 0.99;

    int cut_layer = 20;

    bool norm_all = true;

    if(myid==0){
        char adam_para[1024];
        sprintf(adam_para, "%s/parameters.txt",model->imagedir);            
        FILE *fpout = fopen(adam_para,"w");   
        fprintf(fpout,"beta1 = %e, beta2 = %e\n",beta1,beta2);
        fprintf(fpout,"alpha_1 = %e, alpha = %e\n",alpha_1,alpha);    
        fprintf(fpout,"epsilon = %e\n",epsilon);                
        fprintf(fpout,"gama = %e\n",gama);      
        fprintf(fpout,"cut_layer = %d\n",cut_layer);        
        fprintf(fpout,"norm_all = %d\n",norm_all);                           
        fclose(fpout);
    }


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
				break;
			}

            if( sno<model->minshot||sno>model->maxshot ) {

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

//enter single shot imaging process
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


            char misfitdir[1024];
            sprintf(misfitdir, "%s/misfit.txt",model->imagedir);            
            FILE *fpout = fopen(misfitdir,"a");
            fprintf(fpout,"i = %d,misfit = %e,rms_misfit = %f\n",iter,array_misfit[0][iter],array_misfit[1][iter]);
            fclose(fpout);

            for(int i=0;i<model->allnx;i++)
                for(int j=0;j<cut_layer;j++)
                {
                    pp_grad[i*model->allnz+j] = 0;
                    b_pp_grad[i*model->allnz+j] = 0;
                    b_pp_grad2[i*model->allnz+j] = 0;
                    ps_grad[i*model->allnz+j] = 0;
                    b_ps_grad[i*model->allnz+j] =0;
                    b_ps_grad2[i*model->allnz+j] =0;   
                }


            for(int i=0;i<model->allnx;i++)
                for(int j=0;j<model->allnz;j++)
                {

                    sum_pp_grad[i*model->allnz+j] = pp_grad[i*model->allnz+j] + b_pp_grad[i*model->allnz+j] + b_pp_grad2[i*model->allnz+j];
                    sum_ps_grad[i*model->allnz+j] = ps_grad[i*model->allnz+j] + b_ps_grad[i*model->allnz+j] + b_ps_grad2[i*model->allnz+j];   

                }



            if(norm_all){
                normalize_all(sum_pp_grad, sum_ps_grad, model->allnz, model->allnx);     
            }
            else{
                normalize(sum_pp_grad, model->allnz, model->allnx);         
                normalize(sum_ps_grad, model->allnz, model->allnx);     
            }

            if(iter==0){
                alpha_k = alpha_1;
            }
            else{
                alpha_k = alpha*pow(gama,iter);                
            }

            k = iter;

            for(int i=0;i<model->allnx;i++)
                for(int j=0;j<model->allnz;j++)
                {

                    if(iter==0){
                        s_k_p[i*model->allnz+j] = sum_pp_grad[i*model->allnz+j];
                        s_k_s[i*model->allnz+j] = sum_ps_grad[i*model->allnz+j];

                        r_k_p[i*model->allnz+j] = (sum_pp_grad[i*model->allnz+j]*sum_pp_grad[i*model->allnz+j]);
                        r_k_s[i*model->allnz+j] = (sum_ps_grad[i*model->allnz+j]*sum_ps_grad[i*model->allnz+j]);

                        s_kmod_p[i*model->allnz+j] = s_k_p[i*model->allnz+j];                    
                        s_kmod_s[i*model->allnz+j] = s_k_s[i*model->allnz+j];                    
                        r_kmod_p[i*model->allnz+j] = r_k_p[i*model->allnz+j];   
                        r_kmod_s[i*model->allnz+j] = r_k_s[i*model->allnz+j];                          

                        delta_mp[i*model->allnz+j] = alpha_k*(s_kmod_p[i*model->allnz+j]/(sqrt(r_kmod_p[i*model->allnz+j]) + epsilon));
                        delta_ms[i*model->allnz+j] = alpha_k*(s_kmod_s[i*model->allnz+j]/(sqrt(r_kmod_s[i*model->allnz+j]) + epsilon));

                    }
                    else{
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
                    }

                    image_pp[i*model->allnz+j] = image_pp[i*model->allnz+j] - delta_mp[i*model->allnz+j];       //update image
                    image_ps[i*model->allnz+j] = image_ps[i*model->allnz+j] - delta_ms[i*model->allnz+j];                         

                }

            memcpy(s_k0_p,s_k_p,model->allnx*model->allnz*sizeof(float));
            memcpy(s_k0_s,s_k_s,model->allnx*model->allnz*sizeof(float));
            memcpy(r_k0_p,r_k_p,model->allnx*model->allnz*sizeof(float));
            memcpy(r_k0_s,r_k_s,model->allnx*model->allnz*sizeof(float));


            char imagepath_pp[1024];
            FILE *fpp = NULL;

            sprintf(imagepath_pp, "%s/sum_grad_pp_%d_%d_%d.dat",model->imagedir,iter,model->allnx,model->allnz);
            fpp = fopen(imagepath_pp,"wb");
            fwrite(pp_grad,model->allnx*model->allnz*sizeof(float),1,fpp);
            fclose(fpp);  

            sprintf(imagepath_pp, "%s/sum_grad_ps_%d_%d_%d.dat",model->imagedir,iter,model->allnx,model->allnz);
            fpp = fopen(imagepath_pp,"wb");
            fwrite(ps_grad,model->allnx*model->allnz*sizeof(float),1,fpp);
            fclose(fpp);   


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

        }

	    MPI_Bcast(image_pp, model->allnx*model->allnz, MPI_FLOAT, 0, groupcomm);
	    MPI_Bcast(image_ps, model->allnx*model->allnz, MPI_FLOAT, 0, groupcomm);
	    MPI_Barrier(MPI_COMM_WORLD);
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



int main(int argc,char *argv[]){

	char parfn[1024];
    int i,j,isx,isz;
	strcpy(parfn,argv[1]);
	modelpar model;
	char fvelp[1024];
	char fvels[1024];      
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
	fscanf(fp,"%f %f %d %d %d %d %d %d %d",&dx,&dz,&pml,&allnx,&allnz,&minshot,&maxshot);      
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

	ricker1(nt,f0,dt,model.sou);

    model.record_z = new float[nt*ntr_pre];
    model.record_x = new float[nt*ntr_pre];
	memset(model.record_z,0,sizeof(float)*nt*ntr_pre);
	memset(model.record_x,0,sizeof(float)*nt*ntr_pre);

    model.gc = new int[ntr_pre];

    nop = 4;

    int coordinate_scale = 1;
    int record_left_in_v = 0;    
	MPI_Barrier(MPI_COMM_WORLD);
	init_modelparameters(&model,dx,dz,dt,minshot,maxshot,nx,nz,ns,nxpml,nzpml,allnx,allnz,scale,pml,nt,nop,ntr_pre,light,rbc,cpu_mem,iointerval,all_left);  //init parameters
    lsrtm_all_iteration(myid,np,sy,gy,status,table,&model,image_pp,image_ps,pp_grad,ps_grad,illumination,maxiter,record_left_in_v);     //enter iteration

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	free1float(model.velp);
	free1float(model.vels);    
	free1float(model.sou);
	free2int(table);

    delete[] model.record_z;
    delete[] model.record_x; 

    delete[] model.gc;
    delete[] image_pp;
    delete[] image_ps;  
    delete[] image_pp_m;
    delete[] image_ps_m;          
    delete[] pp_grad;
    delete[] ps_grad;
	return 0;

}











