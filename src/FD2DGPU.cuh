#include "cuda_runtime.h"
#include "cublas_v2.h"
#include <iostream>
#include "curand_kernel.h"
#include<cufft.h>
using namespace std;
#define PI 3.14159265358979
#define BLOCKDIMX 32
#define BLOCKDIMY 32

#define CUFFT_CALL(call) {cufftResult result = call; if (result != CUFFT_SUCCESS) { \
    std::cerr << "CUFFT error: " << result << " at line " << __LINE__ << std::endl; \
    exit(EXIT_FAILURE); }}


__constant__ float coeff2[6]={0,1225.0/1024.0,-245.0/3072.0,49.0/5120.0,-5.0/7168.0,0};


    void normalize(float* matrix, int length, int tracenum);


__global__ void update_vx(float *txx,float *txz,float *vx,float *vx_x,float *vx_z,\
                const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
				const int nop,float *rho,float *dampx,float *dampz,int direction);
__global__ void update_vz(float *tzz,float *txz,float *vz,float *vz_x,float *vz_z,\
                const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
				const int nop,float *rho,float *dampx,float *dampz,int direction);
__global__ void update_txx_tzz(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *vx,float *vz,\
    const int nxpml, const int nzpml,const float dt,const float dx,const float dz,const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction);
__global__ void update_txz(float *txz,float *txz_x,float *txz_z,float *vx,float *vz,const int nxpml, const int nzpml,\
	const float dt,const float dx,const float dz,const int nop,float *miu,float* dampx,float *dampz,int direction);

__global__ void update_vx_vz(float *txx,float *tzz,float *txz,float *vx,float *vx_x,float *vx_z,float *vz,float *vz_x,float *vz_z,\
                const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
				const int nop,float *rho,float *dampx,float *dampz,int direction);
__global__ void update_txx_tzz_txz(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *txz,float *txz_x,float *txz_z,float *vx,float *vz,\
    const int nxpml, const int nzpml,const float dt,const float dx,const float dz,const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction);
__global__ void update_txx_tzz_txz_grad(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *txz,float *txz_x,float *txz_z,float *vx,float *vz,\
    float *vx_gx,float *vz_gz,float *vx_gz,float *vz_gx,const int nxpml, const int nzpml,const float dt,const float dx,const float dz,const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction,int pml);

__global__ void add_born_source_components(float *vx,float *vz,float *b_txx_x,float *b_txx_z, float *b_tzz_x,float *b_tzz_z,float *b_txz_x,float *b_txz_z,float *delta_mp,float *delta_ms,\
    const int nxpml, const int nzpml,const int pml,const float dx,const float dz,const float dt,const int nop,float *lamda,float *miu,float *rho,float *vp,float *vs,int direction);
__global__ void add_born_source_adjoint_components(float *txx,float *tzz,float *txz,float *b_vx, float *b_vz,float *b_vx_x, float *b_vx_z,float *b_vz_x, float *b_vz_z,float *delta_mp,float *delta_ms,\
    const int nxpml, const int nzpml,const int pml,const float dx,const float dz,const float dt,const int nop,float *lamda,float *miu,float *rho,float *vp,float *vs,int direction);    

__global__ void GPUcalculate_elastic_vx_vz_back(float *txx,float *tzz,float *txz,float *vx,float *vx_x,float *vx_z,float *vz,float *vz_x,float *vz_z,float *lamda,float *miu,const int nxpml, const int nzpml,\
    const float dt,const float dx,const float dz,const int nop,float *rho,float *dampx,float *dampz,int direction);
__global__ void GPUcalculate_elastic_txx_tzz_txz_back(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *txz,float *txz_x,float *txz_z,float *vx,float *vz,const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
	const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction);

struct cudasize{
	dim3 grid;
	dim3 block;
};

class FD2DGPU{

public:
	const float dx;
    const float dz;
    const float dt;
    const int nxpml;
    const int nzpml;
    const int nop;
    const int scale;
    const int pml;
    const int nt;
	int nx;
    int nz;
    const int allnx;
    const int allnz;

	struct cudasize Grid_Block;
	struct cudasize record;
    // struct cudasize velcopy;


    FD2DGPU(const float* sou,const float dx1, const float dz1, const float dt1, const int nxpml1, const int nzpml1, \
    const int allnx1,const int allnz1,const int scale1,const int pml1,const int nt1,const int nop1);

	~FD2DGPU();

    void GPUbufferV(float *v);
    void bufferVHtoD(const int ngx_min);
    void calculate_min(int *min_idx);
    void calculate_damp_k_rho();

//    void calculateP();
//    void calculateV();
//    void addS(const int sxnum,const float it);
    void recordD(const int scale,const int pml,const int it);
    float *vel;

//	float *px;
//	float *pz;
//  float *p;
    float *vx;
    float *vz;

    float *source;               //gpu

    float *d_k;             //gpu
    float *d_rho;           //gpu
    float *d_damp;          //gpu

	float *receiver;
    void testcopyV(float *h_v);
    void FD_initial(const int nDim);


    void addS(const int sxnum,const int it); 
    void calculateP();
    void calculateV();


private:
	float *px;
	float *pz;
    float *p; 
    float *vt;        
};


class FD2DGPU_ELASTIC:public FD2DGPU{

public:
//    FD2DGPU_ELASTIC(const float* sou,const int dx1, const int dz1, const float dt1, const int nxpml1, const int nzpml1, \
    const int allnx1,const int allnz1,const int scale1,const int pml1,const int nt1,const int nop1):
//    FD2DGPU(sou,dx1,dz1,dt1,nxpml1,nzpml1,allnx1,allnz1,scale1,pml1,nt1,nop1){};
    FD2DGPU_ELASTIC(const float* sou,const float dx1, const float dz1, const float dt1, const int nxpml1, const int nzpml1, \
    const int allnx1,const int allnz1,const int scale1,const int pml1,const int nt1,const int nop1,int ntr1);
    ~FD2DGPU_ELASTIC();
//  forward related
    int isx;
    int isz;

    int igz;            // the depth of receiver. constant number temporaly

    void elastic_rotation_addS_for(const int it,cudaStream_t stream=0);  
//
//  record related 
    float cmin;     //  left edge of calculation area
	struct cudasize adds_bac;
    float *GPUrecord;
    float *GPUrecord_x;    
    int *GPUgc;
    int ntr;
    void elastic_rotation_addS_bac(const int it,cudaStream_t stream=0,bool stack=false);
    void record_copytoGPU(float* record,float* record_x,cudaStream_t stream=0);
//
// backward variables & image 
    float* theta_bk;
    float* omega_bk;
    float* theta_x_bk;
    float* theta_z_bk;
    float* omega_x_bk;
    float* omega_z_bk;
    float* vpx_bk;
    float* vpz_bk;
    float* vsx_bk;
    float* vsz_bk;
    float* vx_bk;
    float* vz_bk;

    float *vp_bk;                              
    float *vs_bk;

// 
	float *vx_x_bk;
	float *vx_z_bk;
	float *vz_x_bk;
	float *vz_z_bk;
	float *txx_bk;
	float *txz_bk;
	float *txx_x_bk;
	float *txx_z_bk;
	float *txz_x_bk;
	float *txz_z_bk;

	float *tzz_x_bk;
	float *tzz_z_bk;
    float *tzz_bk;

    float *vx_gx;
    float *vx_gz;
    float *vz_gx;
    float *vz_gz;
//    

//
	float *b_txx;
	float *b_txx_x;
	float *b_txx_z;
	float *b_txz;
	float *b_txz_x;
	float *b_txz_z;

    float *b_tzz;
    float *b_tzz_x;
    float *b_tzz_z;

	float *b_vx_x;
	float *b_vx_z;
	float *b_vz_x;
	float *b_vz_z;


//
    void BK_ELASTIC_initial(const int nDim);    
    void image(cudaStream_t stream=0);
    struct cudasize imagesize;
//

// lsrtm
    void GPUbufferM(float *pp,float *ps);
    void buffer_image_extrapolation(const int ngx_min);    
    void rot_exp_calculate_v_born(int direction,bool forward,bool reconstruct,cudaStream_t stream=0,bool rbc=false);
    void rot_exp_calculate_p_born(int direction,bool forward,cudaStream_t stream=0,bool rbc=false);
    void record_scatter(const int it,bool double_scatter=false);   
    void subtract(double alpha,bool doublescatter=false);
    void scatter_add(double alpha);
    void image_gradient(cudaStream_t stream=0);    
    void imagebuffer_resettozero(float *cpu_pp_grad,float *cpu_ps_grad,float *illuminaiton,float *b_cpu_pp_grad,float *b_cpu_ps_grad,float *b_cpu_pp_grad2,float *b_cpu_ps_grad2);
    void LSRTM_initial(const int nDim);        
    float *s_image_pp_m;            // range of single shot image 
    float *s_image_ps_m;    
    float *image_pp_m;              // range of all image
    float *image_ps_m;
    float *image_pp;                // output image
    float *image_ps;
    float *b_theta;
    float *b_omega;
    float *b_theta_x;
    float *b_theta_z;
    float *b_omega_x;
    float *b_omega_z;
    float *b_vpx;
    float *b_vpz;
    float *b_vsx;
    float *b_vsz;
    float *b_vx;
    float *b_vz;
    float *vpx_g;
    float *vpz_g;
    float *vsx_g;
    float *vsz_g;
    float *GPUrecord_scatter_z;
    float *GPUrecord_scatter_x;
    float *numerator_pp;
    float *numerator_ps;
    float *denominator_pp;
	float *allimageGPU_numerator_pp;
	float *allimageGPU_numerator_ps;
	float *allimageGPU_denominator;

	float *single_pp_grad;
	float *single_ps_grad;
	float *pp_gradient;
	float *ps_gradient;

//  lsrtm
    void FD_ELASTIC_initial(const int nDim);
    void calculate_vp_min(int *min_idx);
    void calculate_vs_min(int *min_idx);    
    void calculate_elastic_damp_C();    
    void calculateTxxzz(int direction,bool forward,bool born=false);
    void calculateTxz(int direction,bool forward,bool born=false);
    void calculateVx(int direction,bool forward,bool born=false);
    void calculateVz(int direction,bool forward,bool born=false);
    void elastic_addS(const int sxnum,const int it);      
    void GPUbufferVPVS(float *vp_all,float *vs_all);
    void bufferVpVsHtoD(const int ngx_min); 
    void Helmholtz_VP_VS();
    void Vector_VP_VS();   
    void record_elastic(const int scale,const int pml,const int it,const int gx);
    void record_elastic_VP_VS(const int scale,const int pml,const int it);
    void snapcopy(struct snap snapVPVS,const int nDim);
    void elastic_addS_backward(const int it);

    float *vx_x;                    
    float *vx_z;
    float *vz_x;
    float *vz_z;   

    float *txx;
    float *tzz;
    float *txz;      

    float *txx_x;
    float *txx_z;
    float *tzz_x;
    float *tzz_z;
    float *txz_x;
    float *txz_z;

    float *VP;
    float *VS;

    float *receiver2;
    float *d_dampz;          //gpu


#ifdef VECTORIMG    
    float2* vectorVP;
    float2* vectorVS;

	float2* wavefield;
	float2* cpu_wavefield;
	float2* GPU_FW_wavefield;
#ifdef POYNTING
    float2* wavefield2;
	float2* GPU_FW_wavefield2;    
    float2* vp_source_up;
    float2* vp_source_down;
    float2* vp_receiver_up;
    float2* vp_receiver_down;
#endif
#else
	float* wavefield;
	float* cpu_wavefield;
	float* GPU_FW_wavefield;
#endif


	float *single_image;

	float *single_image_ps;

    float *vp;                              
    float *vs;

// ROT-EXP-EQUATION
    float* theta;
    float* omega;
    float* theta_x;
    float* theta_z;
    float* omega_x;
    float* omega_z;
    float* vpx;
    float* vpz;
    float* vsx;
    float* vsz;
    void elastic_rotation_addS(const int sxnum,const int it);
    void rot_exp_calculate_v(int direction,bool forward,cudaStream_t stream=0,bool rbc=false);
    void rot_exp_calculate_p(int direction,bool forward,cudaStream_t stream=0,bool rbc=false);
// ROT-EXP-EQUATION

// random boundary
void calculate_max(int *max_idx,float *v);
void set_random_boundary();
//

    void normal_forward_and_born(int direction);
//
    void add_virtual_source();
//
    void born_modeling(int direction,int forward_or_born);

    void record_copytoGPU(float* record,float* record_x,int* gc,cudaStream_t stream=0);

    int delta_left;
    int record_left_in_v;
    int ngx_left;
    void pure_born_reconstruct_backward(int direction,int forward);

    float *b_vx_bk;
    float *b_vx_x_bk;
    float *b_vx_z_bk;
    float *b_vz_bk;
    float *b_vz_x_bk;
    float *b_vz_z_bk;
    float *b_txx_bk;
    float *b_txx_x_bk;
    float *b_txx_z_bk;
    float *b_tzz_bk;
    float *b_tzz_x_bk;
    float *b_tzz_z_bk;
    float *b_txz_bk;
    float *b_txz_x_bk;
    float *b_txz_z_bk;

    float *b_vx_gx;
    float *b_vx_gz;    
    float *b_vz_gz;
    float *b_vz_gx;       

    float *b_numerator_pp;
    float *b_numerator_ps;
    float *b_denominator_pp;
	float *b_allimageGPU_numerator_pp;
	float *b_allimageGPU_numerator_ps;
	float *b_allimageGPU_denominator;

	float *b2_vx;
	float *b2_vx_x;
	float *b2_vx_z;
	float *b2_vz;
	float *b2_vz_x;
	float *b2_vz_z;
	float *b2_txx;
	float *b2_txx_x;
	float *b2_txx_z;
	float *b2_tzz;
	float *b2_tzz_x;
	float *b2_tzz_z;
	float *b2_txz;
	float *b2_txz_x;
	float *b2_txz_z;
	float *GPUrecord_scatter2_z;            
	float *GPUrecord_scatter2_x;   

private:

    float *lamda;
    float *miu;
    float *vel_s;


};


void ELSRTM_ROT_EXP_IMAGE_SINGLESHOT(FD2DGPU_ELASTIC &fdtdgpu,const int sxnum,const int nDim,const int myid,float *record,float *record_x,int *gc,\
    bool rbc,bool cpu_mem,int ntt,int iter,char* snapdir,float *image_pp,float *image_ps,float *cpu_pp_grad_old,float *cpu_ps_grad_old,double &misfit3);
