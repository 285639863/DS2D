#include"FD2DGPU.cuh"

__global__ void update_vx(float *txx,float *txz,float *vx,float *vx_x,float *vx_z,\
                const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
				const int nop,float *rho,float *dampx,float *dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
    const int index = ix*nzpml+iz;    
//	__syncthreads();
		float damp1 = 1 - dt*dampx[index]/2;
		float damp2 = 1 + dt*dampx[index]/2;
		float damp3 = 1 - dt*dampz[index]/2;
		float damp4 = 1 + dt*dampz[index]/2;

		float tmp_txx = 0;
		float tmp_txz = 0;
#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{

			 tmp_txx += coeff2[i]*(txx[(ix+i)*nzpml+iz]-txx[(ix-i+1)*nzpml+iz]);
			 tmp_txz += coeff2[i]*(txz[ix*nzpml+(iz+i-1)]-txz[ix*nzpml+(iz-i)]);			
		}
		
		vx_x[index] = (damp1*vx_x[index]+direction*(1.0/rho[index])*(dt/dx)*tmp_txx)/damp2;
		vx_z[index] = (damp3*vx_z[index]+direction*(1.0/rho[index])*(dt/dz)*tmp_txz)/damp4;
		vx[index] = vx_x[index] + vx_z[index];

}



__global__ void update_vz(float *tzz,float *txz,float *vz,float *vz_x,float *vz_z,\
                const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
				const int nop,float *rho,float *dampx,float *dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
    const int index = ix*nzpml+iz;    
//	__syncthreads();
		float damp1 = 1 - dt*dampx[index]/2;
		float damp2 = 1 + dt*dampx[index]/2;
		float damp3 = 1 - dt*dampz[index]/2;
		float damp4 = 1 + dt*dampz[index]/2;

		float tmp_tzz = 0;
		float tmp_txz = 0;
#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{

			 tmp_tzz += coeff2[i]*(tzz[ix*nzpml+(iz+i)]-tzz[ix*nzpml+(iz-i+1)]);
			 tmp_txz += coeff2[i]*(txz[(ix+i-1)*nzpml+iz]-txz[(ix-i)*nzpml+iz]);			
		}
		
		vz_x[index] = (damp1*vz_x[index]+direction*(1.0/rho[index])*(dt/dx)*tmp_txz)/damp2;		
		vz_z[index] = (damp3*vz_z[index]+direction*(1.0/rho[index])*(dt/dz)*tmp_tzz)/damp4;
		vz[index] = vz_x[index] + vz_z[index];

}


__global__ void update_txx_tzz(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *vx,float *vz,\
    const int nxpml, const int nzpml,const float dt,const float dx,const float dz,const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
    const int index = ix*nzpml+iz;
//	__syncthreads();

		float damp1 = 1 - dt*dampx[index]/2;
		float damp2 = 1 + dt*dampx[index]/2;
		float damp3 = 1 - dt*dampz[index]/2;
		float damp4 = 1 + dt*dampz[index]/2;

		float tmp_vx = 0;
		float tmp_vz = 0;

#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{

			 tmp_vx += coeff2[i]*(vx[(ix+i-1)*nzpml+iz]-vx[(ix-i)*nzpml+iz]);
			 tmp_vz += coeff2[i]*(vz[ix*nzpml+(iz+i-1)]-vz[ix*nzpml+(iz-i)]);			
		}
		 txx_x[index] = (damp1*txx_x[index]+direction*(lamda[index]+2*miu[index])*(dt/dx)*tmp_vx)/damp2;
		 txx_z[index] = (damp3*txx_z[index]+direction*lamda[index]*(dt/dz)*tmp_vz)/damp4;
		 tzz_x[index] = (damp1*tzz_x[index]+direction*lamda[index]*(dt/dx)*tmp_vx)/damp2;
		 tzz_z[index] = (damp3*tzz_z[index]+direction*(lamda[index]+2*miu[index])*(dt/dz)*tmp_vz)/damp4;
		 txx[index] = txx_x[index] + txx_z[index];
 		 tzz[index] = tzz_x[index] + tzz_z[index];

}


__global__ void update_txz(float *txz,float *txz_x,float *txz_z,float *vx,float *vz,const int nxpml, const int nzpml,\
	const float dt,const float dx,const float dz,const int nop,float *miu,float* dampx,float *dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
    const int index = ix*nzpml+iz;
//	__syncthreads();

		float damp1 = 1 - dt*dampx[index]/2;
		float damp2 = 1 + dt*dampx[index]/2;
		float damp3 = 1 - dt*dampz[index]/2;
		float damp4 = 1 + dt*dampz[index]/2;

		float tmp_vx = 0;
		float tmp_vz = 0;

#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{

			 tmp_vx += coeff2[i]*(vx[ix*nzpml+(iz+i)]-vx[ix*nzpml+(iz-i+1)]);
			 tmp_vz += coeff2[i]*(vz[(ix+i)*nzpml+iz]-vz[(ix-i+1)*nzpml+iz]);			
		}
		
        txz_x[index] = (damp1*txz_x[index]+direction*(dt/dx)*miu[index]*tmp_vz)/damp2;
        txz_z[index] = (damp3*txz_z[index]+direction*(dt/dz)*miu[index]*tmp_vx)/damp4;
        txz[index] = txz_x[index] + txz_z[index];

}


__global__ void update_vx_vz(float *txx,float *tzz,float *txz,float *vx,float *vx_x,float *vx_z,float *vz,float *vz_x,float *vz_z,\
                const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
				const int nop,float *rho,float *dampx,float *dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
    const int index = ix*nzpml+iz;    
//	__syncthreads();
		float damp1 = 1 - dt*dampx[index]/2;
		float damp2 = 1 + dt*dampx[index]/2;
		float damp3 = 1 - dt*dampz[index]/2;
		float damp4 = 1 + dt*dampz[index]/2;

		float tmp_txx_x = 0;
		float tmp_txz_z = 0;

		float tmp_tzz_z = 0;
		float tmp_txz_x = 0;		
#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			 tmp_txx_x += coeff2[i]*(txx[(ix+i)*nzpml+iz]-txx[(ix-i+1)*nzpml+iz]);
			 tmp_txz_z += coeff2[i]*(txz[ix*nzpml+(iz+i-1)]-txz[ix*nzpml+(iz-i)]);	

			 tmp_tzz_z += coeff2[i]*(tzz[ix*nzpml+(iz+i)]-tzz[ix*nzpml+(iz-i+1)]);
			 tmp_txz_x += coeff2[i]*(txz[(ix+i-1)*nzpml+iz]-txz[(ix-i)*nzpml+iz]);				 		
		}
		
		vx_x[index] = (damp1*vx_x[index]+direction*(1.0/rho[index])*(dt/dx)*tmp_txx_x)/damp2;
		vx_z[index] = (damp3*vx_z[index]+direction*(1.0/rho[index])*(dt/dz)*tmp_txz_z)/damp4;

		vz_x[index] = (damp1*vz_x[index]+direction*(1.0/rho[index])*(dt/dx)*tmp_txz_x)/damp2;		
		vz_z[index] = (damp3*vz_z[index]+direction*(1.0/rho[index])*(dt/dz)*tmp_tzz_z)/damp4;
	
		vx[index] = vx_x[index] + vx_z[index];
		vz[index] = vz_x[index] + vz_z[index];			

}




__global__ void update_txx_tzz_txz(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *txz,float *txz_x,float *txz_z,float *vx,float *vz,\
    const int nxpml, const int nzpml,const float dt,const float dx,const float dz,const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
    const int index = ix*nzpml+iz;
//	__syncthreads();

		float damp1 = 1 - dt*dampx[index]/2;
		float damp2 = 1 + dt*dampx[index]/2;
		float damp3 = 1 - dt*dampz[index]/2;
		float damp4 = 1 + dt*dampz[index]/2;

		float tmp_vx_x = 0;
		float tmp_vz_z = 0;
		float tmp_vx_z = 0;
		float tmp_vz_x = 0;


#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			 tmp_vx_x += coeff2[i]*(vx[(ix+i-1)*nzpml+iz]-vx[(ix-i)*nzpml+iz]);
			 tmp_vz_z += coeff2[i]*(vz[ix*nzpml+(iz+i-1)]-vz[ix*nzpml+(iz-i)]);	

			 tmp_vx_z += coeff2[i]*(vx[ix*nzpml+(iz+i)]-vx[ix*nzpml+(iz-i+1)]);
			 tmp_vz_x += coeff2[i]*(vz[(ix+i)*nzpml+iz]-vz[(ix-i+1)*nzpml+iz]);	

		}
		 txx_x[index] = (damp1*txx_x[index]+direction*(lamda[index]+2*miu[index])*(dt/dx)*tmp_vx_x)/damp2;
		 txx_z[index] = (damp3*txx_z[index]+direction*lamda[index]*(dt/dz)*tmp_vz_z)/damp4;
		 tzz_x[index] = (damp1*tzz_x[index]+direction*lamda[index]*(dt/dx)*tmp_vx_x)/damp2;
		 tzz_z[index] = (damp3*tzz_z[index]+direction*(lamda[index]+2*miu[index])*(dt/dz)*tmp_vz_z)/damp4;
         txz_x[index] = (damp1*txz_x[index]+direction*(dt/dx)*miu[index]*tmp_vz_x)/damp2;
         txz_z[index] = (damp3*txz_z[index]+direction*(dt/dz)*miu[index]*tmp_vx_z)/damp4;
	 
		 txx[index] = txx_x[index] + txx_z[index];
 		 tzz[index] = tzz_x[index] + tzz_z[index];
         txz[index] = txz_x[index] + txz_z[index];	
}





__global__ void update_txx_tzz_txz_grad(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *txz,float *txz_x,float *txz_z,float *vx,float *vz,\
    float *vx_gx,float *vz_gz,float *vx_gz,float *vz_gx,const int nxpml, const int nzpml,const float dt,const float dx,const float dz,const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction,int pml)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
    const int index = ix*nzpml+iz;
//	__syncthreads();

		float damp1 = 1 - dt*dampx[index]/2;
		float damp2 = 1 + dt*dampx[index]/2;
		float damp3 = 1 - dt*dampz[index]/2;
		float damp4 = 1 + dt*dampz[index]/2;

		float tmp_vx_x = 0;
		float tmp_vz_z = 0;
		float tmp_vx_z = 0;
		float tmp_vz_x = 0;


#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			 tmp_vx_x += coeff2[i]*(vx[(ix+i-1)*nzpml+iz]-vx[(ix-i)*nzpml+iz]);
			 tmp_vz_z += coeff2[i]*(vz[ix*nzpml+(iz+i-1)]-vz[ix*nzpml+(iz-i)]);	

			 tmp_vx_z += coeff2[i]*(vx[ix*nzpml+(iz+i)]-vx[ix*nzpml+(iz-i+1)]);
			 tmp_vz_x += coeff2[i]*(vz[(ix+i)*nzpml+iz]-vz[(ix-i+1)*nzpml+iz]);	

		}
		 txx_x[index] = (damp1*txx_x[index]+direction*(lamda[index]+2*miu[index])*(dt/dx)*tmp_vx_x)/damp2;
		 txx_z[index] = (damp3*txx_z[index]+direction*lamda[index]*(dt/dz)*tmp_vz_z)/damp4;
		 tzz_x[index] = (damp1*tzz_x[index]+direction*lamda[index]*(dt/dx)*tmp_vx_x)/damp2;
		 tzz_z[index] = (damp3*tzz_z[index]+direction*(lamda[index]+2*miu[index])*(dt/dz)*tmp_vz_z)/damp4;
         txz_x[index] = (damp1*txz_x[index]+direction*(dt/dx)*miu[index]*tmp_vz_x)/damp2;
         txz_z[index] = (damp3*txz_z[index]+direction*(dt/dz)*miu[index]*tmp_vx_z)/damp4;
	 
		 txx[index] = txx_x[index] + txx_z[index];
 		 tzz[index] = tzz_x[index] + tzz_z[index];
         txz[index] = txz_x[index] + txz_z[index];	

		 vx_gx[iz+ix*nzpml] = tmp_vx_x/dx;
		 vz_gz[iz+ix*nzpml] = tmp_vz_z/dz;

		 vx_gz[iz+ix*nzpml] = tmp_vx_z/dz;
		 vz_gx[iz+ix*nzpml] = tmp_vz_x/dx;

}



// adjoint wavefield
__global__ void GPUcalculate_elastic_vx_vz_back(float *txx,float *tzz,float *txz,float *vx,float *vx_x,float *vx_z,float *vz,float *vz_x,float *vz_z,float *lamda,float *miu,const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
									const int nop,float *rho,float *dampx,float *dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
//	__syncthreads();
		float damp1 = 1 - dt*dampx[iz+ix*nzpml]/2;
		float damp2 = 1 + dt*dampx[iz+ix*nzpml]/2;
		float damp3 = 1 - dt*dampz[iz+ix*nzpml]/2;
		float damp4 = 1 + dt*dampz[iz+ix*nzpml]/2;

		float tmp_txx_x = 0;
		float tmp_txz_z = 0;
		float tmp_tzz_x = 0;

		float tmp_tzz_z = 0;
		float tmp_txz_x = 0;
		float tmp_txx_z = 0;	


#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{		
			tmp_txx_x += coeff2[i]*((lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*txx[(ix+i)*nzpml+iz]-(lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*txx[(ix-i+1)*nzpml+iz]);
			tmp_txz_z += coeff2[i]*(miu[ix*nzpml+(iz+i-1)]*txz[ix*nzpml+(iz+i-1)]-miu[ix*nzpml+(iz-i)]*txz[ix*nzpml+(iz-i)]);				
			tmp_tzz_x += coeff2[i]*((lamda[(ix+i)*nzpml+iz])*tzz[(ix+i)*nzpml+iz]-(lamda[(ix-i+1)*nzpml+iz])*tzz[(ix-i+1)*nzpml+iz]);			
			
			tmp_tzz_z += coeff2[i]*((lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*tzz[ix*nzpml+(iz+i)]-(lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*tzz[ix*nzpml+(iz-i+1)]);
			tmp_txz_x += coeff2[i]*(miu[(ix+i-1)*nzpml+iz]*txz[(ix+i-1)*nzpml+iz]-miu[(ix-i)*nzpml+iz]*txz[(ix-i)*nzpml+iz]);				//TODO: miu and txz are not in a same grid point
			tmp_txx_z += coeff2[i]*((lamda[ix*nzpml+(iz+i)])*txx[ix*nzpml+(iz+i)]-(lamda[ix*nzpml+(iz-i+1)])*txx[ix*nzpml+(iz-i+1)]);					
		}
		
		vx_x[iz+ix*nzpml] = (damp1*vx_x[iz+ix*nzpml]+direction*(1.0/rho[iz+ix*nzpml])*(dt/dx)*tmp_txx_x \
							+direction*(1.0/rho[iz+ix*nzpml])*(dt/dx)*tmp_tzz_x)/damp2;
		vx_z[iz+ix*nzpml] = (damp3*vx_z[iz+ix*nzpml]+direction*(1.0/rho[iz+ix*nzpml])*(dt/dz)*tmp_txz_z)/damp4;					
		
		vz_x[iz+ix*nzpml] = (damp1*vz_x[iz+ix*nzpml]+direction*(1.0/rho[iz+ix*nzpml])*(dt/dx)*tmp_txz_x)/damp2;		
		vz_z[iz+ix*nzpml] = (damp3*vz_z[iz+ix*nzpml]+direction*(1.0/rho[iz+ix*nzpml])*(dt/dz)*tmp_tzz_z \
							+direction*(1.0/rho[iz+ix*nzpml])*(dt/dz)*tmp_txx_z)/damp4;
		
		vz[iz+ix*nzpml] = vz_x[iz+ix*nzpml] + vz_z[iz+ix*nzpml];		
		vx[iz+ix*nzpml] = vx_x[iz+ix*nzpml] + vx_z[iz+ix*nzpml];	
}


__global__ void GPUcalculate_elastic_txx_tzz_txz_back(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *txz,float *txz_x,float *txz_z,float *vx,float *vz,const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
	const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction)
{
	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
//	__syncthreads();

		float damp1 = 1 - dt*dampx[iz+ix*nzpml]/2;
		float damp2 = 1 + dt*dampx[iz+ix*nzpml]/2;
		float damp3 = 1 - dt*dampz[iz+ix*nzpml]/2;
		float damp4 = 1 + dt*dampz[iz+ix*nzpml]/2;

		float tmp_vx_x = 0;
		float tmp_vz_z = 0;

		float tmp_vx_z = 0;
		float tmp_vz_x = 0;

#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			tmp_vx_x += coeff2[i]*(vx[(ix+i-1)*nzpml+iz]-vx[(ix-i)*nzpml+iz]);
			tmp_vz_z += coeff2[i]*(vz[ix*nzpml+(iz+i-1)]-vz[ix*nzpml+(iz-i)]);

			tmp_vx_z += coeff2[i]*(vx[ix*nzpml+(iz+i)]-vx[ix*nzpml+(iz-i+1)]);
			tmp_vz_x += coeff2[i]*(vz[(ix+i)*nzpml+iz]-vz[(ix-i+1)*nzpml+iz]);			
		}
		
		 txx_x[iz+ix*nzpml] = (damp1*txx_x[iz+ix*nzpml]+direction*(dt/dx)*tmp_vx_x)/damp2;
		 tzz_z[iz+ix*nzpml] = (damp3*tzz_z[iz+ix*nzpml]+direction*(dt/dz)*tmp_vz_z)/damp4;
		 txz_x[iz+ix*nzpml] = (damp1*txz_x[iz+ix*nzpml]+direction*(dt/dx)*tmp_vz_x)/damp2;
		 txz_z[iz+ix*nzpml] = (damp3*txz_z[iz+ix*nzpml]+direction*(dt/dz)*tmp_vx_z)/damp4;

		 txx[iz+ix*nzpml] = txx_x[iz+ix*nzpml];
 		 tzz[iz+ix*nzpml] = tzz_z[iz+ix*nzpml];
		 txz[iz+ix*nzpml] = txz_x[iz+ix*nzpml] + txz_z[iz+ix*nzpml];	

}






__global__ void add_born_source_components(float *vx,float *vz,float *b_txx_x,float *b_txx_z, float *b_tzz_x,float *b_tzz_z,float *b_txz_x,float *b_txz_z,float *delta_mp,float *delta_ms,\
    const int nxpml, const int nzpml,const int pml,const float dx,const float dz,const float dt,const int nop,float *lamda,float *miu,float *rho,float *vp,float *vs,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	// if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
	if(iz>=nzpml-pml||ix>=nxpml-pml||iz<pml||ix<pml)return;
    const int index = ix*nzpml+iz;
//add virtual source 

		float tmp_vx_x = 0;
		float tmp_vz_z = 0;
		float tmp_vx_z = 0;
		float tmp_vz_x = 0;
#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			 tmp_vx_x += coeff2[i]*(vx[(ix+i-1)*nzpml+iz]-vx[(ix-i)*nzpml+iz]);
			 tmp_vz_z += coeff2[i]*(vz[ix*nzpml+(iz+i-1)]-vz[ix*nzpml+(iz-i)]);	

			 tmp_vx_z += coeff2[i]*(vx[ix*nzpml+(iz+i)]-vx[ix*nzpml+(iz-i+1)]);
			 tmp_vz_x += coeff2[i]*(vz[(ix+i)*nzpml+iz]-vz[(ix-i+1)*nzpml+iz]);				 		
		}

//use refl vp and vs
			b_txx_x[index] = b_txx_x[index] + direction*dt* ( 2*(lamda[index]+2*miu[index])*(delta_mp[index])*(tmp_vx_x)/dx );

			b_txx_z[index] = b_txx_z[index] + direction*dt* ( 2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*(delta_ms[index]))*(tmp_vz_z)/dz );

			b_tzz_x[index] = b_tzz_x[index] + direction*dt* ( 2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*(delta_ms[index]))*(tmp_vx_x)/dx );

			b_tzz_z[index] = b_tzz_z[index] + direction*dt* ( 2*(lamda[index]+2*miu[index])*(delta_mp[index])*(tmp_vz_z)/dz );

			b_txz_x[index] = b_txz_x[index] + direction*dt* ( 2*miu[index]*(delta_ms[index])*((tmp_vz_x)/dx ) );

			b_txz_z[index] = b_txz_z[index] + direction*dt* ( 2*miu[index]*(delta_ms[index])*((tmp_vx_z)/dz) );		

}



__global__ void add_born_source_adjoint_components(float *txx,float *tzz,float *txz,float *b_vx, float *b_vz,float *b_vx_x, float *b_vx_z,float *b_vz_x, float *b_vz_z,float *delta_mp,float *delta_ms,\
    const int nxpml, const int nzpml,const int pml,const float dx,const float dz,const float dt,const int nop,float *lamda,float *miu,float *rho,float *vp,float *vs,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	// if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
	if(iz>=nzpml-pml||ix>=nxpml-pml||iz<pml||ix<pml)return;
    const int index = ix*nzpml+iz;
//add virtual source 
		float scale =1.0;

		float tmp_txx_x = 0;
		float tmp_txx_z = 0;
		float tmp_tzz_x = 0;
		float tmp_tzz_z = 0;
		float tmp_txz_x = 0;
		float tmp_txz_z = 0;

#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{

// use refl vp and vs 

			 tmp_txx_x += coeff2[i]*(2*(lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*(delta_mp[(ix+i)*nzpml+iz])*txx[(ix+i)*nzpml+iz] - \
			 2*(lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*(delta_mp[(ix-i+1)*nzpml+iz])*txx[(ix-i+1)*nzpml+iz]);

			 tmp_txx_z += coeff2[i]*(2*((lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*(delta_mp[ix*nzpml+(iz+i)])-2*miu[ix*nzpml+(iz+i)]*(delta_ms[ix*nzpml+(iz+i)]))*txx[ix*nzpml+(iz+i)] - \
			 2*((lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*(delta_mp[ix*nzpml+(iz-i+1)])-2*miu[ix*nzpml+(iz-i+1)]*(delta_ms[ix*nzpml+(iz-i+1)]))*txx[ix*nzpml+(iz-i+1)]);	

			 tmp_tzz_x += coeff2[i]*(2*((lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*(delta_mp[(ix+i)*nzpml+iz])-2*miu[(ix+i)*nzpml+iz]*(delta_ms[(ix+i)*nzpml+iz]))*tzz[(ix+i)*nzpml+iz] - \
			 2*((lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*(delta_mp[(ix-i+1)*nzpml+iz])-2*miu[(ix-i+1)*nzpml+iz]*(delta_ms[(ix-i+1)*nzpml+iz]))*tzz[(ix-i+1)*nzpml+iz]);

			 tmp_tzz_z += coeff2[i]*(2*(lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*(delta_mp[ix*nzpml+(iz+i)])*tzz[ix*nzpml+(iz+i)] - \
			 2*(lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*(delta_mp[ix*nzpml+(iz-i+1)])*tzz[ix*nzpml+(iz-i+1)]);	

			 tmp_txz_x += coeff2[i]*(2*miu[(ix+i-1)*nzpml+iz]*(delta_ms[(ix+i-1)*nzpml+iz])*txz[(ix+i-1)*nzpml+iz] - \
			 2*miu[(ix+i-1)*nzpml+iz]*(delta_ms[(ix+i-1)*nzpml+iz])*txz[(ix-i)*nzpml+iz]);

			 tmp_txz_z += coeff2[i]*(2*miu[ix*nzpml+(iz+i-1)]*(delta_ms[ix*nzpml+(iz+i-1)])*txz[ix*nzpml+(iz+i-1)] - \
			 2*miu[ix*nzpml+(iz-i)]*(delta_ms[ix*nzpml+(iz-i)])*txz[ix*nzpml+(iz-i)]);		


		}

		b_vx_x[index] = b_vx_x[index] + direction*dt*((1.0/rho[index])*(tmp_txx_x/dx + tmp_tzz_x/dx ));
		b_vx_z[index] = b_vx_z[index] + direction*dt*((1.0/rho[index])*( tmp_txz_z/dz ));		

		b_vz_x[index] = b_vz_x[index] + direction*dt*((1.0/rho[index])*( tmp_txz_x/dx));
		b_vz_z[index] = b_vz_z[index] + direction*dt*((1.0/rho[index])*(tmp_tzz_z/dz + tmp_txx_z/dz ));		


}




