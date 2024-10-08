#include"FD2DGPU.cuh"
__global__ void born_vx(float *txx,float *txz,float *vx,float *vx_x,float *vx_z,float *b_txx,float *b_txz,float *b_vx,float *b_vx_x,float *b_vx_z,\
                const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
				const int nop,float *rho,float *dampx,float *dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>nzpml-nop||ix>nxpml-nop||iz<nop||ix<nop)return;
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
			tmp_txx += coeff2[i]*(txx[(ix+i-1)*nzpml+iz]-txx[(ix-i)*nzpml+iz]);
			tmp_txz += coeff2[i]*(txz[ix*nzpml+(iz+i-1)]-txz[ix*nzpml+(iz-i)]);
		}
		
		vx_x[index] = (damp1*vx_x[index]+direction*(1.0/rho[index])*(dt/dx)*tmp_txx)/damp2;
		vx_z[index] = (damp3*vx_z[index]+direction*(1.0/rho[index])*(dt/dz)*tmp_txz)/damp4;
		vx[index] = vx_x[index] + vx_z[index];

// born forward
		tmp_txx = 0;
		tmp_txz = 0;
#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			tmp_txx += coeff2[i]*(b_txx[(ix+i-1)*nzpml+iz]-b_txx[(ix-i)*nzpml+iz]);
			tmp_txz += coeff2[i]*(b_txz[ix*nzpml+(iz+i-1)]-b_txz[ix*nzpml+(iz-i)]);
		}
		
		b_vx_x[index] = (damp1*b_vx_x[index]+direction*(1.0/rho[index])*(dt/dx)*tmp_txx)/damp2;
		b_vx_z[index] = (damp3*b_vx_z[index]+direction*(1.0/rho[index])*(dt/dz)*tmp_txz)/damp4;
		b_vx[index] = b_vx_x[index] + b_vx_z[index];

//
}



__global__ void born_vz(float *tzz,float *txz,float *vz,float *vz_x,float *vz_z,float *b_tzz,float *b_txz,float *b_vz,float *b_vz_x,float *b_vz_z,\
                const int nxpml, const int nzpml,const float dt,const float dx,const float dz,\
				const int nop,float *rho,float *dampx,float *dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>nzpml-nop||ix>nxpml-nop||iz<nop||ix<nop)return;
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
			tmp_txz += coeff2[i]*(txz[(ix+i)*nzpml+iz]-txz[(ix-i+1)*nzpml+iz]);
		}
		
		vz_x[index] = (damp1*vz_x[index]+direction*(1.0/rho[index])*(dt/dx)*tmp_txz)/damp2;		
		vz_z[index] = (damp3*vz_z[index]+direction*(1.0/rho[index])*(dt/dz)*tmp_tzz)/damp4;
		vz[index] = vz_x[index] + vz_z[index];

//born forward
		tmp_tzz = 0;
		tmp_txz = 0;
#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			tmp_tzz += coeff2[i]*(b_tzz[ix*nzpml+(iz+i)]-b_tzz[ix*nzpml+(iz-i+1)]);
			tmp_txz += coeff2[i]*(b_txz[(ix+i)*nzpml+iz]-b_txz[(ix-i+1)*nzpml+iz]);
		}
		
		b_vz_x[index] = (damp1*b_vz_x[index]+direction*(1.0/rho[index])*(dt/dx)*tmp_txz)/damp2;		
		b_vz_z[index] = (damp3*b_vz_z[index]+direction*(1.0/rho[index])*(dt/dz)*tmp_tzz)/damp4;
		b_vz[index] = b_vz_x[index] + b_vz_z[index];

}


__global__ void born_txx_tzz(float *txx, float *tzz,float *txx_x,float *txx_z,float *tzz_x,float *tzz_z,float *vx,float *vz,\
    float *b_txx, float *b_tzz,float *b_txx_x,float *b_txx_z,float *b_tzz_x,float *b_tzz_z,float *b_vx,float *b_vz,float *delta_mp,float *delta_ms,float *vp,float *vs,\
    const int nxpml, const int nzpml,const float dt,const float dx,const float dz,const int nop,float *lamda,float *miu,float* dampx,float* dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>nzpml-nop||ix>nxpml-nop||iz<nop||ix<nop)return;
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
			tmp_vx += coeff2[i]*(vx[(ix+i)*nzpml+iz]-vx[(ix-i+1)*nzpml+iz]);
			tmp_vz += coeff2[i]*(vz[ix*nzpml+(iz+i-1)]-vz[ix*nzpml+(iz-i)]);
		}
		 txx_x[index] = (damp1*txx_x[index]+direction*(lamda[index]+2*miu[index])*(dt/dx)*tmp_vx)/damp2;
		 txx_z[index] = (damp3*txx_z[index]+direction*lamda[index]*(dt/dz)*tmp_vz)/damp4;
		 tzz_x[index] = (damp1*tzz_x[index]+direction*lamda[index]*(dt/dx)*tmp_vx)/damp2;
		 tzz_z[index] = (damp3*tzz_z[index]+direction*(lamda[index]+2*miu[index])*(dt/dz)*tmp_vz)/damp4;
		 txx[index] = txx_x[index] + txx_z[index];
 		 tzz[index] = tzz_x[index] + tzz_z[index];

//add virtual source 
        float scale = 1.0;
        // float scale = 1.0*dt; 
		// b_txx_x[index] += scale*2*(lamda[index]+2*miu[index])*(delta_mp[index]/vp[index])*tmp_vx/dx;							//TODO: divide maybe not necessary
		// b_txx_z[index] += scale*2*((lamda[index]+2*miu[index])*(delta_mp[index]/vp[index])-2*miu[index]*delta_ms[index]/vs[index])*tmp_vz/dz;
		// b_tzz_x[index] += scale*2*((lamda[index]+2*miu[index])*(delta_mp[index]/vp[index])-2*miu[index]*delta_ms[index]/vs[index])*tmp_vx/dx;
		// b_tzz_z[index] += scale*2*(lamda[index]+2*miu[index])*(delta_mp[index]/vp[index])*tmp_vz/dz;

		b_txx_x[index] += scale*2*(lamda[index]+2*miu[index])*(delta_mp[index])*tmp_vx/dx;							//TODO: divide maybe not necessary
		b_txx_z[index] += scale*2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*tmp_vz/dz;
		b_tzz_x[index] += scale*2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*tmp_vx/dx;
		b_tzz_z[index] += scale*2*(lamda[index]+2*miu[index])*(delta_mp[index])*tmp_vz/dz;		
//
		tmp_vx = 0;
		tmp_vz = 0;

#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			tmp_vx += coeff2[i]*(b_vx[(ix+i)*nzpml+iz]-b_vx[(ix-i+1)*nzpml+iz]);
			tmp_vz += coeff2[i]*(b_vz[ix*nzpml+(iz+i-1)]-b_vz[ix*nzpml+(iz-i)]);
		}
		 b_txx_x[index] = (damp1*b_txx_x[index]+direction*(lamda[index]+2*miu[index])*(dt/dx)*tmp_vx)/damp2;
		 b_txx_z[index] = (damp3*b_txx_z[index]+direction*lamda[index]*(dt/dz)*tmp_vz)/damp4;
		 b_tzz_x[index] = (damp1*b_tzz_x[index]+direction*lamda[index]*(dt/dx)*tmp_vx)/damp2;
		 b_tzz_z[index] = (damp3*b_tzz_z[index]+direction*(lamda[index]+2*miu[index])*(dt/dz)*tmp_vz)/damp4;
		 b_txx[index] = b_txx_x[index] + b_txx_z[index];
 		 b_tzz[index] = b_tzz_x[index] + b_tzz_z[index];


}


__global__ void born_txz(float *txz,float *txz_x,float *txz_z,float *vx,float *vz,float *b_txz,float *b_txz_x,float *b_txz_z,float *b_vx,float *b_vz,float *delta_ms, float *vs,\ 
    const int nxpml, const int nzpml,const float dt,const float dx,const float dz,const int nop,float *miu,const float* dampx,float *dampz,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	if(iz>nzpml-nop||ix>nxpml-nop||iz<nop||ix<nop)return;
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
			tmp_vz += coeff2[i]*(vz[(ix+i-1)*nzpml+iz]-vz[(ix-i)*nzpml+iz]);
		}
		
        txz_x[index] = (damp1*txz_x[index]+direction*(dt/dx)*miu[index]*tmp_vz)/damp2;
        txz_z[index] = (damp3*txz_z[index]+direction*(dt/dz)*miu[index]*tmp_vx)/damp4;
        txz[index] = txz_x[index] + txz_z[index];

//add  virtual source
        float scale = 1.0;
        // float scale = 1.0*dt;
		// b_txz_x[index] += scale*2*miu[index]*(delta_ms[index]/vs[index])*tmp_vz/dx;
		// b_txz_z[index] += scale*2*miu[index]*(delta_ms[index]/vs[index])*tmp_vx/dz;	

		b_txz_x[index] += scale*2*miu[index]*(delta_ms[index])*tmp_vz/dx;
		b_txz_z[index] += scale*2*miu[index]*(delta_ms[index])*tmp_vx/dz;	
//
		tmp_vx = 0;
		tmp_vz = 0;
#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			tmp_vx += coeff2[i]*(b_vx[ix*nzpml+(iz+i)]-b_vx[ix*nzpml+(iz-i+1)]);
			tmp_vz += coeff2[i]*(b_vz[(ix+i-1)*nzpml+iz]-b_vz[(ix-i)*nzpml+iz]);
		}
		
		  b_txz_x[index] = (damp1*b_txz_x[index]+direction*(dt/dx)*miu[index]*tmp_vz)/damp2;
		  b_txz_z[index] = (damp3*b_txz_z[index]+direction*(dt/dz)*miu[index]*tmp_vx)/damp4;
		  b_txz[index] = b_txz_x[index] + b_txz_z[index];

}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// rewrite born modeling schedue //////
/////////////////////////////////////////////////////////////////////////////////////////////////////
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
//			tmp_txx += coeff2[i]*(txx[(ix+i-1)*nzpml+iz]-txx[(ix-i)*nzpml+iz]);
//			tmp_txz += coeff2[i]*(txz[ix*nzpml+(iz+i-1)]-txz[ix*nzpml+(iz-i)]);
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
//			tmp_tzz += coeff2[i]*(tzz[ix*nzpml+(iz+i)]-tzz[ix*nzpml+(iz-i+1)]);
//			tmp_txz += coeff2[i]*(txz[(ix+i)*nzpml+iz]-txz[(ix-i+1)*nzpml+iz]);
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
//			tmp_vx += coeff2[i]*(vx[(ix+i)*nzpml+iz]-vx[(ix-i+1)*nzpml+iz]);
//			tmp_vz += coeff2[i]*(vz[ix*nzpml+(iz+i-1)]-vz[ix*nzpml+(iz-i)]);
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
//			tmp_vx += coeff2[i]*(vx[ix*nzpml+(iz+i)]-vx[ix*nzpml+(iz-i+1)]);
//			tmp_vz += coeff2[i]*(vz[(ix+i-1)*nzpml+iz]-vz[(ix-i)*nzpml+iz]);
			 tmp_vx += coeff2[i]*(vx[ix*nzpml+(iz+i)]-vx[ix*nzpml+(iz-i+1)]);
			 tmp_vz += coeff2[i]*(vz[(ix+i)*nzpml+iz]-vz[(ix-i+1)*nzpml+iz]);			
		}
		
        txz_x[index] = (damp1*txz_x[index]+direction*(dt/dx)*miu[index]*tmp_vz)/damp2;
        txz_z[index] = (damp3*txz_z[index]+direction*(dt/dz)*miu[index]*tmp_vx)/damp4;
        txz[index] = txz_x[index] + txz_z[index];

}

__global__ void add_born_source(float *vx,float *vz,float *b_txx, float *b_tzz,float *b_txz,float *delta_mp,float *delta_ms,\
    const int nxpml, const int nzpml,const int pml,const float dx,const float dz,const int nop,float *lamda,float *miu,float *rho,float *vp,float *vs,int direction)
{

	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
	// if(iz>=nzpml-nop||ix>=nxpml-nop||iz<nop||ix<nop)return;
	if(iz>=nzpml-pml||ix>=nxpml-pml||iz<pml||ix<pml)return;
    const int index = ix*nzpml+iz;
//add virtual source 
		float scale =1.0;

		// b_txx[index] += scale*2*(lamda[index]+2*miu[index])*(delta_mp[index])*(vx[(ix)*nzpml+iz]-vx[(ix-1)*nzpml+iz])/dx + \
		// 	2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*(vz[ix*nzpml+(iz)]-vz[ix*nzpml+(iz-1)])/dz;

		// b_tzz[index] += scale*2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*(vx[(ix)*nzpml+iz]-vx[(ix-1)*nzpml+iz])/dx + \
		// 	2*(lamda[index]+2*miu[index])*(delta_mp[index])*(vz[ix*nzpml+(iz)]-vz[ix*nzpml+(iz-1)])/dz;

		// b_txz[index] += scale*2*miu[index]*(delta_ms[index])*((vz[(ix+1)*nzpml+iz]-vz[ix*nzpml+iz])/dx + (vx[ix*nzpml+(iz+1)]-vx[ix*nzpml+iz])/dz);


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

// 		b_txx[index] = b_txx[index] + direction*scale*2*(lamda[index]+2*miu[index])*(delta_mp[index])*(tmp_vx_x)/dx + \
// 			2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*(tmp_vz_z)/dz;

// 		b_tzz[index] = b_tzz[index] + direction*scale*2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*(tmp_vx_x)/dx + \
// 			2*(lamda[index]+2*miu[index])*(delta_mp[index])*(tmp_vz_z)/dz;

// 		b_txz[index] = b_txz[index] + direction*scale*2*miu[index]*(delta_ms[index])*((tmp_vz_x)/dx + (tmp_vx_z)/dz);

		b_txx[index] = b_txx[index] + direction*scale*2*(lamda[index]+2*miu[index])*(-delta_mp[index])*(tmp_vx_x)/dx + \
			2*((lamda[index]+2*miu[index])*(-delta_mp[index])-2*miu[index]*(-delta_ms[index]))*(tmp_vz_z)/dz;

		b_tzz[index] = b_tzz[index] + direction*scale*2*((lamda[index]+2*miu[index])*(-delta_mp[index])-2*miu[index]*(-delta_ms[index]))*(tmp_vx_x)/dx + \
			2*(lamda[index]+2*miu[index])*(-delta_mp[index])*(tmp_vz_z)/dz;

		b_txz[index] = b_txz[index] + direction*scale*2*miu[index]*(-delta_ms[index])*((tmp_vz_x)/dx + (tmp_vx_z)/dz);




}


__global__ void add_born_source_adjoint(float *txx,float *tzz,float *txz,float *b_vx, float *b_vz,float *delta_mp,float *delta_ms,\
    const int nxpml, const int nzpml,const int pml,const float dx,const float dz,const int nop,float *lamda,float *miu,float *rho,float *vp,float *vs,int direction)
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

// TODO: lamda miu inner partial differential
		// tmp_txx_x = (txx[(ix+1)*nzpml+iz]-txx[(ix)*nzpml+iz])/(dx);
		// tmp_txx_z = (txx[ix*nzpml+(iz+1)]-txx[ix*nzpml+(iz)])/(dz);			
		// tmp_tzz_x = (tzz[(ix+1)*nzpml+iz]-tzz[(ix)*nzpml+iz])/(dx);
		// tmp_tzz_z = (tzz[ix*nzpml+(iz+1)]-tzz[ix*nzpml+(iz)])/(dz);	 

		// tmp_txz_x = (txz[(ix)*nzpml+iz] - txz[(ix-1)*nzpml+iz])/(dx);
		// tmp_txz_z = (txz[(ix)*nzpml+iz] - txz[ix*nzpml+(iz-1)])/(dz);

		// b_vx[index] = b_vx[index] + direction*(1.0/rho[index])*(2*(lamda[index]+2*miu[index])*(delta_mp[index])*(tmp_txx_x) + \
		// 	2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*(tmp_tzz_x) + \
		// 	2*miu[index]*delta_ms[index]*tmp_txz_z);

		// b_vz[index] = b_vz[index] + direction*(1.0/rho[index])*(2*(lamda[index]+2*miu[index])*(delta_mp[index])*(tmp_tzz_z) + \
		// 	2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*(tmp_txx_z) + \
		// 	2*miu[index]*delta_ms[index]*tmp_txz_x);




// #pragma unroll 4
// 		int i=1;
		
// 		tmp_txx_x += (2*(lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*(delta_mp[(ix+i)*nzpml+iz])*txx[(ix+i)*nzpml+iz] - \
// 		2*(lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*(delta_mp[(ix-i+1)*nzpml+iz])*txx[(ix-i+1)*nzpml+iz]);

// 		tmp_txx_z += (2*((lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*(delta_mp[ix*nzpml+(iz+i)])-2*miu[ix*nzpml+(iz+i)]*delta_ms[ix*nzpml+(iz+i)])*txx[ix*nzpml+(iz+i)] - \
// 		2*((lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*(delta_mp[ix*nzpml+(iz-i+1)])-2*miu[ix*nzpml+(iz-i+1)]*delta_ms[ix*nzpml+(iz-i+1)])*txx[ix*nzpml+(iz-i+1)]);	

// 		tmp_tzz_x += (2*((lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*(delta_mp[(ix+i)*nzpml+iz])-2*miu[(ix+i)*nzpml+iz]*delta_ms[(ix+i)*nzpml+iz])*tzz[(ix+i)*nzpml+iz] - \
// 		2*((lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*(delta_mp[(ix-i+1)*nzpml+iz])-2*miu[(ix-i+1)*nzpml+iz]*delta_ms[(ix-i+1)*nzpml+iz])*tzz[(ix-i+1)*nzpml+iz]);

// 		tmp_tzz_z += (2*(lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*(delta_mp[ix*nzpml+(iz+i)])*tzz[ix*nzpml+(iz+i)] - \
// 		2*(lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*(delta_mp[ix*nzpml+(iz-i+1)])*tzz[ix*nzpml+(iz-i+1)]);	

// 		tmp_txz_x += (2*miu[(ix+i-1)*nzpml+iz]*delta_ms[(ix+i-1)*nzpml+iz]*txz[(ix+i-1)*nzpml+iz] - \
// 		2*miu[(ix+i-1)*nzpml+iz]*delta_ms[(ix+i-1)*nzpml+iz]*txz[(ix-i)*nzpml+iz]);

// 		tmp_txz_z += (2*miu[ix*nzpml+(iz+i-1)]*delta_ms[ix*nzpml+(iz+i-1)]*txz[ix*nzpml+(iz+i-1)] - \
// 		2*miu[ix*nzpml+(iz-i)]*delta_ms[ix*nzpml+(iz-i)]*txz[ix*nzpml+(iz-i)]);			 		
		
// 		b_vx[index] = b_vx[index] + direction*(1.0/rho[index])*(tmp_txx_x/dx + tmp_tzz_x/dx + tmp_txz_z/dz);

// 		b_vz[index] = b_vz[index] + direction*(1.0/rho[index])*(tmp_tzz_z/dz + tmp_txx_z/dz + tmp_txz_x/dx);


// #pragma unroll 4
// 		for(int i=1;i<=nop;i++)
// 		{
// 			 tmp_txx_x += coeff2[i]*(2*(lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*(delta_mp[(ix+i)*nzpml+iz])*txx[(ix+i)*nzpml+iz] - \
// 			 2*(lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*(delta_mp[(ix-i+1)*nzpml+iz])*txx[(ix-i+1)*nzpml+iz]);

// 			 tmp_txx_z += coeff2[i]*(2*((lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*(delta_mp[ix*nzpml+(iz+i)])-2*miu[ix*nzpml+(iz+i)]*delta_ms[ix*nzpml+(iz+i)])*txx[ix*nzpml+(iz+i)] - \
// 			 2*((lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*(delta_mp[ix*nzpml+(iz-i+1)])-2*miu[ix*nzpml+(iz-i+1)]*delta_ms[ix*nzpml+(iz-i+1)])*txx[ix*nzpml+(iz-i+1)]);	

// 			 tmp_tzz_x += coeff2[i]*(2*((lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*(delta_mp[(ix+i)*nzpml+iz])-2*miu[(ix+i)*nzpml+iz]*delta_ms[(ix+i)*nzpml+iz])*tzz[(ix+i)*nzpml+iz] - \
// 			 2*((lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*(delta_mp[(ix-i+1)*nzpml+iz])-2*miu[(ix-i+1)*nzpml+iz]*delta_ms[(ix-i+1)*nzpml+iz])*tzz[(ix-i+1)*nzpml+iz]);

// 			 tmp_tzz_z += coeff2[i]*(2*(lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*(delta_mp[ix*nzpml+(iz+i)])*tzz[ix*nzpml+(iz+i)] - \
// 			 2*(lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*(delta_mp[ix*nzpml+(iz-i+1)])*tzz[ix*nzpml+(iz-i+1)]);	

// 			 tmp_txz_x += coeff2[i]*(2*miu[(ix+i-1)*nzpml+iz]*delta_ms[(ix+i-1)*nzpml+iz]*txz[(ix+i-1)*nzpml+iz] - \
// 			 2*miu[(ix+i-1)*nzpml+iz]*delta_ms[(ix+i-1)*nzpml+iz]*txz[(ix-i)*nzpml+iz]);

// 			 tmp_txz_z += coeff2[i]*(2*miu[ix*nzpml+(iz+i-1)]*delta_ms[ix*nzpml+(iz+i-1)]*txz[ix*nzpml+(iz+i-1)] - \
// 			 2*miu[ix*nzpml+(iz-i)]*delta_ms[ix*nzpml+(iz-i)]*txz[ix*nzpml+(iz-i)]);			 		
// 		}
// 		b_vx[index] = b_vx[index] + direction*(1.0/rho[index])*(tmp_txx_x/dx + tmp_tzz_x/dx + tmp_txz_z/dz);

// 		b_vz[index] = b_vz[index] + direction*(1.0/rho[index])*(tmp_tzz_z/dz + tmp_txx_z/dz + tmp_txz_x/dx);


#pragma unroll 4
		for(int i=1;i<=nop;i++)
		{
			 tmp_txx_x += coeff2[i]*(2*(lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*(-delta_mp[(ix+i)*nzpml+iz])*txx[(ix+i)*nzpml+iz] - \
			 2*(lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*(-delta_mp[(ix-i+1)*nzpml+iz])*txx[(ix-i+1)*nzpml+iz]);

			 tmp_txx_z += coeff2[i]*(2*((lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*(-delta_mp[ix*nzpml+(iz+i)])-2*miu[ix*nzpml+(iz+i)]*(-delta_ms[ix*nzpml+(iz+i)]))*txx[ix*nzpml+(iz+i)] - \
			 2*((lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*(-delta_mp[ix*nzpml+(iz-i+1)])-2*miu[ix*nzpml+(iz-i+1)]*(-delta_ms[ix*nzpml+(iz-i+1)]))*txx[ix*nzpml+(iz-i+1)]);	

			 tmp_tzz_x += coeff2[i]*(2*((lamda[(ix+i)*nzpml+iz]+2*miu[(ix+i)*nzpml+iz])*(-delta_mp[(ix+i)*nzpml+iz])-2*miu[(ix+i)*nzpml+iz]*(-delta_ms[(ix+i)*nzpml+iz]))*tzz[(ix+i)*nzpml+iz] - \
			 2*((lamda[(ix-i+1)*nzpml+iz]+2*miu[(ix-i+1)*nzpml+iz])*(-delta_mp[(ix-i+1)*nzpml+iz])-2*miu[(ix-i+1)*nzpml+iz]*(-delta_ms[(ix-i+1)*nzpml+iz]))*tzz[(ix-i+1)*nzpml+iz]);

			 tmp_tzz_z += coeff2[i]*(2*(lamda[ix*nzpml+(iz+i)]+2*miu[ix*nzpml+(iz+i)])*(-delta_mp[ix*nzpml+(iz+i)])*tzz[ix*nzpml+(iz+i)] - \
			 2*(lamda[ix*nzpml+(iz-i+1)]+2*miu[ix*nzpml+(iz-i+1)])*(-delta_mp[ix*nzpml+(iz-i+1)])*tzz[ix*nzpml+(iz-i+1)]);	

			 tmp_txz_x += coeff2[i]*(2*miu[(ix+i-1)*nzpml+iz]*(-delta_ms[(ix+i-1)*nzpml+iz])*txz[(ix+i-1)*nzpml+iz] - \
			 2*miu[(ix+i-1)*nzpml+iz]*(-delta_ms[(ix+i-1)*nzpml+iz])*txz[(ix-i)*nzpml+iz]);

			 tmp_txz_z += coeff2[i]*(2*miu[ix*nzpml+(iz+i-1)]*(-delta_ms[ix*nzpml+(iz+i-1)])*txz[ix*nzpml+(iz+i-1)] - \
			 2*miu[ix*nzpml+(iz-i)]*(-delta_ms[ix*nzpml+(iz-i)])*txz[ix*nzpml+(iz-i)]);			 		
		}
		b_vx[index] = b_vx[index] + direction*(1.0/rho[index])*(tmp_txx_x/dx + tmp_tzz_x/dx + tmp_txz_z/dz);

		b_vz[index] = b_vz[index] + direction*(1.0/rho[index])*(tmp_tzz_z/dz + tmp_txx_z/dz + tmp_txz_x/dx);




}


// __global__ void add_born_source(float *vx,float *vz,float *b_txx, float *b_tzz,float *b_txz,float *delta_mp,float *delta_ms,\
//     const int nxpml, const int nzpml,const float dx,const float dz,const int nop,float *lamda,float *miu,float *rho,float *vp,float *vs)
// {

// 	const int iz = blockIdx.x*blockDim.x+threadIdx.x;
// 	const int ix = blockIdx.y*blockDim.y+threadIdx.y;
// 	if(iz>nzpml-nop||ix>nxpml-nop||iz<nop||ix<nop)return;
//     const int index = ix*nzpml+iz;

// //add virtual source 
// 		float scale =1.0;

// 		float tmp_vx_z = 0;
// 		float tmp_vz_x = 0;
// 		float tmp_vx_x = 0;
// 		float tmp_vz_z = 0;

// #pragma unroll 4
// 		for(int i=1;i<=nop;i++)
// 		{
// 			tmp_vx_x += coeff2[i]*(vx[(ix+i)*nzpml+iz]-vx[(ix-i+1)*nzpml+iz]);
// 			tmp_vz_z += coeff2[i]*(vz[ix*nzpml+(iz+i-1)]-vz[ix*nzpml+(iz-i)]);			
// 			tmp_vx_z += coeff2[i]*(vx[ix*nzpml+(iz+i)]-vx[ix*nzpml+(iz-i+1)]);
// 			tmp_vz_x += coeff2[i]*(vz[(ix+i-1)*nzpml+iz]-vz[(ix-i)*nzpml+iz]);
// 		}


// 		b_txx[index] += scale*2*(lamda[index]+2*miu[index])*(delta_mp[index])*(tmp_vx_x)/dx + \
// 			2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*(tmp_vz_z)/dz;

// 		b_tzz[index] += scale*2*((lamda[index]+2*miu[index])*(delta_mp[index])-2*miu[index]*delta_ms[index])*(tmp_vx_x)/dx + \
// 			2*(lamda[index]+2*miu[index])*(delta_mp[index])*(tmp_vz_z)/dz;

// 		b_txz[index] += scale*2*miu[index]*(delta_ms[index])*((tmp_vz_x)/dx + (tmp_vx_z)/dz);

// }







