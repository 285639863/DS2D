
#include "laplace.h"
void laplace_filter_2d(int nz,  int nx, float *dat)
{
        float coea = 0.5, coeb = 0.25 ;
        int   n1 = 4, n2 = 1, n3 = 2;

        float *tmp1,*tmp2,*tmp0;

        tmp1 = new float[(nx+4)*(nz+4)];
        tmp2 = new float[(nx+4)*(nz+4)];
        tmp0 = new float[(nx+4)*(nz+4)];
     
        for(int ix =0; ix<nx; ix++)
           for(int iz =0; iz<nz; iz++){
              tmp1[(ix+2)*nz+(iz+2)] = dat[ix*nz+iz];
              tmp0[(ix+2)*nz+(iz+2)] = dat[ix*nz+iz];
        }

        for(int i=0; i<n1; i++){
           for(int ix=0; ix<nx; ix++)
               for(int iz=0; iz<nz; iz++){
                  tmp2[(ix+2)*nz+(iz+2)]=tmp1[(ix+2)*nz+(iz+2)]*coea+tmp1[(ix+2)*nz+(iz+1)]*coeb+tmp1[(ix+2)*nz+(iz+3)]*coeb;
           }
                     
           for(int ix=0; ix<nx; ix++)
               for(int iz=0; iz<nz; iz++){
                  tmp1[(ix+2)*nz+(iz+2)]=tmp2[(ix+2)*nz+(iz+2)]*coea+tmp2[(ix+1)*nz+(iz+2)]*coeb+tmp2[(ix+3)*nz+(iz+2)]*coeb;
           }
    
        }
        

        for(int ix=0; ix<nx; ix++)
               for(int iz=0; iz<nz; iz++){
                  tmp2[(ix+2)*nz+(iz+2)] = tmp0[(ix+2)*nz+(iz+2)] - tmp1[(ix+2)*nz+(iz+2)];
        }


        for(int i=0; i<n3; i++){
           for(int ix=0; ix<nx; ix++)
               for(int iz=0; iz<nz; iz++){
                  tmp1[(ix+2)*nz+(iz+2)]=tmp2[(ix+2)*nz+(iz+2)]*coea+tmp2[(ix+1)*nz+(iz+2)]*coeb+tmp2[(ix+3)*nz+(iz+2)]*coeb;
           }
         
        }

        for(int ix =0; ix<nx; ix++)
           for(int iz =0; iz<nz; iz++){
              dat[ix*nz+iz] = tmp1[(ix+2)*nz+(iz+2)];
        }

        
        delete[] tmp1;
        delete[] tmp2;
        delete[] tmp0;
}

void cutM(float *x, int nz, int nx,int zi, int ze, int xi, int xe)
{

        for (int ix=0; ix<xi; ++ix)
            for (int iz=0; iz<nz; ++iz)
            x[ix*nz+iz] = 0.f;

        for (int ix=xe; ix<nx; ++ix)
            for (int iz=0; iz<nz; ++iz)
            x[ix*nz+iz] = 0.f;

        for (int ix=0; ix<nx; ++ix)
            for (int iz=0; iz<zi; ++iz)
            x[ix*nz+iz] = 0.f;

        for (int ix=0; ix<nx; ++ix)
            for (int iz=ze; iz<nz; ++iz)
            x[ix*nz+iz] = 0.f;

}



void laplace_filter_2d(int nz,  int nx, int order, float dx, float *dat)
{
        float *tmp1;
        float *tmp2;
        tmp1 = new float[(nx+4)*(nz+4)]{};
        tmp2 = new float[nx*nz]{};

        for(int ix =0; ix<nx; ix++)
           for(int iz =0; iz<nz; iz++){
              tmp1[ix*nz+iz] = dat[ix*nz+iz];
        }


        for(int ix =2; ix<nx-2; ix++)
           for(int iz =2; iz<nz-2; iz++){
              tmp2[ix*nz+iz] = (1.0/dx/dx)*(tmp1[(ix+1)*nz+iz]+tmp1[(ix-1)*nz+iz]-2*tmp1[ix*nz+iz]+tmp1[ix*nz+iz+1]+tmp1[ix*nz+iz-1]-2*tmp1[ix*nz+iz]);
        }   

        for(int ix =0; ix<nx; ix++)
           for(int iz =0; iz<nz; iz++){
              dat[ix*nz+iz] = tmp2[ix*nz+iz];
        }

        delete[] tmp1;
        delete[] tmp2;
}