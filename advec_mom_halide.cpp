#include <Halide.h>
#include <cstdio>
#include "f_post_vol.h"
#include "f_pre_vol.h"
#include "ftocmacros.h"

using namespace Halide;

void advec_mom_kernel_halide(int *xmin,int *xmax,int *ymin,int *ymax,
                      double *xvel1,
                      double *yvel1,
                      double *mass_flux_x,
                      double *vol_flux_x,
                      double *mass_flux_y,
                      double *vol_flux_y,
                      double *volume,
                      double *density1,
                      double *node_flux,
                      double *node_mass_post,
                      double *node_mass_pre,
                      double *advec_vel,
                      double *mom_flux,
                      double *pre_vol,
                      double *post_vol,
                      double *celldx,
                      double *celldy,
                         int *whch_vl,
                         int *swp_nmbr,
                         int *drctn) {

  Buffer b_volume(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)volume);
  Buffer b_vol_flux_y(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)vol_flux_y);
  Buffer b_post_vol(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)post_vol);
  Buffer b_pre_vol(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)pre_vol);
  Buffer b_vol_flux_x(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)vol_flux_x);
//  Buffer b_mass_flux_x(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)mass_flux_x);

  f_post_vol(b_volume.raw_buffer(), b_vol_flux_y.raw_buffer(), b_post_vol.raw_buffer()); 
  f_pre_vol(b_vol_flux_x.raw_buffer(), b_vol_flux_y.raw_buffer(), b_volume.raw_buffer(), b_pre_vol.raw_buffer());

  printf("ok\n");

  int x_min=*xmin;
  int x_max=*xmax;
  int y_min=*ymin;
  int y_max=*ymax;
  int which_vel=*whch_vl;
  int sweep_number=*swp_nmbr;
  int direction=*drctn;

  int j,k,mom_sweep;
  int upwind,donor,downwind,dif;
  double sigma,wind,width;
  double vdiffuw,vdiffdw,auw,adw,limiter;

  double *vel1;

  if(which_vel==1){
    vel1=xvel1;
  } else{
    vel1=yvel1;
  }

  mom_sweep=direction+2*(sweep_number-1); 

#pragma omp parallel
{
   if(direction==1) {
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min-2;j<=x_max+2;j++) {
        node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                            *(mass_flux_x[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j+1,k-1,x_max+5,x_min-2,y_min-2)]
                                                             +mass_flux_x[FTNREF2D(j+1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min-1;j<=x_max+2;j++) {
        node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                                 *(density1[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min-1;j<=x_max+2;j++) {
        node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                          -node_flux[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]+node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma omp for private(upwind,downwind,donor,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,j)
    for (k=y_min;k<=y_max+1;k++) {
      for (j=x_min-1;j<=x_max+1;j++) {
        if(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]<0.0){
          upwind=j+2;
          donor=j+1;
          downwind=j;
          dif=donor;
        }
        else{
          upwind=j-1;
          donor=j;
          downwind=j+1;
          dif=upwind;
        }
        sigma=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]);
        width=celldx[FTNREF1D(j,x_min-2)];
        vdiffuw=vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(upwind,k  ,x_max+5,x_min-2,y_min-2)];
        vdiffdw=vel1[FTNREF2D(downwind,k  ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)];
        limiter=0.0;
        if(vdiffuw*vdiffdw>0.0){
          auw=fabs(vdiffuw);
          adw=fabs(vdiffdw);
          wind=1.0;
          if(vdiffdw<=0.0) wind=-1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldx[FTNREF1D(dif,x_min-2)])/6.0,MIN(auw,adw));
        }
        advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=vel1[FTNREF2D(donor,k  ,x_max+5,x_min-2,y_min-2)]+(1.0-sigma)*limiter;
        mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                           *node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }

#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=(vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        *node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        +mom_flux[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]
                                                        -mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                                        /node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
  }
  else if(direction==2){
#pragma omp for private(j)
    for (k=y_min-2;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                            *(mass_flux_y[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j-1,k+1,x_max+4,x_min-2,y_min-2)]
                                                             +mass_flux_y[FTNREF2D(j  ,k+1,x_max+4,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min-1;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=0.25
                                                                 *(density1[FTNREF2D(j  ,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j  ,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k-1,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k-1,x_max+5,x_min-2,y_min-2)]
                                                                  +density1[FTNREF2D(j-1,k  ,x_max+4,x_min-2,y_min-2)]
                                                                  *post_vol[FTNREF2D(j-1,k  ,x_max+5,x_min-2,y_min-2)]);
      }
    }
#pragma omp for private(j)
    for (k=y_min-1;k<=y_max+2;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                          -node_flux[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]+node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma omp for private(upwind,downwind,donor,dif,sigma,width,limiter,vdiffuw,vdiffdw,auw,adw,wind,j)
    for (k=y_min-1;k<=y_max+1;k++) {
      for (j=x_min;j<=x_max+1;j++) {
        if(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]<0.0){
          upwind=k+2;
          donor=k+1;
          downwind=k;
          dif=donor;
        }
        else{
          upwind=k-1;
          donor=k;
          downwind=k+1;
          dif=upwind;
        }
        sigma=fabs(node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])/(node_mass_pre[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]);
        width=celldy[FTNREF1D(k,y_min-2)];
        vdiffuw=vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,upwind,x_max+5,x_min-2,y_min-2)];
        vdiffdw=vel1[FTNREF2D(j  ,downwind ,x_max+5,x_min-2,y_min-2)]-vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)];
        limiter=0.0;
        if(vdiffuw*vdiffdw>0.0){
          auw=fabs(vdiffuw);
          adw=fabs(vdiffdw);
          wind=1.0;
          if(vdiffdw<=0.0) wind=-1.0;
          limiter=wind*MIN(width*((2.0-sigma)*adw/width+(1.0+sigma)*auw/celldy[FTNREF1D(dif,y_min-2)])/6.0,MIN(auw,adw));
        }
        advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=vel1[FTNREF2D(j  ,donor,x_max+5,x_min-2,y_min-2)]+(1.0-sigma)*limiter;
        mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=advec_vel[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                           *node_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }
#pragma omp for private(j)
    for (k=y_min;k<=y_max+1;k++) {
#pragma ivdep
      for (j=x_min;j<=x_max+1;j++) {
        vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]=(vel1[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        *node_mass_pre[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)]
                                                        +mom_flux[FTNREF2D(j  ,k-1,x_max+5,x_min-2,y_min-2)]
                                                        -mom_flux[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)])
                                                        /node_mass_post[FTNREF2D(j  ,k  ,x_max+5,x_min-2,y_min-2)];
      }
    }

  }

 }


}



