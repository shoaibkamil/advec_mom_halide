#include <Halide.h>
#include <cstdio>
#include "advec_mom_halide_gen.h"
#include "ftocmacros.h"

using namespace Halide;

void advec_mom_alloc_buffers(int *xmin,int *xmax,int *ymin,int *ymax,
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
                         int *drctn,

                      Buffer **b_volume,
                      Buffer **b_vol_flux_y,
                      Buffer **b_post_vol,
                      Buffer **b_pre_vol,
                      Buffer **b_vol_flux_x,
                      Buffer **b_mass_flux_x,
                      Buffer **b_node_flux,
                      Buffer **b_density1,
                      Buffer **b_node_mass_post,
                      Buffer **b_node_mass_pre,
                      Buffer **b_vel1,
                      Buffer **b_celldx,
                      Buffer **b_advec_vel,
                      Buffer **b_mom_flux,


                      Buffer **b_out_post_vol,
                      Buffer **b_out_pre_vol,
                      Buffer **b_out_node_flux,
                      Buffer **b_out_node_mass_post,
                      Buffer **b_out_node_mass_pre,
                      Buffer **b_out_advec_vel,
                      Buffer **b_out_mom_flux,
                      Buffer **b_out_vel1) {



  double *vel1;

  if(*whch_vl==1){
    vel1=xvel1;
  } else{
    vel1=yvel1;
  }
  *b_volume = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)volume);
  *b_vol_flux_y = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)vol_flux_y);
  *b_post_vol = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)post_vol);
  *b_pre_vol = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)pre_vol);
  *b_vol_flux_x = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)vol_flux_x);
  *b_mass_flux_x = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)mass_flux_x);
  *b_node_flux = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)node_flux);
  *b_density1 = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)density1);
  *b_node_mass_post = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)node_mass_post);
  *b_node_mass_pre = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)node_mass_pre);
  *b_vel1 = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)vel1);
  *b_celldx = new Buffer(Float(64), *xmax-*xmin, 0, 0, 0, (uint8_t*)celldx);
  *b_advec_vel = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)advec_vel);
  *b_mom_flux = new Buffer(Float(64), *xmax-*xmin, *ymax-*ymin, 0, 0, (uint8_t*)mom_flux);


  *b_out_post_vol = new Buffer(Float(64), *xmax-*xmin-4, *ymax-*ymin-4, 0, 0, (uint8_t*)post_vol);
  *b_out_pre_vol = new Buffer(Float(64), *xmax-*xmin-4, *ymax-*ymin-4, 0, 0, (uint8_t*)pre_vol);
  *b_out_node_flux = new Buffer(Float(64), *xmax-*xmin-4, *ymax-*ymin-4, 0, 0, (uint8_t*)node_flux);
  *b_out_node_mass_post = new Buffer(Float(64), *xmax-*xmin-4, *ymax-*ymin-4, 0, 0, (uint8_t*)node_mass_post);
  *b_out_node_mass_pre = new Buffer(Float(64), *xmax-*xmin-4, *ymax-*ymin-4, 0, 0, (uint8_t*)node_mass_pre);
  *b_out_advec_vel = new Buffer(Float(64), *xmax-*xmin-4, *ymax-*ymin-4, 0, 0, (uint8_t*)advec_vel);
  *b_out_mom_flux = new Buffer(Float(64), *xmax-*xmin-4, *ymax-*ymin-4, 0, 0, (uint8_t*)mom_flux);
  *b_out_vel1 = new Buffer(Float(64), *xmax-*xmin-4, *ymax-*ymin-4, 0, 0, (uint8_t*)vel1);

  (*b_out_post_vol)->set_min(2,2);
  (*b_out_pre_vol)->set_min(2,2);
  (*b_out_node_flux)->set_min(2,2);
  (*b_out_node_mass_post)->set_min(2,2);
  (*b_out_node_mass_pre)->set_min(2,2);
  (*b_out_advec_vel)->set_min(2,2);
  (*b_out_mom_flux)->set_min(2,2);
  (*b_out_vel1)->set_min(2,2);


}

void advec_mom_kernel_halide(
                      Buffer *b_volume,
                      Buffer *b_vol_flux_y,
                      Buffer *b_post_vol,
                      Buffer *b_pre_vol,
                      Buffer *b_vol_flux_x,
                      Buffer *b_mass_flux_x,
                      Buffer *b_node_flux,
                      Buffer *b_density1,
                      Buffer *b_node_mass_post,
                      Buffer *b_node_mass_pre,
                      Buffer *b_vel1,
                      Buffer *b_celldx,
                      Buffer *b_advec_vel,
                      Buffer *b_mom_flux,


                      Buffer *b_out_post_vol,
                      Buffer *b_out_pre_vol,
                      Buffer *b_out_node_flux,
                      Buffer *b_out_node_mass_post,
                      Buffer *b_out_node_mass_pre,
                      Buffer *b_out_advec_vel,
                      Buffer *b_out_mom_flux,
                      Buffer *b_out_vel1) {

  advec_mom_halide_gen(b_volume->raw_buffer(), b_vol_flux_y->raw_buffer(),
      b_vol_flux_x->raw_buffer(), b_mass_flux_x->raw_buffer(),
      b_density1->raw_buffer(), b_celldx->raw_buffer(), b_vel1->raw_buffer(),
      b_post_vol->raw_buffer(), b_pre_vol->raw_buffer(),
      b_node_flux->raw_buffer(), b_node_mass_post->raw_buffer(),
      b_node_mass_pre->raw_buffer(), b_advec_vel->raw_buffer(),
      b_mom_flux->raw_buffer(),

      b_out_post_vol->raw_buffer(), b_out_pre_vol->raw_buffer(),
      b_out_node_flux->raw_buffer(), b_out_node_mass_post->raw_buffer(),
      b_out_node_mass_pre->raw_buffer(), b_out_advec_vel->raw_buffer(),
      b_out_mom_flux->raw_buffer(), b_out_vel1->raw_buffer());
  return;


}



