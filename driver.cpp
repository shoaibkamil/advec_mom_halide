#include "driver.h"

void create_grid(double** ptr, unsigned int size) {
  *ptr = (double*) malloc(sizeof(double)*size);
}


// fill a grid with a random double between 1 and 100.
void fill_grid_rand(double* ptr, unsigned int size) {
  for (int i=0; i<size; i++)
    ptr[i] = 1.0 + ((double)rand() / RAND_MAX) * 99.0;
}

void create_all_grids(double **xvel1,
                      double **yvel1,
                      double **mass_flux_x,
                      double **vol_flux_x,
                      double **mass_flux_y,
                      double **vol_flux_y,
                      double **volume,
                      double **density1,
                      double **node_flux,
                      double **node_mass_post,
                      double **node_mass_pre,
                      double **advec_vel,
                      double **mom_flux,
                      double **pre_vol,
                      double **post_vol,
                      double **celldx,
                      double **celldy,
                      unsigned int size)
{
  create_grid(xvel1, size);
  create_grid(yvel1, size);
  create_grid(mass_flux_x, size);
  create_grid(vol_flux_x, size);
  create_grid(mass_flux_y, size);
  create_grid(vol_flux_y, size);
  create_grid(volume, size);
  create_grid(density1, size);
  create_grid(node_flux, size);
  create_grid(node_mass_post, size);
  create_grid(node_mass_pre, size);
  create_grid(advec_vel, size);
  create_grid(mom_flux, size);
  create_grid(pre_vol, size);
  create_grid(post_vol, size);
  create_grid(celldx, size);
  create_grid(celldy, size);

  fill_grid_rand(*xvel1, size);
  fill_grid_rand(*yvel1, size);
  fill_grid_rand(*mass_flux_x, size);
  fill_grid_rand(*vol_flux_x, size);
  fill_grid_rand(*mass_flux_y, size);
  fill_grid_rand(*vol_flux_y, size);
  fill_grid_rand(*volume, size);
  fill_grid_rand(*density1, size);
  fill_grid_rand(*node_flux, size);
  fill_grid_rand(*node_mass_post, size);
  fill_grid_rand(*node_mass_pre, size);
  fill_grid_rand(*advec_vel, size);
  fill_grid_rand(*mom_flux, size);
  fill_grid_rand(*pre_vol, size);
  fill_grid_rand(*post_vol, size);
  fill_grid_rand(*celldx, size);
  fill_grid_rand(*celldy, size);

}


int main(int argc, char* argv[]) {

  // set integer parameters
  int xmin = 1;
  int ymin = 1;
  int xmax = 2048;
  int ymax = 2048;
  int which_vl = 1;
  int swp_nmbr = 1;
  int drctn = 1;

  // allocate arrays

  // don't know the actual size here
  // TODO: figure out size
  int size = 2*(xmax-xmin)*(ymax-ymin);

  double *xvel1;
  double *yvel1;
  double *mass_flux_x;
  double *vol_flux_x;
  double *mass_flux_y;
  double *vol_flux_y;
  double *volume;
  double *density1;
  double *node_flux;
  double *node_mass_post;
  double *node_mass_pre;
  double *advec_vel;
  double *mom_flux;
  double *pre_vol;
  double *post_vol;
  double *celldx;
  double *celldy;
  Buffer *b_volume;
  Buffer *b_vol_flux_y;
  Buffer *b_post_vol;
  Buffer *b_pre_vol;
  Buffer *b_vol_flux_x;
  Buffer *b_mass_flux_x;
  Buffer *b_node_flux;
  Buffer *b_density1;
  Buffer *b_node_mass_post;
  Buffer *b_node_mass_pre;
  Buffer *b_vel1;
  Buffer *b_celldx;
  Buffer *b_advec_vel;
  Buffer *b_mom_flux;


  Buffer *b_out_post_vol;
  Buffer *b_out_pre_vol;
  Buffer *b_out_node_flux;
  Buffer *b_out_node_mass_post;
  Buffer *b_out_node_mass_pre;
  Buffer *b_out_advec_vel;
  Buffer *b_out_mom_flux;
  Buffer *b_out_vel1;


  create_all_grids(&xvel1,
                    &yvel1,
                    &mass_flux_x,
                    &vol_flux_x,
                    &mass_flux_y,
                    &vol_flux_y,
                    &volume,
                    &density1,
                    &node_flux,
                    &node_mass_post,
                    &node_mass_pre,
                    &advec_vel,
                    &mom_flux,
                    &pre_vol,
                    &post_vol,
                    &celldx,
                    &celldy,
                    size);

  advec_mom_alloc_buffers(&xmin,&xmax,&ymin,&ymax,
                      xvel1,
                      yvel1,
                      mass_flux_x,
                      vol_flux_x,
                      mass_flux_y,
                      vol_flux_y,
                      volume,
                      density1,
                      node_flux,
                      node_mass_post,
                      node_mass_pre,
                      advec_vel,
                      mom_flux,
                      pre_vol,
                      post_vol,
                      celldx,
                      celldy,
                      &which_vl,
                      &swp_nmbr,
                      &drctn,

                      &b_volume,
                      &b_vol_flux_y,
                      &b_post_vol,
                      &b_pre_vol,
                      &b_vol_flux_x,
                      &b_mass_flux_x,
                      &b_node_flux,
                      &b_density1,
                      &b_node_mass_post,
                      &b_node_mass_pre,
                      &b_vel1,
                      &b_celldx,
                      &b_advec_vel,
                      &b_mom_flux,


                      &b_out_post_vol,
                      &b_out_pre_vol,
                      &b_out_node_flux,
                      &b_out_node_mass_post,
                      &b_out_node_mass_pre,
                      &b_out_advec_vel,
                      &b_out_mom_flux,
                      &b_out_vel1);


  // warm up icache
  advec_mom_kernel_c_(&xmin,&xmax,&ymin,&ymax,
                      xvel1,
                      yvel1,
                      mass_flux_x,
                      vol_flux_x,
                      mass_flux_y,
                      vol_flux_y,
                      volume,
                      density1,
                      node_flux,
                      node_mass_post,
                      node_mass_pre,
                      advec_vel,
                      mom_flux,
                      pre_vol,
                      post_vol,
                      celldx,
                      celldy,
                      &which_vl,
                      &swp_nmbr,
                      &drctn);

  advec_mom_kernel_halide(
                      b_volume,
                      b_vol_flux_y,
                      b_post_vol,
                      b_pre_vol,
                      b_vol_flux_x,
                      b_mass_flux_x,
                      b_node_flux,
                      b_density1,
                      b_node_mass_post,
                      b_node_mass_pre,
                      b_vel1,
                      b_celldx,
                      b_advec_vel,
                      b_mom_flux,


                      b_out_post_vol,
                      b_out_pre_vol,
                      b_out_node_flux,
                      b_out_node_mass_post,
                      b_out_node_mass_pre,
                      b_out_advec_vel,
                      b_out_mom_flux,
                      b_out_vel1);


  long long start, end;

  start = rdtsc();

  for (int i=0; i<5; i++) {

  advec_mom_kernel_c_(&xmin,&xmax,&ymin,&ymax,
                      xvel1,
                      yvel1,
                      mass_flux_x,
                      vol_flux_x,
                      mass_flux_y,
                      vol_flux_y,
                      volume,
                      density1,
                      node_flux,
                      node_mass_post,
                      node_mass_pre,
                      advec_vel,
                      mom_flux,
                      pre_vol,
                      post_vol,
                      celldx,
                      celldy,
                      &which_vl,
                      &swp_nmbr,
                      &drctn);
  }
  end = rdtsc();
  printf("Took %lld cycles\n", end-start);
  start = rdtsc();

  for (int i=0; i<5; i++) {
  advec_mom_kernel_halide(
                      b_volume,
                      b_vol_flux_y,
                      b_post_vol,
                      b_pre_vol,
                      b_vol_flux_x,
                      b_mass_flux_x,
                      b_node_flux,
                      b_density1,
                      b_node_mass_post,
                      b_node_mass_pre,
                      b_vel1,
                      b_celldx,
                      b_advec_vel,
                      b_mom_flux,


                      b_out_post_vol,
                      b_out_pre_vol,
                      b_out_node_flux,
                      b_out_node_mass_post,
                      b_out_node_mass_pre,
                      b_out_advec_vel,
                      b_out_mom_flux,
                      b_out_vel1);

  }
  end = rdtsc();
  printf("Took %lld cycles\n", end-start);

  return 0;

}
