#include <cstdlib>
#include <cstdio>
#include <Halide.h>

#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}

#elif defined(__x86_64__)

static __inline__ unsigned long long rdtsc(void)
{
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#endif

void advec_mom_kernel_c_(int *xmin,int *xmax,int *ymin,int *ymax,
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
                         int *drctn);

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
                         int *drctn);


void create_grid(double** ptr, unsigned int size) {
  *ptr = (double*) malloc(sizeof(double)*size);
}

void gen_halide_func();

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
  int xmax = 960;
  int ymax = 960;
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

  gen_halide_func();

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


  long long start, end;

  start = rdtsc();

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

  end = rdtsc();
  printf("Took %lld cycles\n", end-start);
  start = rdtsc();
  advec_mom_kernel_halide(&xmin,&xmax,&ymin,&ymax,
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


  end = rdtsc();
  printf("Took %lld cycles\n", end-start);

  return 0;

}
