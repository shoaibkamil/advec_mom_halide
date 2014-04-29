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

using namespace Halide;

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
                      Buffer *b_out_vel1);

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
                      Buffer **b_out_vel1);


