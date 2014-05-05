#include <Halide.h>
#include <iostream>

using namespace Halide;

/*
 * We only currently generate the case where mom_sweep=2 and direction=1.
 */
int main(void) {

  Var j("j"), k("k");

  ImageParam volume(Float(64), 2);                      // read-only
  ImageParam vol_flux_y(Float(64), 2);                  // read-only
  ImageParam vol_flux_x(Float(64), 2);                  // read-only
  ImageParam mass_flux_x(Float(64), 2);                 // read-only
  ImageParam density1(Float(64), 2);                    // read-only
  ImageParam celldx(Float(64), 1);                      // read-only
  ImageParam vel1(Float(64), 2);                        // read-write
  ImageParam post_vol(Float(64), 2);                    // read-write
  ImageParam pre_vol(Float(64), 2);                     // read-write
  ImageParam node_flux(Float(64), 2);                   // read-write
  ImageParam node_mass_post(Float(64), 2);              // read-write
  ImageParam node_mass_pre(Float(64), 2);               // read-write
  ImageParam advec_vel(Float(64), 2);                   // read-write
  ImageParam mom_flux(Float(64), 2);                    // read-write

  /* Not sure if right place. */
  post_vol.set_min(2,2);
  pre_vol.set_min(2,2);
  node_flux.set_min(2,2);
  node_mass_pre.set_min(2,2);
  node_mass_post.set_min(2,2);
  advec_vel.set_min(2,2);
  mom_flux.set_min(2,2);
  vel1.set_min(2,2);
  volume.set_min(2,2);
  vol_flux_x.set_min(2,2);
  vol_flux_y.set_min(2,2);
  density1.set_min(2,2);
  mass_flux_x.set_min(2,2);


  Func f_post_vol("f_post_vol");
  Func f_pre_vol("f_pre_vol");
  Func f_node_flux("f_node_flux");
  Func f_node_mass_post("f_node_mass_post");
  Func f_node_mass_pre("f_node_mass_pre");
  Func f_advec_vel("f_advec_vel");
  Func f_mom_flux("f_mom_flux");
  Func f_vel1("f_vel1");


  /* 
   * The following are all straightforward translations from the Fortran 
   * We use expressions and functions both, so we name e_* for expressions
   * f_* for functions.  Everything else directly accesses arrays or is an
   * intermediate quantity.
   */

  Expr e_post_vol = volume(j,k) + vol_flux_y(j,k+1) - vol_flux_y(j,k);
  f_post_vol(j,k) = e_post_vol;


  Expr e_pre_vol = f_post_vol(j,k) + vol_flux_x(j+1,k) - vol_flux_x(j,k);
  f_pre_vol(j,k) = e_pre_vol;

  Expr e_node_flux = 0.25f * (mass_flux_x(j,k-1)
                         + mass_flux_x(j,k)
                         + mass_flux_x(j+1,k-1)
                         + mass_flux_x(j+1,k));
  f_node_flux(j,k) = e_node_flux;

  Expr e_node_mass_post = 0.25f * (density1(j,k-1) * f_post_vol(j,k-1)
                              + density1(j,k) * f_post_vol(j,k)
                              + density1(j-1,k-1) * f_post_vol(j-1,k-1)
                              + density1(j-1,k) * f_post_vol(j-1,k));
  f_node_mass_post(j,k) = e_node_mass_post;

  Expr e_node_mass_pre = f_node_mass_post(j,k) - f_node_flux(j-1,k) + f_node_flux(j,k);
  f_node_mass_pre(j,k) = e_node_mass_pre;




  Expr upwind = select(f_node_flux(j,k) < 0.0f, j+2, j-1);
  Expr donor = select(f_node_flux(j,k) < 0.0f, j+1, j);
  Expr downwind = select(f_node_flux(j,k) < 0.0f, j, j+1);
  Expr dif = select(f_node_flux(j,k) < 0.0f, donor, upwind);

  Expr sigma = abs(f_node_flux(j,k)) / f_node_mass_pre(donor, k);
  Expr width = celldx(j);
  Expr vdiffuw = vel1(donor,k) - vel1(upwind,k);
  Expr vdiffdw = vel1(downwind,k) - vel1(donor,k);
  
  Expr auw = abs(vdiffuw);
  Expr adw = abs(vdiffdw);
  Expr wind = select(vdiffdw <= 0.0f, -1.0f, 1.0f);

  
  Expr limiter = select(vdiffuw*vdiffdw > 0.0f, wind*min(width*((2.0f-sigma)*adw/width+(1.0f+sigma)*auw/celldx(dif))/6.0f,
                                                         min(auw,adw)),
                                                cast(Float(64), 0.0f));

  Expr e_advec_vel = vel1(donor,k) + (1.0f-sigma)*limiter;
  f_advec_vel(j,k) = e_advec_vel;

  Expr e_mom_flux = f_advec_vel(j,k) * f_node_flux(j,k);
  f_mom_flux(j,k) = e_mom_flux;

  Expr e_vel1 = (vel1(j,k) * f_node_mass_pre(j,k)
             + f_mom_flux(j-1, k)
             - f_mom_flux(j,k)) / f_node_mass_post(j,k);


  /*
   * End of all expressions.  Now, we build the argument vector.
   */
  std::vector<Argument> args;
  args.push_back(volume);
  args.push_back(vol_flux_y);
  args.push_back(vol_flux_x);
  args.push_back(mass_flux_x);
  args.push_back(density1);
  args.push_back(celldx);
  args.push_back(vel1);
  args.push_back(post_vol);
  args.push_back(pre_vol);
  args.push_back(node_flux);
  args.push_back(node_mass_post);
  args.push_back(node_mass_pre);
  args.push_back(advec_vel);
  args.push_back(mom_flux);

  /*
   * Because there are multiple output arrays, we build a tuple,
   * which Halide uses as output for such a function.
   */
  std::vector<Expr> rhs;
  rhs.push_back(e_post_vol);
  rhs.push_back(e_pre_vol);
  rhs.push_back(e_node_flux);
  rhs.push_back(e_node_mass_post);
  rhs.push_back(e_node_mass_pre);
  rhs.push_back(e_advec_vel);
  rhs.push_back(e_mom_flux);
  rhs.push_back(e_vel1);

  Func advec_mom("advec_mom");
  advec_mom(j,k) = Tuple(rhs);
  AUTOTUNE_HOOK(advec_mom);

  advec_mom.vectorize(j, 4);
  advec_mom.parallel(k);
  BASELINE_HOOK(advec_mom);

  /*
   * Generate and save a file we can link against later.
   */
  advec_mom.compile_to_file("advec_mom_halide_gen", args);


  return 0;
}
