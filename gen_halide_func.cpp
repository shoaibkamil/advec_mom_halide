#include <Halide.h>

using namespace Halide;

void gen_halide_func() {

  Var j("j"), k("k");

  // Inputs
/*  Image<double> volume;
  Image<double> vol_flux_y;
  Image<double> vol_flux_x;
  Image<double> mass_flux_x;
  Image<double> density1;
  Image<double> celldx;
  Image<double> vel1;

  Func post_vol("post_vol");
  Func pre_vol("pre_vol");
  Func node_flux("node_flux");
  Func node_mass_post("node_mass_post");
  Func node_mass_pre("node_mass_pre");
  Func advec_vel("advec_vel");
  Func mom_flux("mom_flux");
*/
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

  Func f_post_vol("f_post_vol");
  Func f_pre_vol("f_pre_vol");
  Func f_node_flux("f_node_flux");
  Func f_node_mass_post("f_node_mass_post");
  Func f_node_mass_pre("f_node_mass_pre");
  Func f_advec_vel("f_advec_vel");
  Func f_mom_flux("f_mom_flux");
  Func f_vel1("f_vel1");

  // if mom_sweep == 1
  f_post_vol(j, k) = volume(j,k) + vol_flux_y(j,clamp(k+1, 0, vol_flux_y.height()-1)) - vol_flux_y(j,k);

  f_pre_vol(j,k) = f_post_vol(j,k) + vol_flux_x(clamp(j+1, 0, vol_flux_x.width()-1),k) - vol_flux_x(j,k);

  f_post_vol.compile_to_file("f_post_vol", volume, vol_flux_y);
  f_pre_vol.compile_to_file("f_pre_vol", vol_flux_x, vol_flux_y, volume);
  

  // if direction == 1
  Func mass_flux_x_clamped("mass_flux_x_clamped");
  Func density1_clamped("density1_clamped");
  Func f_post_vol_clamped("f_post_vol_clamped");
  Func f_node_mass_post_clamped("f_node_mass_post_clamped");
  Func f_node_flux_clamped("f_node_flux_clamped");

  mass_flux_x_clamped(j,k) = mass_flux_x(clamp(j, 1, mass_flux_x.width()-1), clamp(k, 1, mass_flux_x.height()-1));
  density1_clamped(j,k) = density1(clamp(j, 1, density1.width()-1), clamp(k, 1, density1.height()-1));
  //f_post_vol_clamped(j,k) = f_post_vol(clamp(j, 1, f_post_vol.width()-1), clamp(k, 1, f_post_vol.height()-1));
  //f_node_mass_post_clamped(j,k) = f_node_mass_post(clamp(j, 1, f_node_mass_post.width()-1), clamp(k, 1, f_node_mass_post.height()-1));
  //f_node_flux_clamped(j,k) = f_node_flux(clamp(j, 1, f_node_flux.width()-1), clamp(k, 1, f_node_flux.height()-1));

  f_node_flux(j,k) = 0.25f * (mass_flux_x_clamped(j,k-1)
                         + mass_flux_x_clamped(j,k)
                         + mass_flux_x_clamped(j+1,k-1)
                         + mass_flux_x_clamped(j+1,k));
  f_node_mass_post(j,k) = 0.25f * (density1_clamped(j,k-1) * f_post_vol(j,k-1)
                              + density1_clamped(j,k) * f_post_vol(j,k)
                              + density1_clamped(j-1,k-1) * f_post_vol(j-1,k-1)
                              + density1_clamped(j-1,k) * f_post_vol(j-1,k));
  f_node_mass_pre(j,k) = f_node_mass_post(j,k) - f_node_flux(j-1,k) + f_node_flux(j,k);

  f_node_flux.compile_to_file("f_node_flux", mass_flux_x);
  f_node_mass_post.compile_to_file("f_node_mass_post", density1, volume, vol_flux_y);
  f_node_mass_pre.compile_to_file("f_node_mass_pre", density1, volume, vol_flux_y, mass_flux_x);


  return;
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

  f_advec_vel(j,k) = vel1(donor,k) + (1.0f-sigma)*limiter;
  f_mom_flux(j,k) = f_advec_vel(j,k) * f_node_flux(j,k);

  f_vel1(j,k) = (vel1(j,k) * f_node_mass_pre(j,k)
             + f_mom_flux(j-1, k)
             - f_mom_flux(j,k)) / f_node_mass_post(j,k);

  f_advec_vel.compile_to_file("f_advec_vel", vel1, celldx, mass_flux_x, density1, volume);


}
