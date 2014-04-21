all: 
	g++ -I ../Halide/include driver.cpp gen_halide_func.cpp advec_mom_kernel_c.c advec_mom_halide.cpp -o advec_mom -L ../Halide/bin -lHalide f_pre_vol.o f_post_vol.o -lpthread

fast:
	g++ -static -O3 -fopenmp -I ../Halide/include driver.cpp gen_halide_func.cpp advec_mom_kernel_c.c advec_mom_halide.cpp -o advec_mom -L ../Halide/bin -lHalide f_pre_vol.o f_post_vol.o -lpthread `llvm-config --ldflags --libs`


