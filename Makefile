all: 
	g++ -I ../Halide/include driver.cpp  advec_mom_kernel_c.c advec_mom_halide.cpp -o advec_mom -L ../Halide/bin -lHalide f_pre_vol.o f_post_vol.o -lpthread

fast:
	g++ -O3 -fopenmp -I ../Halide/include driver.cpp  advec_mom_kernel_c.c advec_mom_halide.cpp -o advec_mom -L ../Halide/bin -lHalide f_pre_vol.o f_post_vol.o -lpthread `llvm-config --ldflags --libs`

staticfast:
	g++ -static -O3 -fopenmp -I ../Halide/include driver.cpp  advec_mom_kernel_c.c advec_mom_halide.cpp -o advec_mom -L ../Halide/bin -lHalide f_pre_vol.o f_post_vol.o -lpthread `llvm-config --ldflags --libs`

gen: gen_halide_func.cpp
	g++ -I ../Halide/include gen_halide_func.cpp -L ../Halide/bin -o gen -lHalide
