HALIDE_DIR=../Halide
INCLUDES=-I${HALIDE_DIR}/include
LIBS=-L${HALIDE_DIR}/bin -lHalide -lpthread `llvm-config --ldflags --libs`

fast:
	g++ -O3 -fopenmp ${INCLUDES} driver.cpp  advec_mom_kernel_c.c advec_mom_halide.cpp -o advec_mom advec_mom_halide_gen.o ${LIBS}

staticfast:
	g++ -static -O3 -fopenmp ${INCLUDES} driver.cpp  advec_mom_kernel_c.c advec_mom_halide.cpp -o advec_mom advec_mom_halide_gen.o ${LIBS}

gen: gen_halide_func.cpp
	g++ ${INCLUDES} gen_halide_func.cpp -o gen ${LIBS}
