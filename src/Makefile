# Intel compiler for AVX512 (snmem)
#CXX = icpc
#CXXFLAGS = -Ofast -mtune=native -ffast-math -mavx -mavx2 -mfma -mavx512f -I. -std=c++17 -DUSE_AVX512 -Wall -Wno-unused-but-set-variable -m64 -I${MKLROOT}/include
#LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -lpthread -lm -ldl


# g++-9 for AVX2 (laptop) with mkl
#CXX = g++-9
#CXXFLAGS = -Ofast -mtune=native -ffast-math -mavx -mavx2 -mfma -I. -std=c++17 -DUSE_AVX2 -Wall -Wno-unused-but-set-variable -m64 -I${MKLROOT}/include
#LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpapi -lpthread -lm -ldl


# g++-9 for AVX2 (laptop) with libxsmm
#CXX = g++-9
##CXXFLAGS = -O0 -g -mtune=native -ffast-math -mavx -mavx2 -mfma -I. -std=c++17 -DUSE_AVX2 -DUSE_LIBXSMM -Wall -Wno-unused-but-set-variable -m64 -I/home/rjh/Downloads/libxsmm/include
#CXXFLAGS = -Ofast -mtune=native -ffast-math -mavx -mavx2 -mfma -I. -std=c++17 -DUSE_AVX2 -DUSE_LIBXSMM -Wall -Wno-unused-but-set-variable -m64 -I/home/rjh/Downloads/libxsmm/include
#LDFLAGS =  -L/home/rjh/Downloads/libxsmm/lib -lxsmm -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpapi -lpthread -lm -ldl


# g++-9 for AVX2 (haswell seawulf) with libxsmm
#MKLROOT=/gpfs/software/intel/parallel-studio-xe/2018/mkl
#CXX = g++
#CXXFLAGS = -Ofast -mtune=native -ffast-math -mavx -mavx2 -mfma -I. -std=c++17 -DUSE_AVX2 -DUSE_LIBXSMM -Wall -Wno-unused-but-set-variable -m64 -I/gpfs/home/rharrison/libxsmm/include
#LDFLAGS =  -L/gpfs/home/rharrison/libxsmm/lib -lxsmm -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

# g++-9 for AVX512 (snmem seawulf) with libxsmm
MKLROOT=/gpfs/software/intel/parallel-studio-xe/2018/mkl
CXX = g++
CXXFLAGS = -Ofast -mtune=native -ffast-math -mavx -mavx2 -mavx512f -mfma -I. -std=c++17 -DUSE_LIBXSMM -Wall -Wno-unused-but-set-variable -m64 -I/gpfs/home/rharrison/libxsmm/include
LDFLAGS =  -L/gpfs/home/rharrison/libxsmm/lib -lxsmm -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

# Intel compiler for AVX2 (laptop)
#CXX = icpc
#CXXFLAGS = -Ofast -mtune=native -ffast-math -march=core-avx2 -I. -std=c++17 -DUSE_AVX2 -Wall -Wno-unused-but-set-variable -m64
#LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpapi -lpthread -lm -ldl


# ARM compiler for graviton
#CXXFLAGS = -Ofast -ffp-contract=fast -ffast-math -Wall -mcpu=native -I. -std=c++17 -I/opt/arm/armpl-19.2.0_Cortex-A72_Ubuntu-16.04_arm-hpc-compiler_19.2_aarch64-linux/include -DUSE_NEON
#CXX = armclang++
#LDFLAGS=-L/opt/arm/armpl-19.2.0_Cortex-A72_Ubuntu-16.04_arm-hpc-compiler_19.2_aarch64-linux/lib -larmpl -lamath

test:	test.o mTxmq.o timerstuff.h
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

mTxmq.h mTxmq.cc:	gen.py gen_avx2.py gen_avx512.py gen_neon.py
	python gen.py > mTxmq.cc

test.o:	mTxmq.h 

asm:	mTxmq.cc
	$(CXX) $(CXXFLAGS) -S -c mTxmq.cc

clean:
	/bin/rm -rf test *.o mTxmq.cc mTxmq.h *~ __pycache__ 
