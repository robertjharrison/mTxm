TARGET_ARCHITECTURE=neon
TUNING_DATA=tune-graviton.txt

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
#MKLROOT=/gpfs/software/intel/parallel-studio-xe/2018/mkl
#CXX = g++
#CXXFLAGS = -Ofast -mtune=native -ffast-math -mavx -mavx2 -mavx512f -mfma -I. -std=c++17 -DUSE_LIBXSMM -Wall -Wno-unused-but-set-variable -m64 -I/gpfs/home/rharrison/libxsmm/include
#LDFLAGS =  -L/gpfs/home/rharrison/libxsmm/lib -lxsmm -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

# Intel compiler for AVX2 (laptop)
#CXX = icpc
#CXXFLAGS = -Ofast -mtune=native -ffast-math -march=core-avx2 -I. -std=c++17 -DUSE_AVX2 -Wall -Wno-unused-but-set-variable -m64
#LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpapi -lpthread -lm -ldl


# ARM compiler for graviton
CXXFLAGS = -g -Ofast -ffp-contract=fast -ffast-math -Wall -mcpu=native -I. -std=c++17 -I/opt/arm/armpl-19.2.0_Cortex-A72_Ubuntu-16.04_arm-hpc-compiler_19.2_aarch64-linux/include -DUSE_NEON
CXX = armclang++
LDFLAGS=-L/opt/arm/armpl-19.2.0_Cortex-A72_Ubuntu-16.04_arm-hpc-compiler_19.2_aarch64-linux/lib -larmpl -lamath
PYTHON=python3

# g++ compiler for graviton
#CXXFLAGS = -Ofast -ffast-math -Wall -mtune=native -mcpu=native -I. -std=c++17 -I/opt/arm/armpl-19.2.0_Cortex-A72_Ubuntu-16.04_arm-hpc-compiler_19.2_aarch64-linux/include 
#CXX = g++
#LDFLAGS=-L/opt/arm/armpl-19.2.0_Cortex-A72_Ubuntu-16.04_arm-hpc-compiler_19.2_aarch64-linux/lib -larmpl -lamath
#PYTHON=python3

TARGETS = tester tuner bencher

all:	$(TARGETS)

tester:	tester.o mTxmq.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

tuner:	tuner.o mTxmq.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

bencher:	bencher.o mTxmq.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

mTxmq.h mTxmq.cc:	gen.py avx2.py avx512.py neon.py tune.py
	$(PYTHON) gen.py $(TARGET_ARCHITECTURE) $(TUNING_DATA) > mTxmq.cc

mTxmq.o:	mTxmq.cc mTxmq.h

tester.o tuner.o bencher.o:	mTxmq.h 

asm:	mTxmq.cc
	$(CXX) $(CXXFLAGS) -S -c mTxmq.cc

clean:
	/bin/rm -rf $(TARGETS) *.o mTxmq.cc mTxmq.h *~ __pycache__  mTxmq.s
