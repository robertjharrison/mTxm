#if defined(USE_LIBXSMM)
#include <libxsmm.h>
#define dgemm_ libxsmm_dgemm
#elif defined(USE_ARMPL)
#else
using armpl_int_t = int;
extern "C" void dgemm_(const char *transa, const char *transb, const armpl_int_t *m, const armpl_int_t *n, const armpl_int_t *k, const double *alpha, const double *a, const armpl_int_t *lda, const double *b, const armpl_int_t *ldb, const double *beta, double *c, const armpl_int_t *ldc, ... );
#endif

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <list>
#include <sched.h>

#include "mTxmq.h"

double DFLOPS=20e9;

#include "timerstuff.h"

using namespace std;

static const int ALIGNMENT=4096; //128;

double ran()
{
    static unsigned int seed = 76521;
    seed = seed *1812433253 + 12345;
    return ((double) (seed & 0x7fffffff)) * 4.6566128752458e-10;
}

void ran_fill(int n, double *a) {
    while (n--) *a++ = ran();
}


void mTxm_blas(int dimi, int dimj, int dimk, double* c, const double* a, const double*b, bool ttt=false) {
#if defined(USE_LIBXSMM)
    double one = 1.0;
    if (ttt) {
        libxsmm_dgemm("n","t",&dimi,&dimj,&dimk,&one,a,&dimi,b,&dimj,&one,c,&dimi);
    }  
    else {
        libxsmm_dgemm("n","t",&dimj,&dimi,&dimk,&one,b,&dimj,a,&dimi,&one,c,&dimj);
    }
#else
    double one = 1.0;
    if (ttt) {
        dgemm_("n","t",&dimi,&dimj,&dimk,&one,a,&dimi,b,&dimj,&one,c,&dimi,1,1);
    }  
    else {
        dgemm_("n","t",&dimj,&dimi,&dimk,&one,b,&dimj,a,&dimi,&one,c,&dimj,1,1);
    }
#endif
}

void mTxm(int dimi, int dimj, int dimk,
          double* __restrict__ c, const double* a, const double* b, bool transposeC=false) {
    if (transposeC) {
        for (int k=0; k<dimk; ++k) {
            for (int j=0; j<dimj; ++j) {
                for (int i=0; i<dimi; ++i) {
                    c[j*dimi+i] += a[k*dimi+i]*b[k*dimj+j];
                }
            }
        }
    }
    else {
        for (int k=0; k<dimk; ++k) {
            for (int j=0; j<dimj; ++j) {
                for (int i=0; i<dimi; ++i) {
                    c[i*dimj+j] += a[k*dimi+i]*b[k*dimj+j];
                }
            }
        }
    }
}



//#pragma GCC optimize "1"
template <typename funcT>
void timer(const char* s, funcT mTxmq, int ni, int nj, int nk, double *a, double *b, double *c) {
    double fastest=0.0, fastest_blas=0.0;

    double scale = 1e-9;
    double target = 1e-3; // was 1e-3
    
    double nflop = 2.0*ni*nj*nk;
    double est = 1e-6 + nflop/DFLOPS;
    int ntrial=10; // was 100

    // Ajust repetition count to hit target duration
    int nloop=std::max(1,int(target/est));
    for (int rep=0; rep<2; rep++) {
      //{ //while (true) {
      PerfData perf;
      for (int loop=0; loop<nloop; ++loop) {
	mTxmq(ni,nj,nk,c,a,b,-1,16);
      }
      perf.stop();
      double used = perf.cpu_time();
      nloop = std::max(1,int(nloop*target/used));
      // if (used > 2*target && nloop>2) nloop /= 2;
      // else if (used < 0.5*target) nloop *= 2;
      // else break;
    }					

    for (int t=0; t<ntrial; t++) {
        double rate;
        PerfData perf;
        for (int loop=0; loop<nloop; ++loop) {
            mTxmq(ni,nj,nk,c,a,b,-1,16);
        }
        perf.stop();
        double used = perf.cpu_time()/nloop;
        rate = scale*nflop/used;
        if (rate > fastest) fastest = rate;
    }
    
    for (int t=0; t<ntrial; t++) {
        double rate;
        PerfData perf;
        for (int loop=0; loop<nloop; ++loop) {
            mTxm_blas(nj,ni,nk,c,b,a,false);
        }
        perf.stop();
        double used = perf.cpu_time()/nloop;
        rate = scale*nflop/used;
        if (rate > fastest_blas) fastest_blas = rate;
    }
    
    printf("%20s %3d %3d %3d %8.2f %8.2f\n",s, ni,nj,nk, fastest, fastest_blas);
}

template <typename funcT>
void timer(funcT mTxmq) {
    const int len = 32768;
    
    double *a, *b, *c, *d;  
    if(posix_memalign((void **) &a, ALIGNMENT, len*sizeof(double))) exit(1);
    if(posix_memalign((void **) &b, ALIGNMENT, len*sizeof(double))) exit(1);
    if(posix_memalign((void **) &c, ALIGNMENT, len*sizeof(double))) exit(1);
    if(posix_memalign((void **) &d, ALIGNMENT, len*sizeof(double))) exit(1);
    ran_fill(len,a);
    ran_fill(len,b);
    ran_fill(len,c);
    ran_fill(len,d);
    
    printf("%20s %3s %3s %3s %8s (GF/s)\n", "type", "M", "N", "K", "LOOP");
    
    timer("(144,12)T*(12,12)", mTxmq, 144, 12, 12,a,b,c);
    timer("(12,12)T*(12,144)", mTxmq, 12, 144, 12,a,b,c);
    for (int ni=1; ni<=64; ni+=1) timer("(m,m)T*(m,m)", mTxmq, ni,ni,ni,a,b,c);
    for (int m=2; m<=32; m+=1) timer("(m,m*m)T*(m,m)", mTxmq, m*m,m,m,a,b,c);
    for (int m=2; m<=32; m+=1) timer("(m,m)T*(m,m*m)", mTxmq, m,m*m,m,a,b,c);
    for (int m=1; m<=20; m+=1) timer("(20,20*20)T*(20,m)", mTxmq, 20*20,m,20,a,b,c);
    for (int m=1; m<=20; m+=1) timer("(20,m)T*(20,20*20)", mTxmq, m,20*20,20,a,b,c);
    
    free(a);
    free(b);
    free(c);
    free(d);
}


int main() {

#ifdef USE_LIBXSMM
    libxsmm_init();
#endif

    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(1,&mask);
    if (sched_setaffinity(0, sizeof(mask), &mask) == -1) {
        perror("system error message");
        std::cout << "ThreadBase: set_affinity: Could not set cpu affinity" << std::endl;
    }
    auto num = std::chrono::high_resolution_clock::period::num;
    auto den = std::chrono::high_resolution_clock::period::den;
    std::cout << "period " << double(num)/den << std::endl;
    {
      PerfData p;
      double freq = p.frequency();
      std::cout << "frequency (GHz) " << freq << std::endl;
    }

    try {
      timer(&mTxmqG);
    }
    catch (const char* c) {
      std::cout << "exception " << c << std::endl;
      return 1;
    }
    
    return 0;
}
