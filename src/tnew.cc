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


void printmat(const char* msg, const double* p, int dimi, int dimj) {
    printf("%s\n", msg);
    for (int i=0; i<dimi; i++) {
        for (int j=0; j<dimj; j++) {
            printf("%11.7f ", *p++);
        }
        printf("\n");
    }
}

template <typename funcT>
void testone(funcT mTxmq, int ni, int nj, int nk, bool debugging, bool trace) {
  double *a, *b, *c, *d, *e;  
  if (debugging || trace) std::cout << " testone " << ni << " " << nj << " " << nk << std::endl;
  if(posix_memalign((void **) &a, ALIGNMENT, nk*ni*sizeof(double))) throw "a";
  if(posix_memalign((void **) &b, ALIGNMENT, nk*nj*sizeof(double))) throw "b";
  if(posix_memalign((void **) &c, ALIGNMENT, ni*nj*sizeof(double))) throw "c";
  if(posix_memalign((void **) &d, ALIGNMENT, ni*nj*sizeof(double))) throw "d";
  if(posix_memalign((void **) &e, ALIGNMENT, ni*nj*sizeof(double))) throw "e";
  ran_fill(nk*ni, a);
  ran_fill(nk*nj, b);
  
  for (int i=0; i<ni*nj; ++i) d[i] = c[i] = e[i] = 0.0;
  
  if (debugging) {
    printmat("a", a, nk, ni);
    printmat("b", b, nk, nj);
  }

  mTxm (ni,nj,nk,c,a,b,false);

  set_kernel_trace(trace);
  if (trace) std::cout << "normal\n";
  mTxmq(ni,nj,nk,d,a,b,-1,16);
  // if (trace) std::cout << "transpose\n";
  // mTxmq(nj,ni,nk,e,b,a,true);
  set_kernel_trace(false);
  
  if (debugging) {
    printmat("c", c, ni, nj);
    printmat("d", d, ni, nj);
    printmat("e", e, ni, nj);
  }
  
  for (int i=0; i<ni*nj; ++i) {
    double err = abs(d[i]-c[i]);
    if (err > 1e-13) {
      printf("test_mtxmq: error %d %d %d (%d,%d) %.6e %.6e %.2e\n",ni,nj,nk,i/nj,i-(i/nj)*nj,c[i],d[i],err);
      throw "f";
    }
    // err = abs(e[i]-c[i]);
    // if (err > 1e-13) {
    //   printf("test_mtxmq:Terror %d %d %d (%d,%d) %.6e %.6e %.2e\n",ni,nj,nk,i/nj,i-(i/nj)*nj,c[i],e[i],err);
    //   exit(1);
    // }
  }
  if (debugging) exit(0);
  
  free(a);
  free(b);
  free(c);
  free(d);
  free(e);
}

template <typename funcT>
void tester(funcT mTxmq) {
    bool smalltest=true;
    const int nimax=!smalltest ? 128 : 64;
    const int njmax=!smalltest ? 128 : 64;
    const int nkmax=!smalltest ? 128 : 64;
    
    int nistart=1, njstart=1, nkstart=1;
    bool debugging = (nistart!=1 || njstart!=1 || nkstart!=1);
    
    for (int ni=nistart; ni<=nimax; ni+=1) {
      if (debugging) cout << "ni " << ni << endl;
      for (int nj=njstart; nj<=njmax; nj+=1) {
	if (debugging) cout << "  nj " << nj << endl;
	for (int nk=nkstart; nk<=nkmax; nk+=1) {
	  if (debugging) cout << "    nk " << nk << endl;
	  testone(mTxmq, ni, nj, nk, debugging, debugging);
        }
      }
    }
    printf("... OK!\n");
    
}

//#pragma GCC optimize "1"
template <typename funcT>
void timer(const char* s, funcT mTxmq, int ni, int nj, int nk, double *a, double *b, double *c) {
    double fastest=0.0, fastest_blas=0.0;

    double scale = 1e-9;
    double target = 1e-3;
    
    double nflop = 2.0*ni*nj*nk;
    double est = 1e-6 + nflop/DFLOPS;
    int ntrial=100;

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
    
    printf("%20s %3s %3s %3s %8s (GF/s)\n", "type", "M", "N", "K", "LOOP");
    // int itiles[] = {1,2,3,4,5,6,7,8};
    // for (int nj=2; nj<=12; nj+=2) {
    //     for (int ni : itiles) {
    //         int nk = 12;
    //         timer("(16,M)T*(16,N)", mTxmq, ni,nj,nk,a,b,c);
    //     }
    // }
    
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

template <typename funcT>
void search(funcT f, int ni, int nj, int nk) {
    const int len = 32768;
    
    double *a, *b, *c, *d;  
    if(posix_memalign((void **) &a, ALIGNMENT, len*sizeof(double))) exit(1);
    if(posix_memalign((void **) &b, ALIGNMENT, len*sizeof(double))) exit(1);
    if(posix_memalign((void **) &c, ALIGNMENT, len*sizeof(double))) exit(1);
    if(posix_memalign((void **) &d, ALIGNMENT, len*sizeof(double))) exit(1);

    double scale = 1e-9;
    double target = 1e-3;
    
    double nflop = 2.0*ni*nj*nk;
    double est = 1e-6 + nflop/2e10;
    int ntrial=100;

    // Ajust repetition count to hit target duration
    int nloop=std::max(1,int(target/est));
    for (int rep=0; rep<2; rep++) {
      PerfData perf;
      for (int loop=0; loop<nloop; ++loop) {
	f(ni,nj,nk,c,a,b,-1,-1);
      }
      perf.stop();
      double used = perf.cpu_time();
      nloop = std::max(1,int(nloop*target/used));
    }					

    struct stats {
        double rate;
        int it, jt;
        stats(double rate, int it, int jt) : rate(rate), it(it), jt(jt) {}
        bool operator<(const stats& other) {return rate<other.rate;}
    };
    std::list<stats> results;

    const int W=REGISTER_WIDTH;
    const int NR=NUMBER_OF_REGISTERS;
    const int maxjtile=MAX_JTILE;
    for (int jtile=std::min(W,nj); jtile<=std::min(maxjtile,nj); jtile++) {
        int jr = (jtile-1)/W+1;
        int maxitile = std::min(MAX_ITILE,std::min((NR-1)/jr-1 + 1,ni));
        for (int itile=std::max(1,maxitile-3); itile<=maxitile; itile++) {
            double fastest=0.0;
            for (int t=0; t<ntrial; t++) {
                double rate;
                PerfData perf;
                for (int loop=0; loop<nloop; ++loop) {
                    f(ni,nj,nk,c,a,b,itile,jtile);
                }
                perf.stop();
                double used = perf.cpu_time()/nloop;
                rate = scale*nflop/used;
                if (rate > fastest) fastest = rate;
            }
            results.push_back(stats(fastest, itile, jtile));
	    //std::cout << ni << " " << nj << " " << nk  << "     " << itile << " " <<  jtile << " " << fastest << std::endl;
        }
    }
    results.sort();
    printf("%4d %4d %4d", ni, nj, nk);
    size_t n = std::min(size_t(3),results.size());
    for (size_t i=0; i<n; i++) {
        stats s = results.back(); results.pop_back();
	printf("    %2d %2d %8.2f", s.it, s.jt, s.rate);
    }
    printf("\n");
    free(a);
    free(b);
    free(c);
    free(d);
}    
    

#include <sched.h>
#include <stdio.h>

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
    
    //tester(&mTxmqG);
    timer(&mTxmqG);

    auto num = std::chrono::high_resolution_clock::period::num;
    auto den = std::chrono::high_resolution_clock::period::den;
    std::cout << "period " << double(num)/den << std::endl;

    //testone(&mTxmqG,36,36,36,false, true);
    //tester(&mTxmqG);
 

    // std::cout << "square\n";
    // for (int n=1; n<129; n++) search(&mTxmqG,n,n,n);

    // std::cout << "square k=4\n";
    // for (int n=1; n<129; n++) search(&mTxmqG,n,n,4);

    // std::cout << "square k=8\n";
    // for (int n=1; n<65; n++) search(&mTxmqG,n,n,8);

    // std::cout << "small non-square k=8\n";
    // for (int n=1; n<33; n++) 
    //   for (int m=1; m<33; m++) 
    // 	search(&mTxmqG,n,m,8);

    std::cout << "(n,n*n)T(n,n)\n";
    for (int n=1; n<33; n++) search(&mTxmqG,n*n,n,n);

    // std::cout << "(n,n*n)T(n,m)\n";
    // for (int n=1; n<33; n++) {
    //   for (int m=1; m<=n; m++) {
    // 	search(&mTxmqG,n*n,m,n);
    //   }
    // }

    // std::cout << "(m,n*n)T(m,n)\n";
    // for (int n=1; n<33; n++) {
    //   for (int m=1; m<=n; m++) {
    // 	search(&mTxmqG,m,n*n,n);
    //   }
    // }

    
    return 0;
}
