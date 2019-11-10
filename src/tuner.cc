#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <list>
#include <sched.h>

#include "mTxmq.h"

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

void search(int ni, int nj, int nk) {
    const int len = 32768*2;
    
    double *a, *b, *c;  
    if(posix_memalign((void **) &a, ALIGNMENT, len*sizeof(double))) exit(1);
    if(posix_memalign((void **) &b, ALIGNMENT, len*sizeof(double))) exit(1);
    if(posix_memalign((void **) &c, ALIGNMENT, len*sizeof(double))) exit(1);
    ran_fill(len,a);
    ran_fill(len,b);
    ran_fill(len,c);

    double target = 1e-3; // target time for kernel in seconds
    
    double nflop = 2.0*ni*nj*nk;
    double est = 1e-6 + nflop/2e10;
    int ntrial=10;

    // Ajust repetition count to hit target duration
    int nloop=std::max(1,int(target/est));
    for (int rep=0; rep<2; rep++) {
      PerfData perf;
      for (int loop=0; loop<nloop; ++loop) {
	mTxmqG(ni,nj,nk,c,a,b,-1,-1);
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
    for (int jtile=std::min(W,nj); jtile<=std::min(MAX_JTILE,nj); jtile++) {
        int jr = (jtile-1)/W+1;
	// this must match corresponding loop in gen.py
        int maxitile = std::min(MAX_ITILE,std::min((NR-1)/jr-1 + 2,ni));
        //for (int itile=std::max(1,maxitile-3); itile<=maxitile; itile++) {
	for (int itile=1; itile<=maxitile; itile++) {
	  //std::cout << "search " << itile << " " << jtile << std::endl;
            double fastest=0.0;
            for (int t=0; t<ntrial; t++) {
                double rate;
                PerfData perf;
                for (int loop=0; loop<nloop; ++loop) {
                    mTxmqG(ni,nj,nk,c,a,b,itile,jtile);
                }
                perf.stop();
                double used = perf.cpu_time()/nloop;
                rate = 1e-9*nflop/used;
                if (rate > fastest) fastest = rate;
            }
            results.push_back(stats(fastest, itile, jtile));
	    //std::cout << ni << " " << nj << " " << nk  << "     " << itile << " " <<  jtile << " " << fastest << std::endl;
        }
    }
    results.sort();
    printf("%4d %4d %4d", ni, nj, nk);
    double best = results.back().rate;
    size_t n = results.size();
    for (size_t i=0; i<n; i++) {
        stats s = results.back(); results.pop_back();
	if (s.rate>0.95*best)
	  printf("      %2d %2d %8.2f", s.it, s.jt, s.rate);
    }
    printf("\n");  fflush(stdout);
    free(a);
    free(b);
    free(c);
}    
    


int main() {

    // Attempt get more consisent timing by binding thread to a core
    cpu_set_t mask;
    CPU_ZERO(&mask);
    CPU_SET(1,&mask);
    if (sched_setaffinity(0, sizeof(mask), &mask) == -1) {
        perror("system error message");
        std::cout << "ThreadBase: set_affinity: Could not set cpu affinity" << std::endl;
    }

    try {
      const int maxsq = 64;
      std::cout << "square " << maxsq << std::endl;
      for (int n=1; n<=maxsq; n++) search(n,n,n);
      
      const int maxmmm = 32;
      std::cout << "(n,n*n)T(n,n) " << maxmmm << std::endl;
      for (int n=1; n<=maxmmm; n++) search(n*n,n,n);

      const int maxnsq = 32;
      std::cout << "small non-square k=8 " << maxnsq << std::endl;
      for (int n=1; n<=maxnsq; n++) 
        for (int m=1; m<=maxnsq; m++) 
	  search(n,m,8);

      std::cout << "(16,16*16)T(16,n) " << maxmmm << std::endl;
      for (int n=1; n<=maxmmm; n++) search(16*16,n,16);

      std::cout << "(16,n)T(16,16*16) " << maxmmm << std::endl;
      for (int n=1; n<=maxmmm; n++) search(n,16*16,16);

      // std::cout << "(20,20*20)T(20,n)\n";
      //  for (int n=1; n<21; n++) search(400,n,20);

       // std::cout << "square k=4\n";
       // for (int n=1; n<129; n++) search(n,n,4);

       // std::cout << "square k=8\n";
       // for (int n=1; n<65; n++) search(n,n,8);

       // std::cout << "(n,n*n)T(n,m)\n";
       // for (int n=1; n<33; n++) {
       //   for (int m=1; m<=n; m++) {
       // 	search(n*n,m,n);
       //   }
       // }

       // std::cout << "(m,n*n)T(m,n)\n";
       // for (int n=1; n<33; n++) {
       //   for (int m=1; m<=n; m++) {
       // 	search(m,n*n,n);
       //   }
       // }
    }
    catch (const char* c) {
      std::cout << "exception " << c << std::endl;
      return 1;
    }
    
    return 0;
}
