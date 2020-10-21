#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <list>

#include "mTxmq.h"

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

int main() {


    try {
        //testone(&mTxmqG,400,1,20,false, true);
      tester(&mTxmqG);
    }
    catch (const char* c) {
      std::cout << "exception " << c << std::endl;
      return 1;
    }
    
    return 0;
}
