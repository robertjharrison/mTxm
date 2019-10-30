##############################################
# This for AVX2

REGISTER_TYPE="__m256d"
DATA_TYPE="double"
REGISTER_WIDTH = 4
NUMBER_OF_REGISTERS = 16
MASK_IN_REGISTER = 0 # should be 1 ??
TARGET_JTILE = 12
TARGET_ITILE = 4
MAX_JTILE = 20
MAX_ITILE = 14 + 2

def zero(register):
    print("%s=_mm256_setzero_pd(); " % register,end="")

def load(register,ptr,is_incomplete):
    if is_incomplete:
        print("%s=_mm256_maskload_pd(%s,mask); " % (register,ptr),end="")
    else:
        print("%s=_mm256_loadu_pd(%s); " % (register,ptr),end="")

def store(register,ptr,is_incomplete):
    if is_incomplete:
        print("_mm256_maskstore_pd(%s, mask, %s); " % (ptr,register),end="")
    else:
        print("_mm256_storeu_pd(%s, %s); " % (ptr,register),end="")
        
# if needed declare and initialize the mask for incomplete j-register operation
def decl_mask(indent,jtile):
    rem=jtile%REGISTER_WIDTH
    if rem==0:
        pass
    elif rem==1:
        print(indent,"__m256i mask=_mm256_set_epi64x( 0,0,0,-1);")
    elif rem==2:
        print(indent,"__m256i mask=_mm256_set_epi64x( 0,0,-1,-1);")
    else:
        print(indent,"__m256i mask=_mm256_set_epi64x( 0,-1,-1,-1);")

def broadcast(indent,register,ptr):
    print(indent,"%s=_mm256_broadcast_sd(%s); "%(register,ptr),end="")

def fma(a,b,c):
    print("%s=_mm256_fmadd_pd(%s,%s,%s); "%(c,a,b,c),end="")

def mul(a,b,c):
    print("%s=_mm256_mul_pd(%s,%s); "%(c,a,b),end="")

def print_include():
    print("#include <immintrin.h>")
    print("#include <iostream>")
    print("#include <cassert>")

def print_tuner():
    print('''
void tune(int ni, int nj, int nk, int& itile, int& jtile) {
  if (ni==nj && nj==nk) { // SQUARE
    if (nj>20) {
      if (((nj&7)&&!(nj&3)) || !(nj%12)) {itile=4; jtile=12;}
      else {itile=6; jtile=8;}
    }
    else if (nj<=7){itile=jtile=nj;}
    else {
      static const int ilist[] = {4, 3, 4, 4, 4, 4, 5, 5, 6, 4, 6, 4, 4};
      static const int jlist[] = {8, 9,10,11,12,13, 7, 8, 8,12, 8,12,12};
      itile=ilist[nj-8];
      jtile=jlist[nj-8];
    }
    if (kernel_trace) printf("tuner_avx2: square: ni=%d nj=%d nk=%d itile=%d jtile=%d\\n", ni, nj, nk, itile, jtile);
  }
  else if (ni==nj*nj && nj==nk) { // (m,m*m)T*(m,m)
    if (nj>20) {itile=6; jtile=8;}
    else {
      static const int ilist[] = {-1, 1, 2, 9,12, 5, 6, 5, 4, 4, 4, 4, 4, 6, 6, 5, 6, 4, 4, 4, 4};
      static const int jlist[] = {-1, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 8, 8, 8, 8,12,12,12,12};
      itile=ilist[nj];
      jtile=jlist[nj];
    }
    if (kernel_trace) printf("tuner_avx2: (m,m*m)T(m,m): ni=%d nj=%d nk=%d itile=%d jtile=%d\\n", ni, nj, nk, itile, jtile);
  }
  else { // GENERAL
    /*
      Not quite doing this yet but it explains most of the data
      
      0) maximize occupancy and efficiency
           njreg  maxitile   jtile      nloadins  nflopins   eff
             1     14        1-4          it+1      it       jtile/4
             2     6(7)      5-8         2*(it+1)  2*it      jtile/8
             3     4(5)      9-12        3*(it+1)  3*it      jtile/12
             4     2(3)     12-16        ...                 jtile/16
             5     2        17-20        ...                 jtile/20

             and need at least 4 cycles per iter

             So from (large) itile perspective all are equal efficiency, so just
             look at joccupancy to compute efficiency ... also applies to store

             njreg=(jtile-1)//REGISTER_WIDTH+1 or (NUMBER_OF_REGISTERS-1)//(itile+1)
             itile=(NUMBER_OF_REGISTERS-1-njreg)//njreg

      1) for same occ/eff prefer to make remj=0 

      2) otherwise try make remi=0

      3) 1 extra itile when () above due to GCC optimization
    */

    if (nj<=20) {      // small nj any ni
      jtile = nj;
      if (ni>4) {jtile=std::min(8,jtile);}
      int njreg = (jtile-1)/REGISTER_WIDTH + 1;
      itile = std::min((NUMBER_OF_REGISTERS-1-njreg)/njreg,ni);
    }
    else if (ni<=3) {    // small ni large nj
      itile=ni;
      jtile=20;
    }
    else if (ni<=5) {    // small ni large nj
      itile=ni;
      jtile=12;
    }
    else if (ni<=7) {    // small ni large nj
      itile=ni;
      jtile=8;
    }
    else {          // large both
      itile=6;
      jtile=8;
    }
    if (kernel_trace) printf("tuner_avx2: general: ni=%d nj=%d nk=%d itile=%d jtile=%d\\n", ni, nj, nk, itile, jtile);
  }
}
''')
