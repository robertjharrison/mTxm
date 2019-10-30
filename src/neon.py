##############################################
# This for AVX2

REGISTER_TYPE="float64x2_t"
DATA_TYPE="double"
REGISTER_WIDTH = 2
NUMBER_OF_REGISTERS = 32
MASK_IN_REGISTER = 0
MAX_ITILE = 20
MAX_JTILE = 14
TARGET_ITILE = 14
TARGET_JTILE = 4

def zero(register):
    print("%s=vld1q_dup_f64(&dzero); " % (register),end="")

def load(register,ptr,is_incomplete):
    if is_incomplete:
        print("%s=vld1q_dup_f64(%s); " % (register,ptr),end="")
    else:
        print("%s=vld1q_f64(%s); " % (register,ptr),end="")

def store(register,ptr,is_incomplete):
    if is_incomplete:
        print("vst1q_lane_f64(%s,%s,0); " % (ptr,register),end="")
    else:
        print("vst1q_f64(%s,%s); " % (ptr,register),end="")
        
# if needed declare and initialize the mask for incomplete j-register operation
def decl_mask(indent,jtile):
    pass

def broadcast(indent,register,ptr):
    print(indent,"%s=vld1q_dup_f64(%s); "%(register,ptr),end="")

def fma(a,b,c):
    print("%s=vfmaq_f64(%s,%s,%s); "%(c,c,a,b),end="")

def mul(a,b,c):
    print("%s=vmulq_f64(%s,%s); "%(c,a,b),end="")

def print_include():
    print("#include <arm_neon.h>")
    print("#include <iostream>")
    print("#include <cassert>")
    print("static const double dzero=0.0;")

def print_tuner():
    print('''
void tune(int ni, int nj, int nk, int& itile, int& jtile) {
    if (ni==nj) { // square
       if (ni<=7) {itile=jtile=ni;}
       else if (ni<14) {
          static const int ilist[] = { 8, 9, 5,11,12,13};
          static const int jlist[] = { 4, 6,10, 4, 4, 4};
          itile=ilist[ni-8];
          jtile=jlist[ni-8];
       }
       else {itile=14; jtile=4;}
       if (kernel_trace) printf("tuner_neon: square: ni=%d nj=%d nk=%d itile=%d jtile=%d\\n", ni, nj, nk, itile, jtile);
    }
    else if (ni==nj*nj) { // 
       if (nj<=7) {
          static const int ilist[] = {-1,1,4,9,8,5,9,5};
          itile = ilist[nj];
          jtile = nj;
       }
       else {itile=14; jtile=4;}
       if (kernel_trace) printf("tuner_neon: (m,m*m)T(m,m): ni=%d nj=%d nk=%d itile=%d jtile=%d\\n", ni, nj, nk, itile, jtile);
    }
    else {
       if (ni<=14) {
          itile = ni;
          jtile = std::min(MAX_JTILE,(NUMBER_OF_REGISTERS-1)/(itile+1));
          jtile = std::min(nj,jtile);
       }
       else {
          itile = 14;
          jtile = std::min(nj,4);
       }
       if (kernel_trace) printf("tuner_neon: general: ni=%d nj=%d nk=%d itile=%d jtile=%d\\n", ni, nj, nk, itile, jtile);
    }
}
''')
