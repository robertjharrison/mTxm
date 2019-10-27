##############################################
# This for AVX2

REGISTER_TYPE="float64x2_t"
DATA_TYPE="double"
REGISTER_WIDTH = 2
NUMBER_OF_REGISTERS = 16
MASK_IN_REGISTER = 0
MAX_JTILE = 14
MAX_ITILE = 14
TARGET_JTILE = 8

def zero(register):
    print("%s=vsub_f64(%s,%s); " % (register,register,register),end="")

def load(register,ptr,is_incomplete):
    if is_incomplete:
        print("%s=vld1q_dup_f64(%s); " % (register,ptr),end="")
    else:
        print("%s=vld1q_f64(%s); " % (register,ptr),end="")

def store(register,ptr,is_incomplete):
    if is_incomplete:
        print("vst1q_lane_f64(%s,%s,0); " % (ptr,register),end="")
    else:
        print("vst1q_f64(%,%s); " % (ptr,register),end="")
        
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

def print_tuner():
    print('''
void tune(int ni, int nj, int nk, int& itile, int& jtile) {
    if (nj<=MAX_JTILE) {
       jtile=nj;
    }
    else {
       jtile=TARGET_JTILE;
    }
    int njreg = (jtile-1)/REGISTER_WIDTH + 1;
    itile = std::min((NUMBER_OF_REGISTERS-1)/njreg-1,ni);
    itile = std::min(MAX_ITILE,itile);
    if (kernel_trace) printf("tuner_avx512: general: ni=%d nj=%d nk=%d itile=%d jtile=%d\\n", ni, nj, nk, itile, jtile);
}
''')
