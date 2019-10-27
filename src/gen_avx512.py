##############################################
# This for AVX2

REGISTER_TYPE="__m512d"
DATA_TYPE="double"
REGISTER_WIDTH = 8
NUMBER_OF_REGISTERS = 32
MASK_IN_REGISTER = 0
MAX_JTILE = 32
MAX_ITILE = 30
TARGET_JTILE = 24

def zero(register):
    print("%s=_mm512_setzero_pd(); " % register,end="")

def load(register,ptr,is_incomplete):
    if is_incomplete:
        print("%s=_mm512_mask_loadu_pd(%s,mask,%s); " % (register,register,ptr),end="")
    else:
        print("%s=_mm512_loadu_pd(%s); " % (register,ptr),end="")

def store(register,ptr,is_incomplete):
    if is_incomplete:
        print("_mm512_mask_storeu_pd(%s, mask, %s); " % (ptr,register),end="")
    else:
        print("_mm512_storeu_pd(%s, %s); " % (ptr,register),end="")
        
# if needed declare and initialize the mask for incomplete j-register operation
def decl_mask(indent,jtile):
    rem=jtile%REGISTER_WIDTH
    if rem:
        print("const unsigned char mask=%d;" % ((1<<rem)-1))

def broadcast(indent,register,ptr):
    print(indent,"%s=_mm512_broadcastsd_pd(_mm_load_sd(%s));; "%(register,ptr),end="") # optimized to vbroadcastsd (%r8,%r11,8), %zmm13

def fma(a,b,c):
    print("%s=_mm512_fmadd_pd(%s,%s,%s); "%(c,a,b,c),end="")

def mul(a,b,c):
    print("%s=_mm512_mul_pd(%s,%s); "%(c,a,b),end="")

def print_include():
    print("#include <immintrin.h>")
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
