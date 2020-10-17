##############################################
# This for SVE with 512bit register (not VLA)

REGISTER_TYPE="vec"
DATA_TYPE="double"
REGISTER_WIDTH = 8
NUMBER_OF_REGISTERS = 32
MASK_IN_REGISTER = 0
MAX_JTILE = 48
MAX_ITILE = 30
TARGET_ITILE = 4   
TARGET_JTILE = 48
# 6,32 or I=4,J=48 also OK ... symmetric since 32/8=4 6*8=48 

def zero(register):
    print("%s=svdup_n_f64(0.0); " % register,end="")

def load(register,ptr,is_incomplete):
    if is_incomplete:
        print("%s=svld1_f64(mask,%s); " % (register,ptr),end="") 
    else:
        print("%s=svld1_f64(everything,%s); " % (register,ptr),end="")

def store(register,ptr,is_incomplete):
    if is_incomplete:    
        print("svstnt1_f64(mask, %s, %s); " % (ptr,register),end="")
    else:
        print("svstnt1_f64(everything, %s, %s); " % (ptr,register),end="")
        
# if needed declare and initialize the mask for incomplete j-register operation
def decl_mask(indent,jtile):
    rem=jtile%REGISTER_WIDTH
    print(indent,"pred everything=svptrue_b64();")
    if rem:
        print(indent,"pred mask=svptrue_pat_b64(svpattern(%s));" % rem)

def broadcast(indent,register,ptr):
    print(indent,"%s=svdup_n_f64(*(%s)); "%(register,ptr),end="") # optimized to vbroadcastsd (%r8,%r11,8), %zmm13

def fma(a,b,c):
    print("%s=svmad_f64_x(everything,%s,%s,%s); "%(c,a,b,c),end="")

def mul(a,b,c):
    print("%s=svmul_f64_x(everything,%s,%s); "%(c,a,b),end="")

def print_include():
    print("#include <arm_sve.h>")
    print("#include <iostream>")
    print("#include <cassert>")
    print("typedef svfloat64_t vec __attribute__((arm_sve_vector_bits(512)));")
    print("typedef svbool_t pred __attribute__((arm_sve_vector_bits(512)));")

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
    if (kernel_trace) printf("tuner_sve512: general: ni=%d nj=%d nk=%d itile=%d jtile=%d\\n", ni, nj, nk, itile, jtile);
}
''')
