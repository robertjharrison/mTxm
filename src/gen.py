from math import *

#from gen_avx2 import *
from gen_avx512 import *


#################################################

    
def intceil(x):
    return int(ceil(x))

tab="   ";

# Name of b register from j register index            
def b(j):
    return "bkj%2.2d"%j

# Name of c register from i tile index and j register index
def c(i,j):
    return "ci%2.2dj%2.2d"%(i,j)

# Name of kernel
def kernel_name(itile,jtile):
    return "kernel_i%2.2d_j%2.2d"%(itile,jtile)

def print_opening(itile,jtile):
    print("\nstatic void %s(const int dimi,const int dimj,const int dimk,%s* __restrict__ c,int incC,const %s* a,int incA,const %s* b,int incB) {" % (kernel_name(itile,jtile),DATA_TYPE,DATA_TYPE,DATA_TYPE))
    indent=tab
    print(indent,"const int itile=%d;"%itile)
    print(indent,"const int jtile=%d;"%jtile)
    print(indent,'if (kernel_trace) std::cout << "%s: dimi=" << dimi << " dimj=" << dimj << " dimk=" << dimk << " incA=" << incA << " incB=" << incB << " incC=" << incC << std::endl;' % kernel_name(itile,jtile))
    print(indent,"calls[itile][jtile]+=2.0*dimi*dimj*dimk;")
    
    print(indent,'__asm__("/*start %s*/");'%kernel_name(itile,jtile))

# declare a vector of registers over j index
def decl_jreg(indent,basename,njreg):
    print(indent,REGISTER_TYPE," ",end="")
    for j in range(njreg):
        print("%s%2.2d"%(basename,j),end="")
        if (j == (njreg-1)):
            print(";")
        else:
            print(",",end="")

def print_nounroll():
    print('''
#ifdef __clang__
#pragma clang loop unroll_count(1)
#elif defined(__INTEL_COMPILER)
// must be before GNUC since intel also defines this?
#pragma nounroll
#elif defined(__GNUC__)
#pragma GCC unroll 1
#else
#error "failed to detect compiler to stop loop unrolling"            
#endif''')
            
def gen(itile,jtile):
    njreg=intceil(jtile/REGISTER_WIDTH)

    # open function kernel_itile_jtile
    print_opening(itile,jtile)
    indent=tab
    decl_mask(indent,jtile);

    # loop over j tiles
    print(indent,"for (int j=0; j<dimj; j+=jtile,c+=jtile,b+=jtile) {")
    indent=2*tab
    
    # loop over i tiles
    print(indent,"%s* __restrict__ pcij=c;"%DATA_TYPE)
    print(indent,"for (int i=0; i<dimi; i+=itile,pcij+=itile*incC) {")
    indent=3*tab
    print(indent,"const %s* pbkj=b;"%DATA_TYPE)
    print(indent,"const %s* paki=a+i;"%DATA_TYPE)

    # declare registers
    print(indent,"%s aki;" % REGISTER_TYPE)
    decl_jreg(indent,"bkj",njreg);

    for i in range(itile):
        decl_jreg(indent,"ci%2.2dj"%i,njreg)

    ## # zero C registers --- eliminated by peeling first loop of k
    ## for i in range(itile):
    ##     print(indent+" ",end="")
    ##     for j in range(njreg):
    ##         zero(c(i,j))
    ##     print("")

    # peel off first loop of k replacing fma with mul
    print(indent+" ",end="")
    for j in range(njreg):
        is_incomplete = j==(njreg-1)and(jtile%REGISTER_WIDTH)        
        load(b(j),"pbkj+%d"%(j*REGISTER_WIDTH),is_incomplete);
    print("");
    for i in range(itile):
        broadcast(indent,"aki","paki+%2d"%i)
        for j in range(njreg):
            mul("aki",b(j),c(i,j))
        print("")
    print(indent,"pbkj+=incB,paki+=incA;");
    
    # loop over remainder of k
    #### for when we peel off front loop print(indent,"pbkj+=incB,paki+=incA;")
    print_nounroll()
    print(indent,"for (int k=1; k<dimk; k++,pbkj+=incB,paki+=incA) {")
    indent=4*tab

    # load b
    print(indent+" ",end="")
    for j in range(njreg):
        is_incomplete = j==(njreg-1)and(jtile%REGISTER_WIDTH)        
        load(b(j),"pbkj+%d"%(j*REGISTER_WIDTH),is_incomplete);
    print("");

    # finally do some computation
    for i in range(itile):
        broadcast(indent,"aki","paki+%2d"%i)
        for j in range(njreg):
            fma("aki",b(j),c(i,j))
        print("")
    indent=3*tab
    print(indent,"}")

    # store C
    for i in range(itile):
        print(indent+" ",end="")
        for j in range(njreg):
            is_incomplete = j==(njreg-1)and(jtile%REGISTER_WIDTH)        
            store(c(i,j),"pcij+incC*%2d+%d"%(i,j*REGISTER_WIDTH),is_incomplete)
        print("")

    indent=2*tab
    print(indent,"}")
    indent = tab
    print(indent,"}")
    print(indent,'__asm__("/*end %s*/");'%kernel_name(itile,jtile))
    print("}")

print_include()
print("static bool kernel_trace = false;")
print("void set_kernel_trace(bool value){kernel_trace=value;}")
print("using kernelT = void (*)(const int dimi,const int dimj,const int dimk,double* __restrict__ c,int incC,const double* a,int incA,const double* b,int incB);")
print("static const int MAX_ITILE=%d, MAX_JTILE=%d;" % (MAX_ITILE,MAX_JTILE))
print("static const int TARGET_JTILE=%d;"%TARGET_JTILE)
print("static const int REGISTER_WIDTH=%d;"%REGISTER_WIDTH)
print("static const int MASK_IN_REGISTER=%d;"%MASK_IN_REGISTER)
print("static const int NUMBER_OF_REGISTERS=%d;"%NUMBER_OF_REGISTERS)

print("static double calls[MAX_ITILE+1][MAX_JTILE+1] = {{0.0}};")
print("static kernelT dispatch[MAX_ITILE+1][MAX_JTILE+1] = {{nullptr}};")


# When generating the kernels permit use a of few more registers than actually
# exist since sometimes the compiler can do something creative like swtiching
# to use data directly from memory

# This loop nest must match the one immediately below
for jtile in range(1,MAX_JTILE+1):
    jr = (jtile-1)//REGISTER_WIDTH + 1
    maxitile = min((NUMBER_OF_REGISTERS-1)//jr-1 +2,MAX_ITILE)
    for itile in range(1,maxitile+1):
        gen(itile,jtile)

print("\nstatic void init_dispatch() {")
for jtile in range(1,MAX_JTILE+1):
    jr = (jtile-1)//REGISTER_WIDTH + 1
    maxitile = min((NUMBER_OF_REGISTERS-1)//jr-1 +2,MAX_ITILE)
    for itile in range(1,maxitile+1):
        print(tab,"dispatch[%2d][%2d]=&%s;"%(itile,jtile,kernel_name(itile,jtile)))
print("}")

print('''
static inline kernelT kernel(int itile, int jtile) {
  if (itile<0 || itile>MAX_ITILE) {std::cerr<<itile<<std::endl; throw "bad itile dispatching kernel";}
  if (jtile<0 || jtile>MAX_JTILE) {std::cerr<<itile<<std::endl; throw "bad jtile dispatching kernel";}
  kernelT f = dispatch[itile][jtile];
  if (!f) {std::cerr << itile << " " << jtile << std::endl; throw "kernel not found dispatching kernel";}
  return f;
}
''')

print_tuner()

print('''
void mTxmqG(int dimi, int dimj, int dimk, double* __restrict__ c, const double* a, const double* b, int itile=-1, int jtile=-1) {
  static bool initialized = false;
  if (!initialized) {
    init_dispatch();
    initialized=true;
  }

  const int incA = dimi;
  const int incB = dimj;
  const int incC = dimj;

  if (itile==-1 || jtile==-1) tune(dimi,dimj,dimk,itile,jtile);

  const int cachesize = 4096;
  int itile_cache = (cachesize-jtile*dimk-itile*jtile)/dimk;
  itile_cache = std::min(dimi,std::max(itile,itile_cache));
    
  if (itile_cache<dimi) {
    itile_cache -= (itile_cache%itile); // make it a multiple of target register tile size
  }

  int njtile = dimj/jtile; // no. of full j-tiles
  int nj0 = njtile*jtile;  // size of full j-block
  int nj1 = dimj - nj0;    // size of j-remainder

  // Loop thru cache blocks
  int ni = dimi;
  while (ni) {
    itile_cache = std::min(ni,itile_cache);
    itile = std::min(ni,itile);
    int nitile = itile_cache/itile;    // number of full (register) i-tiles per cache block
    int ni0 = nitile*itile;            // size of full i-block
    int ni1 = itile_cache-ni0;         // size of i-remainder
    
    (*kernel(std::min(ni0,itile),jtile))(ni0, nj0, dimk, c, incC, a, incA, b, incB);
    if (nj1) {
        (*kernel(std::min(ni0,itile),nj1))(ni0, nj1, dimk, c+nj0, incC, a, incA, b+nj0, incB);
    }
    if (ni1) {
        (*kernel(ni1,jtile))(ni1, nj0, dimk, c+ni0*incC, incC, a+ni0, incA, b, incB);
    }
    if (ni1 && nj1) {
        (*kernel(std::min(ni1,itile),nj1))(ni1, nj1, dimk, c+ni0*incC+nj0, incC, a+ni0, incA, b+nj0, incB);
    }

    a += itile_cache;
    c += itile_cache*incC;
    ni -= itile_cache;
  }
}
''')


f = open("mtxmq.h","w")
f.write("extern void set_kernel_trace(bool);\n")
f.write("extern void mTxmqG(int dimi, int dimj, int dimk, double* __restrict__ c, const double* a, const double* b, int itile=-1, int jtile=-1);\n")
f.write("static const int MAX_ITILE=%d, MAX_JTILE=%d;\n" % (MAX_ITILE,MAX_JTILE))
f.write("static const int TARGET_JTILE=%d;\n"%TARGET_JTILE)
f.write("static const int REGISTER_WIDTH=%d;\n"%REGISTER_WIDTH)
f.write("static const int MASK_IN_REGISTER=%d;\n"%MASK_IN_REGISTER)
f.write("static const int NUMBER_OF_REGISTERS=%d;\n"%NUMBER_OF_REGISTERS)
f.close()
