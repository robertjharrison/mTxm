import sys
import os
import math
import importlib

import tune

def print_help():
    print("""
Usage: python3 gen.py target [tuning data file name]")

    target.py     Name of the python module defining the target architecture
                  (i.e., the name of a python script omitting the trailing '.py').
                  Example targets include: 'neon', 'avx2', 'avx512', 'sve'.

    tuning data   Optionally provide the name of a file containing tuning
                  data generated by running the tuner executable
                  (e.g., 'tune_graviton.txt', 'tune_haswell.txt', 'tune_skylake.txt')
""")

if len(sys.argv) == 2:
    target = sys.argv[1]
    tuningdata = "none"
    #print("Generating for target '%s' with default tuner"%target)
elif len(sys.argv) == 3:
    target = sys.argv[1]
    tuningdata = sys.argv[2]
    #print("Generating for target '%s' with tuning data from '%s'"%(target,tuningdata))
else:
    print_help()
    sys.exit(1)

arch =  importlib.import_module(target)


#################################################

def intceil(x):
    ''' Ceiling of floating point value as an integer '''
    return int(math.ceil(x))

tab="   ";

def b(j):
    '''Name of b register from j register index '''
    return "bkj%2.2d"%j

def c(i,j):
    ''' Name of c register from i tile index and j register index'''
    return "ci%2.2dj%2.2d"%(i,j)

def kernel_name(itile,jtile):
    ''' Name of kernel from i/j tile sizes '''
    return "kernel_i%2.2d_j%2.2d"%(itile,jtile)

def print_opening(itile,jtile):
    ''' Print the stuff at the beginning of a kernel ''' 
    print("\nstatic void %s(const int dimi,const int dimj,const int dimk,%s* __restrict__ c,int incC,const %s* a,int incA,const %s* b,int incB) {" % (kernel_name(itile,jtile),arch.DATA_TYPE,arch.DATA_TYPE,arch.DATA_TYPE))
    indent=tab
    print(indent,"const int itile=%d;"%itile)
    print(indent,"const int jtile=%d;"%jtile)
    print(indent,'if (kernel_trace) std::cout << "%s: dimi=" << dimi << " dimj=" << dimj << " dimk=" << dimk << " incA=" << incA << " incB=" << incB << " incC=" << incC << std::endl;' % kernel_name(itile,jtile))
    print(indent,"calls[itile][jtile]+=2.0*dimi*dimj*dimk;")
    
def decl_jreg(indent,basename,njreg):
    ''' declare a vector of registers over j index '''
    print(indent,arch.REGISTER_TYPE," ",end="")
    for j in range(njreg):
        print("%s%2.2d"%(basename,j),end="")
        if (j == (njreg-1)):
            print(";")
        else:
            print(",",end="")


def print_nounroll():
    ''' Print pragma to prevent loop unrolling that is usually detrimental to performance '''
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
    ''' Generate a kernel for given itile and jtile sizes'''
    njreg=intceil(jtile/arch.REGISTER_WIDTH)

    # open function kernel_itile_jtile
    print_opening(itile,jtile)
    indent=tab
    arch.decl_mask(indent,jtile);

    # loop over j tiles
    print(indent,"for (int j=0; j<dimj; j+=jtile,c+=jtile,b+=jtile) {")
    indent=2*tab
    
    # loop over i tiles
    print(indent,"%s* __restrict__ pcij=c;"%arch.DATA_TYPE)
    print(indent,"for (int i=0; i<dimi; i+=itile,pcij+=itile*incC) {")
    indent=3*tab
    print(indent,"const %s* pbkj=b;"%arch.DATA_TYPE)
    print(indent,"const %s* paki=a+i;"%arch.DATA_TYPE)

    # declare registers
    print(indent,"%s aki;" % arch.REGISTER_TYPE)
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
        is_incomplete = j==(njreg-1)and(jtile%arch.REGISTER_WIDTH)        
        arch.load(b(j),"pbkj+%d"%(j*arch.REGISTER_WIDTH),is_incomplete);  ############# HERE don't +j if j==0
    print("");
    for i in range(itile):
        arch.broadcast(indent,"aki","paki+%2d"%i)                         ############# HERE don't +i if i==0
        for j in range(njreg):
            arch.mul("aki",b(j),c(i,j))
        print("")
    print(indent,"pbkj+=incB,paki+=incA;");
    
    print(indent,'__asm__("/*start %s inner loop*/");'%kernel_name(itile,jtile))
    # loop over remainder of k
    #### for when we peel off front loop print(indent,"pbkj+=incB,paki+=incA;")
    print_nounroll()
    print(indent,"for (int k=1; k<dimk; k++,pbkj+=incB,paki+=incA) {")
    indent=4*tab

    # load b
    print(indent+" ",end="")
    for j in range(njreg):
        is_incomplete = j==(njreg-1)and(jtile%arch.REGISTER_WIDTH)        
        arch.load(b(j),"pbkj+%d"%(j*arch.REGISTER_WIDTH),is_incomplete);
    print("");

    # finally do some computation
    for i in range(itile):
        arch.broadcast(indent,"aki","paki+%2d"%i)
        for j in range(njreg):
            arch.fma("aki",b(j),c(i,j))
        print("")
    indent=3*tab
    print(indent,"}")
    print(indent,'__asm__("/*end %s*/");'%kernel_name(itile,jtile))
    

    # store C
    for i in range(itile):
        print(indent+" ",end="")
        for j in range(njreg):
            is_incomplete = j==(njreg-1)and(jtile%arch.REGISTER_WIDTH)        
            arch.store(c(i,j),"pcij+incC*%2d+%d"%(i,j*arch.REGISTER_WIDTH),is_incomplete)
        print("")

    indent=2*tab
    print(indent,"}")
    indent = tab
    print(indent,"}")
    print("}")

arch.print_include()
print("static bool kernel_trace = false;")
print("void set_kernel_trace(bool value){kernel_trace=value;}")
print("using kernelT = void (*)(const int dimi,const int dimj,const int dimk,double* __restrict__ c,int incC,const double* a,int incA,const double* b,int incB);")
print("static const int MAX_ITILE=%d, MAX_JTILE=%d;" % (arch.MAX_ITILE,arch.MAX_JTILE))
print("static const int TARGET_ITILE=%d;"%arch.TARGET_ITILE)
print("static const int TARGET_JTILE=%d;"%arch.TARGET_JTILE)
print("static const int REGISTER_WIDTH=%d;"%arch.REGISTER_WIDTH)
print("static const int MASK_IN_REGISTER=%d;"%arch.MASK_IN_REGISTER)
print("static const int NUMBER_OF_REGISTERS=%d;"%arch.NUMBER_OF_REGISTERS)

print("static double calls[MAX_ITILE+1][MAX_JTILE+1] = {{0.0}};")
print("static kernelT dispatch[MAX_ITILE+1][MAX_JTILE+1] = {{nullptr}};")


# When generating the kernels permit use a of few more registers than actually
# exist since sometimes the compiler can do something creative like swtiching
# to use data directly from cache ... must match loop in search

# This loop nest must match the one immediately below
for jtile in range(1,arch.MAX_JTILE+1):
    jr = (jtile-1)//arch.REGISTER_WIDTH + 1
    maxitile = min((arch.NUMBER_OF_REGISTERS-1)//jr-1 +2,arch.MAX_ITILE)
    for itile in range(1,maxitile+1):
        gen(itile,jtile)

print("\nstatic void init_dispatch() {")
for jtile in range(1,arch.MAX_JTILE+1):
    jr = (jtile-1)//arch.REGISTER_WIDTH + 1
    maxitile = min((arch.NUMBER_OF_REGISTERS-1)//jr-1 +2,arch.MAX_ITILE)
    for itile in range(1,maxitile+1):
        print(tab,"dispatch[%2d][%2d]=&%s;"%(itile,jtile,kernel_name(itile,jtile)))
print("}")

print('''
static inline kernelT kernel(int itile, int jtile) {
  if (itile<0 || itile>MAX_ITILE) {std::cerr<<itile<<std::endl; throw "bad itile dispatching kernel";}
  if (jtile<0 || jtile>MAX_JTILE) {std::cerr<<jtile<<std::endl; throw "bad jtile dispatching kernel";}
  kernelT f = dispatch[itile][jtile];
  if (!f) {std::cerr << itile << " " << jtile << std::endl; throw "kernel not found dispatching kernel";}
  return f;
}
''')


if tuningdata == 'none':
    arch.print_tuner()
else:
    tune.print_tuner(arch,tuningdata)

print('''
void mTxmqG(int dimi, int dimj, int dimk, double* __restrict__ c, int incC, const double* a, int incA, const double* b, int incB, int itile=-1, int jtile=-1) {
  static bool initialized = false;
  if (!initialized) {
    init_dispatch();
    initialized=true;
  }

  //const int incA = dimi;
  //const int incB = dimj;
  //const int incC = dimj;

  if (itile==-1 || jtile==-1) tune(dimi,dimj,dimk,itile,jtile);

  const int cachesize = %(CACHE_SIZE)s;
  int itile_cache = (cachesize-jtile*dimk-itile*jtile)/dimk;
  itile_cache = std::min(dimi,std::max(itile,itile_cache));
    
  if (itile_cache<dimi) {
    itile_cache -= (itile_cache%%itile); // make it a multiple of target register tile size
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
''' % {'CACHE_SIZE':arch.CACHE_SIZE})


f = open("mTxmq.h","w")
f.write("extern void set_kernel_trace(bool);\n")
f.write("extern void mTxmqG(int dimi, int dimj, int dimk, double* __restrict__ c, int incC, const double* a, int incA, const double* b, int incB, int itile=-1, int jtile=-1);\n")
f.write("namespace mTxmqdetail {\n")
f.write("  static const int MAX_ITILE=%d, MAX_JTILE=%d;\n" % (arch.MAX_ITILE,arch.MAX_JTILE))
f.write("  static const int TARGET_ITILE=%d;\n"%arch.TARGET_ITILE)
f.write("  static const int TARGET_JTILE=%d;\n"%arch.TARGET_JTILE)
f.write("  static const int REGISTER_WIDTH=%d;\n"%arch.REGISTER_WIDTH)
f.write("  static const int MASK_IN_REGISTER=%d;\n"%arch.MASK_IN_REGISTER)
f.write("  static const int NUMBER_OF_REGISTERS=%d;\n"%arch.NUMBER_OF_REGISTERS)
f.write("}\n")
f.close()
