import os
import sys

tab = "  "

def parse(f,nlines):
    ''' For now just extract the first value '''
    lines = []
    for i in range(nlines): lines.append(f.readline())
    data = {}
    for line in lines:
        d = line.split()
        ni,nj,nk = int(d[0]),int(d[1]),int(d[2])
        data[(ni,nj,nk)] = (int(d[3]),int(d[4]),float(d[5]))
    return data

def print_list(indent,name,values):
    n = len(values)
    print(indent*tab,"static const unsigned char %s[%d] = {0,"%(name,n+1),end="")
    for i in range(n):
        if values[i]>255:
            print("ERROR --- value (%d) will not fit into char"%values[i],file=sys.stderr)
            sys.exit(1)
        print(values[i],end="")
        if i != n-1:
            print(",",end="")
    print("};")

def print_sq_tuner(indent,n,data):
    print(indent*tab,"if (ni==nj && nj==nk && ni<=%d) {"%n);
    indent += 1
    print_list(indent,"ilist",[data[i,i,i][0] for i in range(1,n+1)])
    print_list(indent,"jlist",[data[i,i,i][1] for i in range(1,n+1)])
    print(indent*tab,"itile=ilist[ni]; jtile=jlist[ni];")
    indent -= 1
    print(indent*tab,"}")

def print_mmm_tuner(indent,n,data):
    print(indent*tab,"else if (ni==nj*nj && nj==nk && nj<%d) {" % n);
    indent += 1
    print_list(indent,"ilist",[data[i*i,i,i][0] for i in range(1,n+1)])
    print_list(indent,"jlist",[data[i*i,i,i][1] for i in range(1,n+1)])
    print(indent*tab,"itile=ilist[nj]; jtile=jlist[nj];")
    indent -= 1
    print(indent*tab,"}")

def print_small_nonsq_tuner(indent,n,data):
    print(indent*tab,"else if (ni<=%d || nj<=%d) {"%(n,n));
    indent += 1
    print_list(indent,"ilist",[data[i,j,8][0] for i in range(1,n+1) for j in range(1,n+1)])
    print_list(indent,"jlist",[data[i,j,8][1] for i in range(1,n+1) for j in range(1,n+1)])
    print(indent*tab,"ni=std::min(ni,%d);"%n)
    print(indent*tab,"nj=std::min(nj,%d);"%n)
    print(indent*tab,"itile=ilist[ni+(nj-1)*%d]; jtile=jlist[ni+(nj-1)*%d];"%(n,n))
    indent -= 1
    print(indent*tab,"}")

def do_square(indent,f):
    line = f.readline()
    n = int(line.split()[1])
    data = parse(f,n)
    print_sq_tuner(indent,n,data)

def do_mmm(indent,f):
    line = f.readline()
    n = int(line.split()[1])
    data = parse(f,n)
    print_mmm_tuner(indent,n,data)

def do_small_nonsq(indent,f):
    line = f.readline()
    n = int(line.split()[3])
    data = parse(f,n*n)
    print_small_nonsq_tuner(indent,n,data)

def print_tuner(arch,fname):
    if not os.path.exists(fname):
        print("ERROR: tuner: file '%s' does exist"%fname, file=sys.stderr)
        sys.exit(1)

    try:
        f = open(fname,'r')
    except:
        print("tuner: error opening file '%s' for reading" % fname)
        
    print("void tune(int ni, int nj, int nk, int& itile, int& jtile) {")
    indent = 1
    do_square(indent,f)
    do_mmm(indent,f)
    do_small_nonsq(indent,f)
    print("%selse {\n%sitile=TARGET_ITILE;jtile=TARGET_JTILE;\n%s}" % (indent*tab,(indent+1)*tab,indent*tab))
    print('%sif (kernel_trace) printf("tuner: ni=%%d nj=%%d nk=%%d itile=%%d jtile=%%d\\n", ni, nj, nk, itile, jtile);'%tab)
    print("}")

