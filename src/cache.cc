#include <iostream>
#include "timerstuff.h"

/***

    Try to estimate the associativity of L1 cache by timing how long
    it takes to read from elements separated by the cachesize.  There
    should be a sharp increase in time past the level of
    associativity.

    Intel SkyLake i7-7600U --- documented as 8-way

    1   2.60014 1.07687
    2   2.68136 0.522122
    3   1.82105 0.512525
    4   1.36736 0.511935
    5   1.08757 0.514911
    6   0.911601 0.51192
    7   0.783775 0.510351
    8   0.685113 0.510865
    9   0.462292 0.672975  <----
    10  0.343269 0.815686
    11  0.271632 1.03081
    12  0.223672 1.25184

    

 ***/

int main() {
    const int CACHESIZE = 32768/sizeof(int);
    volatile int c[CACHESIZE*12] = {0};
    const int nloop = 200000000;
    const double freq = 2.8e9;
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "1   " << rate << " " << cyc/1 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "2   " << rate << " " << cyc/2 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "3   " << rate << " " << cyc/3 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "4   " << rate << " " << cyc/4 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
            sum += c[CACHESIZE*4];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "5   " << rate << " " << cyc/5 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
            sum += c[CACHESIZE*4];
            sum += c[CACHESIZE*5];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "6   " << rate << " " << cyc/6 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
            sum += c[CACHESIZE*4];
            sum += c[CACHESIZE*5];
            sum += c[CACHESIZE*6];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "7   " << rate << " " << cyc/7 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
            sum += c[CACHESIZE*4];
            sum += c[CACHESIZE*5];
            sum += c[CACHESIZE*6];
            sum += c[CACHESIZE*7];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "8   " << rate << " " << cyc/8 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
            sum += c[CACHESIZE*4];
            sum += c[CACHESIZE*5];
            sum += c[CACHESIZE*6];
            sum += c[CACHESIZE*7];
            sum += c[CACHESIZE*8];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "9   " << rate << " " << cyc/9 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
            sum += c[CACHESIZE*4];
            sum += c[CACHESIZE*5];
            sum += c[CACHESIZE*6];
            sum += c[CACHESIZE*7];
            sum += c[CACHESIZE*8];
            sum += c[CACHESIZE*9];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "10  " << rate << " " << cyc/10 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
            sum += c[CACHESIZE*4];
            sum += c[CACHESIZE*5];
            sum += c[CACHESIZE*6];
            sum += c[CACHESIZE*7];
            sum += c[CACHESIZE*8];
            sum += c[CACHESIZE*9];
            sum += c[CACHESIZE*10];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "11  " << rate << " " << cyc/10 << std::endl;
    }
    {
        int sum = 0;
        PerfData perf;
        for (int i=0; i<nloop; i++) {
            sum += c[0];
            sum += c[CACHESIZE];
            sum += c[CACHESIZE*2];
            sum += c[CACHESIZE*3];
            sum += c[CACHESIZE*4];
            sum += c[CACHESIZE*5];
            sum += c[CACHESIZE*6];
            sum += c[CACHESIZE*7];
            sum += c[CACHESIZE*8];
            sum += c[CACHESIZE*9];
            sum += c[CACHESIZE*10];
            sum += c[CACHESIZE*11];
        }
        perf.stop();
        double used = perf.cpu_time();
        double rate = 1e-9*nloop/used;
        double cyc = freq*used/nloop;
        std::cout << "12  " << rate << " " << cyc/10 << std::endl;
    }
    return 0;
}
