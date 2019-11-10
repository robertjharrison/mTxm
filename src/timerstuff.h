#ifndef TIMER_STUFF_H
#define TIMER_STUFF_H

#include <stdlib.h>
#include <stdint.h>
#include <ctime>
#include <ratio>
#include <chrono>

#ifdef __x86_64
/////// Cycle count on x86
// Returns current cycle count for this thread
static inline uint64_t cycle_count() {
    uint64_t x;
    unsigned int a,d;
    __asm__ volatile("rdtsc" : "=a"(a), "=d"(d));
    x = ((uint64_t)a) | (((uint64_t)d)<<32);
    return x;
}
#else
static inline uint64_t cycle_count() {
  return 0;
}
#endif

/////////// Uses protected system counter to access cycle count on ARM
#include <unistd.h>
#include <string.h>
#include <stdio.h>
#include <sys/ioctl.h>
#include <linux/perf_event.h>
#include <asm/unistd.h>

static long perf_event_open(struct perf_event_attr *hw_event, pid_t pid,
			    int cpu, int group_fd, unsigned long flags)
{
  return syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
}

static int fd;

static void init_perf_ctr() {
  struct perf_event_attr pe;
  memset(&pe, 0, sizeof(struct perf_event_attr));
  pe.type = PERF_TYPE_HARDWARE;
  pe.size = sizeof(struct perf_event_attr);
  pe.config = PERF_COUNT_HW_CPU_CYCLES; //PERF_COUNT_HW_INSTRUCTIONS;
  pe.disabled = 1;
  pe.exclude_kernel = 1;
  pe.exclude_hv = 1;
  
  fd = perf_event_open(&pe, 0, -1, -1, 0);
  if (fd == -1) {
    fprintf(stderr, "Error opening leader %llx\n", pe.config);
    exit(1);
  }
}

static void term_perf_ctr() {
  close(fd);
}

static inline void start_perf_ctr() {
  ioctl(fd, PERF_EVENT_IOC_RESET, 0);
  ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
}

static inline long long read_perf_ctr() {
  long long count;
  ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
  read(fd, &count, sizeof(long long));
  return count;
}
///////// end of perf ctr stuff


//////// start high-res timer stuff
using time_point = std::chrono::high_resolution_clock::time_point;

static time_point __start;

static void start_timer() {
  __start = std::chrono::high_resolution_clock::now();
}

static double read_timer() {
  time_point __end = std::chrono::high_resolution_clock::now();  
  std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(__end-__start);
  return time_span.count();
}
/////// end high-res timer stuff

#ifdef USE_PAPI
#include <papi.h>

class PerfData {
    float real_time, proc_time, mflops;
    long long flpins;
    uint64_t cycle_cnt;
    
 public:
    PerfData() {start();}
    void start() {
        if(PAPI_flops( &real_time, &proc_time, &flpins, &mflops)<PAPI_OK)
            throw "PAPI_flops failed";
        cycle_cnt = ::cycle_count();
    }
    void stop() {
        cycle_cnt = ::cycle_count() - cycle_cnt;
        if(PAPI_flops( &real_time, &proc_time, &flpins, &mflops)<PAPI_OK)
            throw "PAPI_flops failed";
    }
    double cpu_time() {return proc_time;}
    double cycle_count() {return double((cycle_cnt<<14)>>14);}
    double flop_count() {return double(flpins);}
    double gflops() {return 1e-3*mflops;}
};
#else
class PerfData {
    time_point __start;
    uint64_t __cycle_cnt;
    double proc_time;
    
 public:
    PerfData() {start();}
    void start() {
        __start = std::chrono::high_resolution_clock::now();            
        __cycle_cnt = ::cycle_count();
    }
    void stop() {
        __cycle_cnt = ::cycle_count() - __cycle_cnt;
        time_point __end = std::chrono::high_resolution_clock::now();  
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(__end-__start);
        proc_time = time_span.count();
    }
    double cpu_time() {return proc_time;}
    double cycle_count() {return double((__cycle_cnt<<14)>>14);}
    double flop_count() {return -1.0;}
    double gflops() {return -1.0;}
    double frequency() {
        time_point start = std::chrono::high_resolution_clock::now();            
        uint64_t cycle_cnt = ::cycle_count();
	int sum=0;
	for (int i=0; i<100000000; i++) sum+=i;
	cycle_cnt = ::cycle_count() - cycle_cnt;
        time_point end = std::chrono::high_resolution_clock::now();  
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
	if (sum==99) std::cout << sum << std::endl;
        return 1e-9*cycle_cnt/time_span.count();
    }	
};
#endif



#endif
