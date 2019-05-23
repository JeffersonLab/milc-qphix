#ifndef _MISC_H_
#define _MISC_H_

#ifdef ENABLE_SEP_SAMPLING
#include "sampling_MIC.h"
#endif

#ifdef ENABLE_STREAMING_STORES
#pragma message "Using STREAM FAC 1"
#define STREAM_FACTOR  1.0
#else
#pragma message "Using STREAM FAC 2"
#define STREAM_FACTOR 2.0
#endif


#define SSC_MARK(mark_value)                  \
		{__asm  mov ebx, mark_value           \
		__asm  _emit 0x64                     \
		__asm  _emit 0x67                     \
		__asm  _emit 0x90}             


#ifdef ENABLE_SEP_SAMPLING
#define START_MARKER  VTResumeSampling()
#define STOP_MARKER   VTPauseSampling()
#elif defined(MCRT)
#define START_MARKER  sim_startPerf()
#define STOP_MARKER   sim_stopPerf()
#elif defined(SDE)
#define START_MARKER  SSC_MARK(0x111)
#define STOP_MARKER   SSC_MARK(0x222)
#else
#define START_MARKER
#define STOP_MARKER
#endif


#ifdef _OPENMP
#ifndef MCRT
#include <omp.h>
#else
// If in simulator define these as externs
extern "C" int omp_get_max_threads();
extern "C" int omp_get_num_threads();
extern "C" int omp_get_thread_num();
// Simulator doesn't support this call but anyway we want to use maximum threads 
// as defined by simulator config so we ignore this call
#define omp_set_num_threads(x)
#endif
#else // _OPENMP
#define omp_get_max_threads() (1)
#define omp_get_num_threads() (1)
#define omp_get_thread_num() (0)
#define omp_set_num_threads(x)
#endif // _OPENMP

#ifdef MCRT
#define _mm_malloc mcrtAlignedMalloc
#define _mm_free mcrtAlignedFree
#endif // MCRT

// This next part is to get large memory pages.
#ifdef USE_LARGE_PAGES

// Allocate and deallocate macros
#include <fcntl.h>
#include <sys/mman.h>
#ifdef MALLOC
#undef MALLOC
#endif
#define MALLOC(size,alignment) \
	mmap(NULL,size,PROT_READ|PROT_WRITE, MAP_ANONYMOUS | MAP_SHARED | MAP_HUGETLB | MAP_POPULATE, -1, 0)
#define FREE(addr) munmap(addr,0)
#else
// Otherwise use regular _mm_malloc
#ifndef MALLOC
#define MALLOC(size,alignment)  _mm_malloc(size,alignment);
#endif
#define FREE(addr) _mm_free(addr)
#endif

#endif // _MISC_H_
