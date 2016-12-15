#ifndef _KS_CONFIG_H_
#define _KS_CONFIG_H_
#if 0
/* Define the global vars and macros that we use in the micro-benchmark */
#define BY 2
#define BZ 2
#define PRECISION 1

#ifdef __AVX512F__
#define VECLEN 16
#define ARCH avx512
#elif defined(__AVX2__)
#define VECLEN 8
#define ARCH avx2
#elif defined(__AVX__)
#define VECLEN 8
#define ARCH avx
#elif defined(__MIC__)
#define VECLEN 16
#define ARCH mic
#else
#error "Either -mavx or -mmic must be used"
#endif

#define CODEGEN_PATH ./
#endif
//#define fptype float
#define XY_PAD  (1)
#define XYZ_PAD (0)
#define MALLOC _mm_malloc

#endif
