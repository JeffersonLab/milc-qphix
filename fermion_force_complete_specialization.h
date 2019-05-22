#ifndef _FERMION_FORCE_COMPLETE_H_
#define _FERMION_FORCE_COMPLETE_H_

#include "qcd_data_types.h"
#include "qphix_int.h"
#include "qphix.h"
#include "su3.h"

#if 1
#undef _mm_vector
#if QPHIX_PrecisionInt == 1
#if VECLEN == 16
#undef ARCH
#ifdef AVX512
#define ARCH avx512
#define _mm_vector __m512
#pragma message "Using ARCH avx512"
#else
#define ARCH mic
#define _mm_vector __m512
#endif
#elif VECLEN == 8
#undef ARCH
#ifdef AVX2
#define ARCH avx2
#define _mm_vector __m256
#warning "using avx2"
#else
#define ARCH avx
#define _mm_vector __m256
#endif
#elif VECLEN == 4
#undef ARCH
#define ARCH sse
#define _mm_vector __m128
#elif VECLEN == 1
#undef ARCH
#define ARCH scalar
#define _mm_vector float
#endif
#elif QPHIX_PrecisionInt == 2
#if VECLEN == 8
#undef ARCH
#ifdef AVX512
#define ARCH avx512
#define _mm_vector __m512d
#pragma message "Using ARCH avx512"
#else
#define ARCH mic
#define _mm_vector __m512d
#endif
#elif VECLEN == 4
#undef ARCH
#ifdef AVX2
#define ARCH avx2
#define _mm_vector __m256d
#else
#define ARCH avx
#define _mm_vector __m256d
#endif
#elif VECLEN == 2
#undef ARCH
#define ARCH sse
#define _mm_vector __m128d
#elif VECLEN == 1
#undef ARCH
#define ARCH scalar
#define _mm_vector double
#endif
#endif //PRECISION

#include "utils.h"

#undef QUOTEME
#define QUOTEME(M)       #M
#undef INCLUDE_FILE1
#define INCLUDE_FILE1(PRE,FPTYPE,VLEN,SLEN,POST)  QUOTEME(PRE ## FPTYPE ## _v ## VLEN ## _s ## SLEN ## POST) 
#undef INCLUDE_FILE2
#define INCLUDE_FILE2(PRE,FPTYPE,VLEN,SLEN,POST)  QUOTEME(PRE ## FPTYPE ## _ ## FPTYPE ## _v ## VLEN ## _s ## SLEN ## POST) 
#undef INCLUDE_FILE_VAR1
#define INCLUDE_FILE_VAR1(PRE,FPTYPE,VLEN,SLEN,POST) INCLUDE_FILE1(PRE,FPTYPE,VLEN,SLEN,POST)
#undef INCLUDE_FILE_VAR2
#define INCLUDE_FILE_VAR2(PRE,FPTYPE,VLEN,SLEN,POST) INCLUDE_FILE2(PRE,FPTYPE,VLEN,SLEN,POST)

	void ff_pack_face_dir_vec_noinline(
                QSU3M *giBase,
                fptype *rBuf,
                int dir
        );

	void ff_unpack_face_dir_vec_noinline(
                QSU3M *goBase,
                fptype *rBuf,
                int dir
        );

	void ff_pack_face_v_dir_vec_noinline(
                QSU3V *giBase,
                fptype *rBuf,
                int dir
        );

	void ff_unpack_face_v_dir_vec_noinline(
                QSU3V *goBase,
                fptype *rBuf,
                int dir
        );
	void ff_pack_face_v3_dir_vec_noinline(
                QSU3V *giBase,
                fptype *rBuf,
                int dir
        );

	void ff_unpack_face_v3_dir_vec_noinline(
                QSU3V *goBase,
                fptype *rBuf,
                int dir
        );

	void ff_pack_face_v3_dir_vec_noinline_special(
                QSU3V *giBase,
                fptype *rBuf,
                int dir
        );

	void ff_unpack_face_v3_dir_vec_noinline_special(
                QSU3V *goBase,
                fptype *rBuf,
                int dir
        );
#endif
#endif
