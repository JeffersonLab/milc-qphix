#ifndef _GAUGE_FORCE_IMP_COMPLETE_H_
#define _GAUGE_FORCE_IMP_COMPLETE_H_

#include "qcd_data_types.h"
#include "qphix_int.h"

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


	/* Improved gauge force calculation body part */
	__attribute__((noinline))
	void gforce_imp_body_vec_noinline
	(
		double epsilonH,
		Hermit *outHBase[8], /* Output momentum addresses */ 
                Gauge *inGBase[8],   /* Input gauge addresses */
		//const char accumulate,
		//int gf_mask,
		const char isBoundary,
		/* Input momentum prefetches */
		/*
		Hermit *hxyBase, 
		Hermit *hpfBase2,
		Hermit *hpfBase3,
		Hermit *hpfBase4,
		const long hprefdist1,
		const long hprefdist2,
		const long hprefdist3,
		const long hprefdist4,
		*/
		Gauge *gxyBase,   /* Input gauge prefetches */
                Gauge *gpfBase2,
                Gauge *gpfBase3,
                Gauge *gpfBase4,
		double kappaS,
		double kappaR,
		double kappaB,
                const long gprefdist1,
                const long gprefdist2,
                const long gprefdist3,
                const long gprefdist4
	);

        /* Update momentum and gauge field body part */
        __attribute__((noinline))
        void update_mg_body_vec_noinline(
                double epsilonU,
                Gauge *gBase,
                Hermit *hBase,
		HermitHelperYZT *hYZTBase,
                const char accumulate,
                const long gprefdist1,
                const long gprefdist2,
                const long gprefdist3,
                const long gprefdist4
        );

	/* Send gauge fields (7-way) face */
        __attribute__((noinline))
        void gauge_pack_face_dir_vec_noinline(
                Gauge *giBase,
                fptype *lBuf,
                fptype *rBuf,
                int dir
        );

	/* Send gauge fields (1-way) face */
	__attribute__((noinline))
	void gauge_pack_face_dir_vec_noinline(
		Gauge *giBase,
		fptype *rBuf,
		int dir
	);

	/* unpack gauge fields (1-way) face */
        __attribute__((noinline))
        void gauge_unpack_face_dir_vec_noinline(
                Gauge *giBase,
                fptype *rBuf,
                int dir
	);

	/* Send momentum (1-way) face */
        __attribute__((noinline))
        void momentum_pack_face_dir_vec_noinline(
                Hermit *hiBase,
                fptype *rBuf,
                int dir
	);

	/* unpack momentum (1-way) face */
        __attribute__((noinline))
        void momentum_unpack_face_dir_vec_noinline(
                Hermit *hiBase,
                fptype *rBuf,
                int dir
        );

        /* unpack gauge force and update momentum */
        __attribute__((noinline))
        void gforce_unpack_face_dir_vec_noinline(
                fptype *hBuf,
                Hermit *outHBase,
                Hermit *hiBase[2],
		HermitHelperYZT *htmp,
                int dir,
		int hsize,
		const bool execm,
		const bool execp
        );

        /* unpack gauge field and pack momentum in face */
        __attribute__((noinline))
        void gauge_unpack_face_pack_gforce_dir_vec_noinline(
                fptype *hBuf[8],
		//fptype *hlBuf[6],
                fptype *rBuf[8],
                fptype *lBuf[8],
                Hermit *outHBase[8],
                Gauge *inGBase[7],
                int dir,
		//int gsize,
		int hsize[4],
		fptype kappaS, fptype kappaR, fptype kappaB,
		fptype epsilonH,
		const char isBoundary
        );

#endif

#endif
