#ifndef _KS_DSLASH_COMPLETE_H_
#define _KS_DSLASH_COMPLETE_H_

#include "qcd_data_types.h"

#if 1
#if QPHIX_PrecisionInt == 1
#if VECLEN == 16
#undef ARCH
#ifdef AVX512
#define ARCH avx512
#pragma message "Using ARCH avx512"
#else
#define ARCH mic
#endif
#elif VECLEN == 8
#undef ARCH
#ifdef AVX2
#define ARCH avx2
#warning "using avx2"
#else
#define ARCH avx
#endif
#elif VECLEN == 4
#undef ARCH
#define ARCH sse
#elif VECLEN == 1
#undef ARCH
#define ARCH scalar
#endif
#elif QPHIX_PrecisionInt == 2
#if VECLEN == 8
#undef ARCH
#ifdef AVX512
#define ARCH avx512
#pragma message "Using ARCH avx512"
#else
#define ARCH mic
#endif
#elif VECLEN == 4
#undef ARCH
#ifdef AVX2
#define ARCH avx2
#else
#define ARCH avx
#endif
#elif VECLEN == 2
#undef ARCH
#define ARCH sse
#elif VECLEN == 1
#undef ARCH
#define ARCH scalar
#endif
#endif //PRECISION

#include "utils.h"

//./mic/dslash_body_,float_float,_v16_s4_12
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

// THESE ARE HERE SO WE CAN LOOK AT THE ASSEMBLER OUTPUT
//
//extern "C" {
#if 0
	__attribute__((noinline))
		void 
		ks_dslash_vec_noinline(
		const KS *xyBase,
		const KS *zbBase,
		const KS *zfBase,
		const KS *tbBase,
		const KS *tfBase,
		KS *oBase,
		const Gauge *gBase,
		const int xbOffs[VECLEN],
		const int xfOffs[VECLEN],
		const int ybOffs[VECLEN],
		const int yfOffs[VECLEN],
		const int offs[VECLEN],
		const int gOffs[VECLEN],
		const int siprefdist1,
		const int siprefdist2,
		const int siprefdist3,
		const int siprefdist4,
		const int gprefdist,
		const int pfyOffs[VECLEN],
		const KS *pfBase2,
		const KS *pfBase3,
		const KS *pfBase4,
		const int x,
		const int y,
		const int t,
		const unsigned int accumulate[8])
	{
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_body_,fptype,VECLEN,,COMPRESSED_SUFFIX)
	}
#endif

	__attribute__((noinline))
		void 
		ks_long_dslash_vec_noinline(
		int nthreads,
		KS *oBase,
		const Gauge18 *gBase,
		const Gauge *gllBase,
		KS *neighs[8],
		KS *neighs3[8],
		/*
		const unsigned int accumulate[8],
		const unsigned int accumulate3[8],
		unsigned int isBoundary[8],
		unsigned int isBoundary3[8],
		*/
		char accumulate,
		char accumulate3,
		char isBoundary,
		char isBoundary3,
		const int ind_next,
		const int x,
		const int y,
		const int t);

	__attribute__((noinline))
	void ks_face_pack_dir(
		const KS *siBase,
		//const int si_prefdist,
		fptype *lBuf,
		fptype *rBuf,
		//const int hsprefdist,
		int dir);

	__attribute__((noinline))
	void ks_dslash_face_unpack18_dir(
		const fptype *lBuf,
		const fptype *rBuf,
		const Gauge18 *gBase,
		KS *oBase,
		//const int hsprefdist,
		//const int gprefdist,
		//const int soprefdist,
		const fptype beta,
		const int x,
		const int y,
		const int t,
		int dir); 

	__attribute__((noinline))
	void ks_dslash_face_unpack_dir(
		const fptype *lBuf,
		const fptype *rBuf,
		const Gauge *gBase,
		KS *oBase,
		//const int hsprefdist,
		//const int gprefdist,
		//const int soprefdist,
		const fptype beta,
		const int x,
		const int y,
		const int t,
		int dir/*,
		int isBoundary3[]*/);

#if 0
	__attribute__((noinline))
	void ks_long_dslash_face_unpack_dir(
		const fptype *inbuf,
		const fptype *inbuf_ll,
		const Gauge18 *gBase,
		const Gauge *gllBase,
		KS *oBase,
		const int gOffs[VECLEN],
		const int gllOffs[VECLEN],
		const int offs[VECLEN],
		const int hsprefdist,
		const int gprefdist,
		const int soprefdist,
		const int x,
		const int y,
		const int t,
		unsigned int mask,
		unsigned int mask_ll,
		int dir,
		int width) 
	{
		if(dir == 0) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_face_unpack_from_back_X_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 1) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_face_unpack_from_forw_X_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 2) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_face_unpack_from_back_Y_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 3) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_face_unpack_from_forw_Y_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 4) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_face_unpack_from_back_Z_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 5) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_face_unpack_from_forw_Z_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 6) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_face_unpack_from_back_T_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 7) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_face_unpack_from_forw_T_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else {
			printf("Invalid dir for unpack boundary\n");
			exit(1);
		}
	}

	__attribute__((noinline))
	void rephase_vec_noinline(
		//const bool flag,
		Gauge *gBase,
		const int gOffs[VECLEN],
		const int gprefdist,
		const int x,
		const int y,
		const int t)
	{
#include INCLUDE_FILE_VAR1(./ARCH/ks_rephase_body_,fptype,VECLEN,,COMPRESSED_SUFFIX)
	}

	// load long links (assuming gauge contains phase)
	__attribute((noinline))
	void load_longlinks_vec_noinline(
		const Gauge *gBase,
		const Gauge *gzfBase,
		const Gauge *gtfBase,
		Gauge *gllBase,
		Gauge *gllzfBase,
		Gauge *glltfBase,
		const int   gOffs[VECLEN],
		const int   gxfOffs[VECLEN],
		const int   gyfOffs[VECLEN],
		const int   giprefdist1,
		const int   giprefdist2,
		const int   gprefdist,
		const int   gpfyOffs[VECLEN],
		const Gauge *gpfBase1,
		const Gauge *gpfBase2,
		const int x,
		const int y,
		const int t,
		int dirmask[4])
	{
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_body_,fptype,VECLEN,,COMPRESSED_SUFFIX)
	}

	__attribute__((noinline))
	void ks_gauge_pack_dim(
		const Gauge *gBase,
		const int   gOffs[VECLEN],
		fptype *outbuf,
		const int gprefdist,
		int dim)
	{
		if(dim == 0) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_face_pack_to_forw_X_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dim == 1) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_face_pack_to_forw_Y_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dim == 2) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_face_pack_to_forw_Z_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dim == 3) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_face_pack_to_forw_T_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else {
			printf("Invalid dir for pack boundary\n");
			exit(1);
		}
	}

	__attribute__((noinline))
	void ks_gauge_unpack_dim(
		const Gauge *gBase,
		Gauge *gllBase,
		const int   gOffs[VECLEN],
		const int   gfOffs[VECLEN],
		const fptype *inbuf,
		fptype *outbuf,
		const int gprefdist,
		const int x,
		const int y,
		const int t,
		int dim)
	{
		if(dim == 0) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_face_unpack_from_forw_X_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dim == 1) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_face_unpack_from_forw_Y_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dim == 2) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_face_unpack_from_forw_Z_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dim == 3) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_load_longlinks_face_unpack_from_forw_T_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else {
			printf("Invalid dir for pack boundary\n");
			exit(1);
		}
	}
#endif
//}
#endif

#endif
