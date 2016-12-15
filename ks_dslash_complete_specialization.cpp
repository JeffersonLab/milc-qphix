#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <immintrin.h>
#include <iostream>

//#include "qcd_data_types.h"
#include "utils.h"
#include "ks_dslash_complete_specialization.h"
#include "utils.h"

	__attribute__((noinline))
		void 
		ks_long_dslash_vec_noinline(
		int nthreads,
		KS *oBase,
		const Gauge18 *gBase,
		const Gauge *gllBase,
		KS *neighs[8],
		KS *neighs3[8],
		char accumulate,//const unsigned int accumulate[8],
		char accumulate3,//const unsigned int accumulate3[8],
		char isBoundary,//unsigned int isBoundary[8],
		char isBoundary3,//unsigned int isBoundary3[8],
		/*
		const unsigned int accumulate[8],
		const unsigned int accumulate3[8],
		unsigned int isBoundary[8],
		unsigned int isBoundary3[8],
		*/
		const int ind_next,
		const int x,
		const int y,
		const int t)
	{
/*
#if VECLEN==8
#warning "VECLEN=8"
#elif VECLEN==16
#warning "VECLEN=16"
#endif

	    printf("ks_long_dslash_vec_noinline: VECLEN = %d\n", VECLEN);
		printf("fptype = %d\n", sizeof(fptype)/sizeof(float));
#if QPHIX_PrecisionInt==1
		printf("QPHIX_PrecisionInt==1\n");
#else
		printf("QPHIX_PrecisionInt==2\n");
#endif
*/
#include INCLUDE_FILE_VAR2(./ARCH/ks_long_dslash_body_,fptype,VECLEN,,COMPRESSED_SUFFIX)
	}

	__attribute__((noinline))
	void ks_face_pack_dir(
		const KS *siBase,
		//const int si_prefdist,
		fptype *lBuf,
		fptype *rBuf,
		//const int hsprefdist,
		int dir)
	{
		if(dir == 0) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_dslash_face_pack_to_back_X_,fptype,VECLEN,,)
		}
		else if(dir == 1) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_dslash_face_pack_to_forw_X_,fptype,VECLEN,,)
		}
		else if(dir == 2) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_dslash_face_pack_to_back_Y_,fptype,VECLEN,,)
		}
		else if(dir == 3) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_dslash_face_pack_to_forw_Y_,fptype,VECLEN,,)
		}
		else if(dir == 4) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_dslash_face_pack_to_back_Z_,fptype,VECLEN,,)
		}
		else if(dir == 5) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_dslash_face_pack_to_forw_Z_,fptype,VECLEN,,)
		}
		else if(dir == 6) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_dslash_face_pack_to_back_T_,fptype,VECLEN,,)
		}
		else if(dir == 7) {
#include INCLUDE_FILE_VAR1(./ARCH/ks_dslash_face_pack_to_forw_T_,fptype,VECLEN,,)
		}
		else {
			printf("Invalid dir for pack boundary\n");
			exit(1);
		}
	}

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
		int dir) 
	{
		if(dir == 0) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_back_X_,fptype,VECLEN,,_18)
		}
		else if(dir == 1) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_forw_X_,fptype,VECLEN,,_18)
		}
		else if(dir == 2) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_back_Y_,fptype,VECLEN,,_18)
		}
		else if(dir == 3) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_forw_Y_,fptype,VECLEN,,_18)
		}
		else if(dir == 4) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_back_Z_,fptype,VECLEN,,_18)
		}
		else if(dir == 5) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_forw_Z_,fptype,VECLEN,,_18)
		}
		else if(dir == 6) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_back_T_,fptype,VECLEN,,_18)
		}
		else if(dir == 7) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_forw_T_,fptype,VECLEN,,_18)
		}
		else {
			printf("Invalid dir for unpack boundary\n");
			exit(1);
		}
	}

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
		int isBoundary3[]*/) 
	{
		if(dir == 0) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_back_X_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 1) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_forw_X_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 2) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_back_Y_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 3) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_forw_Y_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 4) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_back_Z_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 5) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_forw_Z_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 6) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_back_T_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 7) {
#include INCLUDE_FILE_VAR2(./ARCH/ks_dslash_face_unpack_from_forw_T_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else {
			printf("Invalid dir for unpack boundary\n");
			exit(1);
		}
	}

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

