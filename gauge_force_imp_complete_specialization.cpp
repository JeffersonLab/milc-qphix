#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <immintrin.h>
#include <iostream>

#include "gauge_force_imp_complete_specialization.h"

#define RE 0
#define IM 1

static
int index(int dir, int dir2)
{
  int dt = (dir<dir2?dir:dir2);
  if(dt==0) return (dir+dir2 - 1);
  else return (dir+dir2);
}

        /* Improved gauge force calculation body part */
	/* Update momentum body part */
        __attribute__((noinline))
        void gforce_imp_body_vec_noinline
        (
		double epsilonH,
                /* Output momentum address */
		Hermit *outHBase[8],
		/* Input gauge addresses */
		Gauge *inGBase[8],
/*
		Gauge *inGxback,
                Gauge *inGxforw,
                Gauge *inGyback,
                Gauge *inGyforw,
                Gauge *inGzback,
                Gauge *inGzforw,
                Gauge *inGtback,
                Gauge *inGtforw,
*/
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
        )
{
#include INCLUDE_FILE_VAR2(./ARCH/gf_imp_body_,fptype,VECLEN,,COMPRESSED_SUFFIX)
}

	/* Update gauge field body part */
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
	)
{
#include INCLUDE_FILE_VAR2(./ARCH/update_mg_body_,fptype,VECLEN,,COMPRESSED_SUFFIX)
}

	/* Send gauge fields face */
	__attribute__((noinline))
	void gauge_pack_face_dir_vec_noinline(
		Gauge *giBase,
		fptype *lBuf,
		fptype *rBuf,
		int dir
	)
{
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_face_,fptype,VECLEN,,_12)
/*
  if(dir==0) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_to_back_X_,fptype,VECLEN,,_12)
  }
  else if(dir==1) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_to_forw_X_,fptype,VECLEN,,_12)
  }
  else if(dir==2) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_to_back_Y_,fptype,VECLEN,,_12)
  }
  else if(dir==3) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_to_forw_Y_,fptype,VECLEN,,_12)
  }
  else if(dir==4) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_to_back_Z_,fptype,VECLEN,,_12)
  }
  else if(dir==5) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_to_forw_Z_,fptype,VECLEN,,_12)
  }
  else if(dir==6) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_to_back_T_,fptype,VECLEN,,_12)
  }
  else if(dir==7) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_7w_pack_to_forw_T_,fptype,VECLEN,,_12)
  }
*/
}

        __attribute__((noinline))
        void gauge_pack_face_dir_vec_noinline(
                Gauge *giBase,
                fptype *rBuf,
                int dir
        )
{
  if(dir==0) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_pack_to_back_X_,fptype,VECLEN,,_12)
  }
  else if(dir==1) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_pack_to_forw_X_,fptype,VECLEN,,_12)
  }
  else if(dir==2) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_pack_to_back_Y_,fptype,VECLEN,,_12)
  }
  else if(dir==3) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_pack_to_forw_Y_,fptype,VECLEN,,_12)
  }
  else if(dir==4) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_pack_to_back_Z_,fptype,VECLEN,,_12)
  }
  else if(dir==5) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_pack_to_forw_Z_,fptype,VECLEN,,_12)
  }
  else if(dir==6) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_pack_to_back_T_,fptype,VECLEN,,_12)
  }
  else if(dir==7) {
#include INCLUDE_FILE_VAR1(./ARCH/gauge_pack_to_forw_T_,fptype,VECLEN,,_12)
  }
}

        __attribute__((noinline))
        void gauge_unpack_face_dir_vec_noinline(
                Gauge *giBase,
                fptype *rBuf,
                int dir
        )
{

}

	__attribute__((noinline))
	void momentum_pack_face_dir_vec_noinline(
		Hermit *hiBase,
		fptype *rBuf,
		int dir
	)
{
#include INCLUDE_FILE_VAR1(./ARCH/momentum_pack_face_,fptype,VECLEN,,)
}

	__attribute__((noinline))
	void momentum_unpack_face_dir_vec_noinline(
		Hermit *hoBase,
                fptype *rBuf,
                int dir
	)
{
#include INCLUDE_FILE_VAR1(./ARCH/momentum_unpack_face_,fptype,VECLEN,,)
}

	/* unpack gauge force and update momentum */
	__attribute__((noinline))
	void gforce_unpack_face_dir_vec_noinline(
		fptype *hBuf,
		Hermit *hBase,
		Hermit *hiBase[2],
		HermitHelperYZT *hYZTBase,
		int dir,
		int hsize,
		const bool execm,
		const bool execp
	)
{
#include INCLUDE_FILE_VAR1(./ARCH/gforce_unpack_face_,fptype,VECLEN,,)
}

	/* unpack gauge field and pack gauge force in face */
	__attribute__((noinline))
	void gauge_unpack_face_pack_gforce_dir_vec_noinline(
		fptype *hBuf[8],
		fptype *rBuf[8],
		fptype *lBuf[8],
		Hermit *outHBase[8],
		Gauge *inGBase[7],
		int dir,
		int hsize[4], 
		/*int gsize,
		int hsize,*/
		fptype kappaS, fptype kappaR, fptype kappaB,
		fptype epsilonH,
		const char isBoundary
	)
{
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_,fptype,VECLEN,,COMPRESSED_SUFFIX)
#if 0
  if(dir==0) {
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_to_back_X_,fptype,VECLEN,,COMPRESSED_SUFFIX)
  }
  else if(dir==1) {
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_to_forw_X_,fptype,VECLEN,,COMPRESSED_SUFFIX)
  }
  else if(dir==2) {
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_to_back_Y_,fptype,VECLEN,,COMPRESSED_SUFFIX)
  }
  else if(dir==3) {
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_to_forw_Y_,fptype,VECLEN,,COMPRESSED_SUFFIX)
  }
  else if(dir==4) {
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_to_back_Z_,fptype,VECLEN,,COMPRESSED_SUFFIX)
  }
  else if(dir==5) {
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_to_forw_Z_,fptype,VECLEN,,COMPRESSED_SUFFIX)
  }
  else if(dir==6) {
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_to_back_T_,fptype,VECLEN,,COMPRESSED_SUFFIX)
  }
  else if(dir==7) {
#include INCLUDE_FILE_VAR2(./ARCH/gauge_unpack_face_pack_gforce_to_forw_T_,fptype,VECLEN,,COMPRESSED_SUFFIX)
  }
#endif
}

