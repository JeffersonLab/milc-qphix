#include "qcd_data_types.h"

#if 1
#if PRECISION == 1
#if VECLEN == 16
#ifdef AVX512
#define ARCH avx512
#pragma message "Using ARCH avx512"
#else
#define ARCH mic
#endif
#elif VECLEN == 8
#ifdef AVX2
#define ARCH avx2
#else
#define ARCH avx
#endif
#elif VECLEN == 4
#define ARCH sse
#elif VECLEN == 1
#define ARCH scalar
#endif
#elif PRECISION == 2
#if VECLEN == 8
#ifdef AVX512
#define ARCH avx512
#pragma message "Using ARCH avx512"
#else
#define ARCH mic
#endif
#elif VECLEN == 4
#ifdef AVX2
#define ARCH avx2
#else
#define ARCH avx
#endif
#elif VECLEN == 2
#define ARCH sse
#elif VECLEN == 1
#define ARCH scalar
#endif
#endif //PRECISION

#include "utils.h"

//./mic/dslash_plus_body_,float_float,_v16_s4_12
#define QUOTEME(M)       #M
#define INCLUDE_FILE(PRE,FPTYPE,VLEN,SLEN,POST)  QUOTEME(PRE ## FPTYPE ## _ ## FPTYPE ## _v ## VLEN ## _s ## SLEN ## POST) 
#define INCLUDE_FILE_VAR(PRE,FPTYPE,VLEN,SLEN,POST) INCLUDE_FILE(PRE,FPTYPE,VLEN,SLEN,POST)

// THESE ARE HERE SO WE CAN LOOK AT THE ASSEMBLER OUTPUT
//
extern "C" {
	__attribute__((noinline))
		void 
		dslash_plus_vec_noinline(
		Spinor *oBase,
		const Gauge *gBase,
		Spinor* neighs[8],
		unsigned int accumulate[8],
		unsigned int isBoundary[8],
		const int ind_next,
		const fptype coeff_s,
		const fptype coeff_t_f,
		const fptype coeff_t_b)
	{
#include INCLUDE_FILE_VAR(./ARCH/dslash_plus_body_,fptype,VECLEN,,COMPRESSED_SUFFIX)
	}

#if 1
	__attribute__((noinline))
	void face_proj_dir_plus(
		const Spinor *siBase,
		//const int si_prefdist,
		fptype *lBuf,
		fptype *rBuf,
		//const int hsprefdist,
		int dir)
	{
		if(dir == 0) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_pack_to_back_X_plus_,fptype,VECLEN,,)
		}
		else if(dir == 1) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_pack_to_forw_X_plus_,fptype,VECLEN,,)
		}
		else if(dir == 2) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_pack_to_back_Y_plus_,fptype,VECLEN,,)
		}
		else if(dir == 3) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_pack_to_forw_Y_plus_,fptype,VECLEN,,)
		}
		else if(dir == 4) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_pack_to_back_Z_plus_,fptype,VECLEN,,)
		}
		else if(dir == 5) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_pack_to_forw_Z_plus_,fptype,VECLEN,,)
		}
		else if(dir == 6) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_pack_to_back_T_plus_,fptype,VECLEN,,)
		}
		else if(dir == 7) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_pack_to_forw_T_plus_,fptype,VECLEN,,)
		}
		else {
			printf("Invalid dir for pack boundary\n");
			exit(1);
		}
	}

	__attribute__((noinline))
	void face_finish_dir_plus(
		const fptype *lBuf,
		const fptype *rBuf,
		const Gauge *gBase,
		Spinor *oBase,
		//const int hsprefdist,
		//const int gprefdist,
		//const int soprefdist,
		const fptype beta,
		int dir) 
	{
		if(dir == 0) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_unpack_from_back_X_plus_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 1) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_unpack_from_forw_X_plus_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 2) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_unpack_from_back_Y_plus_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 3) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_unpack_from_forw_Y_plus_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 4) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_unpack_from_back_Z_plus_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 5) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_unpack_from_forw_Z_plus_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 6) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_unpack_from_back_T_plus_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else if(dir == 7) {
#include INCLUDE_FILE_VAR(./ARCH/dslash_face_unpack_from_forw_T_plus_,fptype,VECLEN,,COMPRESSED_SUFFIX)
		}
		else {
			printf("Invalid dir for unpack boundary\n");
			exit(1);
		}
	}
#endif
}
#endif
