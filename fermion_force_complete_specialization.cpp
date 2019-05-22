#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <immintrin.h>
#include <iostream>

#include "fermion_force_complete_specialization.h"

#define RE 0
#define IM 1


        __attribute__((noinline))
	void ff_pack_face_dir_vec_noinline(
                QSU3M *giBase,
                fptype *rBuf,
                int dir
        )
{
#include INCLUDE_FILE_VAR1(./ARCH/ff_pack_face_,fptype,VECLEN,,_18)
}

        __attribute__((noinline))
	void ff_unpack_face_dir_vec_noinline(
                QSU3M *goBase,
                fptype *rBuf,
                int dir
        )
{
#include INCLUDE_FILE_VAR1(./ARCH/ff_unpack_face_,fptype,VECLEN,,_18)
}

        __attribute__((noinline))
	void ff_pack_face_v_dir_vec_noinline(
                QSU3V *giBase,
                fptype *rBuf,
                int dir
        )
{
#include INCLUDE_FILE_VAR1(./ARCH/ff_pack_face_v_,fptype,VECLEN,,_18)
}

        __attribute__((noinline))
	void ff_unpack_face_v_dir_vec_noinline(
                QSU3V *goBase,
                fptype *rBuf,
                int dir
        )
{
#include INCLUDE_FILE_VAR1(./ARCH/ff_unpack_face_v_,fptype,VECLEN,,_18)
}

        __attribute__((noinline))
	void ff_pack_face_v3_dir_vec_noinline(
                QSU3V *giBase,
                fptype *rBuf,
                int dir
        )
{
#include INCLUDE_FILE_VAR1(./ARCH/ff_pack_face_v3_,fptype,VECLEN,,_18)
}

        __attribute__((noinline))
	void ff_unpack_face_v3_dir_vec_noinline(
                QSU3V *goBase,
                fptype *rBuf,
                int dir
        )
{
#include INCLUDE_FILE_VAR1(./ARCH/ff_unpack_face_v3_,fptype,VECLEN,,_18)
}

        __attribute__((noinline))
	void ff_pack_face_v3_dir_vec_noinline_special(
                QSU3V *giBase,
                fptype *rBuf,
                int dir
        )
{
#include INCLUDE_FILE_VAR1(./ARCH/ff_pack_face_v3_special_,fptype,VECLEN,,_18)
}

        __attribute__((noinline))
	void ff_unpack_face_v3_dir_vec_noinline_special(
                QSU3V *goBase,
                fptype *rBuf,
                int dir
        )
{
#include INCLUDE_FILE_VAR1(./ARCH/ff_unpack_face_v3_special_,fptype,VECLEN,,_18)
}

