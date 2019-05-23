#ifndef _FF_BOUNDARY_H
#define _FF_BOUNDARY_H

#include "qphix.h"
#include "qphix_internal.h"

#define MYASSERT(cond) if(!(cond)) { printf("Rank %d: %s:%d Assertion failed\n\t%s evaluates false\n", myRank, __FILE__, __LINE__, #cond); exit(1);}

void ff_pack_face(int tid, int num_cores, int threads_per_core, QSU3M *gi, fptype *res, int cb, int dir, int fb);
void ff_unpack_face( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3M *res, int cb, int dir, int fb );
void ff_pack_and_send_boundaries(int tid, QSU3M *gi, int cb, int fb, int d);
void ff_recv_and_unpack_boundaries(int tid, QSU3M *gio, int cb, int fb, int d);
void ff_pack_face_v(int tid, int num_cores, int threads_per_core, QSU3V *gi, fptype *res, int cb, int dir, int fb);
void ff_unpack_face_v( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3V *res, int cb, int dir, int fb );
void ff_pack_and_send_boundaries_v(int tid, QSU3V *gi, int cb, int fb, int d);
void ff_recv_and_unpack_boundaries_v(int tid, QSU3V *gio, int cb, int fb, int d);
void ff_pack_face_v3(int tid, int num_cores, int threads_per_core, QSU3V *gi, fptype *res, int cb, int dir, int fb);
void ff_unpack_face_v3( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3V *res, int cb, int dir, int fb );
void ff_pack_and_send_boundaries_v3(int tid, QSU3V *gi, int cb, int fb, int d);
void ff_recv_and_unpack_boundaries_v3(int tid, QSU3V *gio, int cb, int fb, int d);
void ff_setup_comms_F();
void ff_setup_comms_D();
void ff_destroy_comms_F();
void ff_destroy_comms_D();
#endif
