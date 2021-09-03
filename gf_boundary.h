#ifndef _GF_BOUNDARY_H
#define _GF_BOUNDARY_H

#include "ks_boundary.h"

void pack_gauge_face_dir(int tid, Gauge *gi, int dir, int cb);
void pack_gauge_face_dir(int tid, fptype *gi, int dir, int cb);
void pack_momentum_face_dir(int tid, Hermit *hi, int dir, int cb);
void pack_momentum_face_dir(int tid, fptype *hi, int dir, int cb);
void unpack_gauge_face_dir(int tid, Gauge *go, int dir, int cb);
void unpack_gauge_face_dir(int tid, fptype *go, int dir, int cb);
void unpack_momentum_face_dir(int tid, Hermit *ho, int dir, int cb);
void unpack_momentum_face_dir(int tid, fptype *ho, int dir, int cb);
void gf_pack_and_send_boundaries_u(int tid, Gauge *gi, int cb);
void gf_recv_and_unpack_boundaries_h(int tid, Hermit *hio, HermitHelperYZT *htmp, Gauge *gio, int cb);
void gf_recv_and_unpack_u_and_send_boundaries_h(int tid, Hermit *hio, HermitHelperYZT *htmp, Gauge *gio, fptype kappaS, fptype kappaR, fptype kappaB, fptype epsilonH, int cb);
void gf_print_boundary_timings(unsigned long long t_gf);
void gf_setup_comms_F();
void gf_setup_comms_D();
void gf_destroy_comms_F();
void gf_destroy_comms_D();
#endif
