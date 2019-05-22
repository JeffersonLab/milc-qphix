#include <unistd.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

//#define NO_COMPUTE 1

#include "ks_globals.h"        /* All global values and macros */
#include "ff_boundary.h"
#include "qphix_internal.h"
#include "fermion_force_complete_specialization.h"
#include "misc.h"

#if QPHIX_PrecisionInt == 1
#error "QPHIX_PrecisionInt SP not defined/supported!"
#undef ff_setup_comms
#define ff_setup_comms ff_setup_comms_F
#undef ff_destroy_comms
#define ff_destroy_comms ff_destroy_comms_F
#elif QPHIX_PrecisionInt == 2
#undef ff_setup_comms
#define ff_setup_comms ff_setup_comms_D
#undef ff_destroy_comms
#define ff_destroy_comms ff_destroy_comms_D
#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif
#define MYASSERT(cond) if(!(cond)) { printf("Rank %d: %s:%d Assertion failed\n\t%s evaluates false\n", myRank, __FILE__, __LINE__, #cond); exit(1);}

extern int PadBound;
extern int PadNeigh;
extern char * BoundTable;
extern unsigned int * NeighTable;

#ifndef ENABLE_MPI
/* link & force remappings */
void ff_pack_face(int tid, int num_cores, int threads_per_core, QSU3M *gi, fptype *res, int cb, int dir, int fb){}
void ff_unpack_face( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3M *res, int cb, int dir, int fb ){}
void ff_pack_face_v(int tid, int num_cores, int threads_per_core, QSU3V *gi, fptype *res, int cb, int dir, int fb){}
void ff_unpack_face_v( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3V *res, int cb, int dir, int fb ){}
void ff_pack_face_v3(int tid, int num_cores, int threads_per_core, QSU3V *gi, fptype *res, int cb, int dir, int fb){}
void ff_unpack_face_v3( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3V *res, int cb, int dir, int fb ){}

void ff_pack_and_send_boundaries(int tid, QSU3M *gi, int cb, int fb, int d){}
void ff_recv_and_unpack_boundaries(int tid, QSU3M *gio, int cb, int fb, int d){}
void ff_pack_and_send_boundaries_v(int tid, QSU3V *gi, int cb, int fb, int d){}
void ff_recv_and_unpack_boundaries_v(int tid, QSU3V *gio, int cb, int fb, int d){}
void ff_pack_and_send_boundaries_v3(int tid, QSU3V *gi, int cb, int fb, int d){}
void ff_recv_and_unpack_boundaries_v3(int tid, QSU3V *gio, int cb, int fb, int d){}
void ff_setup_comms_D(){}
void ff_destroy_comms_D(){}
#else // ENABLE_MPI is defined
#ifndef ASSUME_MULTINODE
#include <mpi.h>
#endif

const int nGX = (VECLEN < 2 ? 1 : 2);
const int nGY = (VECLEN < 4 ? 1 : 2);
const int nGZ = (VECLEN < 8 ? 1 : 2);
const int nGT = (VECLEN < 16 ? 1 : 2);

extern int neigh_ranks[8];
extern int myCoord[4];
extern int Ls[4];
static int bufSizeL[4];
static int bufSizeV[4];
static int bufSizeV3[4];
void *fcommsBuf[16];
void *vcommsBuf[16];
void *v3commsBuf[16];

#ifndef ASSUME_MULTINODE
static MPI_Request reqSendsFF[8];
static MPI_Request reqRecvsFF[8];
static MPI_Request reqSendsV[8];
static MPI_Request reqRecvsV[8];
#endif

static int nCommsDirsFF = 0;

void ff_setup_comms_F()
{
	printf("SP FF not supported\n");
	exit(1);
}

void ff_setup_comms_D()
{
        int gDims[4] = {Gx, Gy, Gz, Gt};
        int lDims[4] = {Gx, Gy, Gz, Gt};
	for(int i = 0; i < 4; i++) {
	    lDims[i] = gDims[i] / geometry[i];
	}

        int numProcs = geometry[0] * geometry[1] * geometry[2] * geometry[3];
        MYASSERT(numProcs == nRanks);
        for(int i = 0; i < 4; i++) {
                bufSizeL[i] = 1;
                bufSizeV[i] = 1;
                bufSizeV3[i] = 1;
                for(int j = 0; j < 4; j++) if(j != i) {
		    bufSizeL[i] *= lDims[j]; 
		    bufSizeV[i] *= lDims[j]; 
		    bufSizeV3[i] *= lDims[j]; 
		}
		bufSizeL[i] *= 18/2 /* Links */; 
		bufSizeV[i] *= 6/2 /* Vectors */; 
		bufSizeV3[i] *= 3*6/2 /* 3rd neigh Vectors */; 
	}

        int alignedBufSizeL[4] = {0};
        int totalBufSizeL = 0;
        int alignedBufSizeV[4] = {0};
        int totalBufSizeV = 0;
        int alignedBufSizeV3[4] = {0};
        int totalBufSizeV3 = 0;
        int pgSize = 4096 / sizeof(fptype);
        for(int i = 0; i < 4; i++) {
                alignedBufSizeL[i] = (bufSizeL[i] + (pgSize-1)) & ~(pgSize-1);
                alignedBufSizeV[i] = (bufSizeV[i] + (pgSize-1)) & ~(pgSize-1);
                alignedBufSizeV3[i] = (bufSizeV3[i] + (pgSize-1)) & ~(pgSize-1);
                if(!local_dir[i]) {
                        totalBufSizeL += alignedBufSizeL[i];
                        totalBufSizeV += alignedBufSizeV[i];
                        totalBufSizeV3 += alignedBufSizeV3[i];
                }
        }
        fptype *ftmpBuf = NULL;
        if(totalBufSizeL > 0) {
		ftmpBuf = (fptype*)MALLOC(4 * totalBufSizeL * sizeof(fptype), 4096);
		MYASSERT(ftmpBuf != NULL);
	}
        fptype *vtmpBuf = NULL;
        if(totalBufSizeV > 0) {
		vtmpBuf = (fptype*)MALLOC(4 * totalBufSizeV * sizeof(fptype), 4096);
		MYASSERT(vtmpBuf != NULL);
	}
        fptype *v3tmpBuf = NULL;
        if(totalBufSizeV3 > 0) {
		v3tmpBuf = (fptype*)MALLOC(4 * totalBufSizeV3 * sizeof(fptype), 4096);
		MYASSERT(v3tmpBuf != NULL);
	}
        for(int i = 0; i < 4; i++) {
                if(!local_dir[i]) {
                        nCommsDirsFF += 2;
#ifndef ASSUME_MULTINODE
                        for(int j = 0; j < 4; j++) {
                                fcommsBuf[4*i+j] = (void*)ftmpBuf;
				ftmpBuf += alignedBufSizeL[i];
                                vcommsBuf[4*i+j] = (void*)vtmpBuf;
				vtmpBuf += alignedBufSizeV[i];
                                v3commsBuf[4*i+j] = (void*)v3tmpBuf;
				v3tmpBuf += alignedBufSizeV3[i];
                        }
#endif
                }
                else {
                        fcommsBuf[4*i] = 0x00;
                        fcommsBuf[4*i+1] = 0x00;
                        fcommsBuf[4*i+2] = 0x00;
                        fcommsBuf[4*i+3] = 0x00;
                        vcommsBuf[4*i] = 0x00;
                        vcommsBuf[4*i+1] = 0x00;
                        vcommsBuf[4*i+2] = 0x00;
                        vcommsBuf[4*i+3] = 0x00;
                        v3commsBuf[4*i] = 0x00;
                        v3commsBuf[4*i+1] = 0x00;
                        v3commsBuf[4*i+2] = 0x00;
                        v3commsBuf[4*i+3] = 0x00;
                }
        }
#ifndef ASSUME_MULTINODE
        for(int d = 0; d < 8; d++) {
                reqSendsFF[d] = reqRecvsFF[d] = MPI_REQUEST_NULL;
                reqSendsV[d] = reqRecvsV[d] = MPI_REQUEST_NULL;
        }
#endif
}

void ff_pack_face(int tid, int num_cores, int threads_per_core, QSU3M *gi, fptype *res, int cb, int dir, int fb)
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;

        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pktsize = VECLEN;
	if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	MYASSERT(18*pktsize*npkts == bufSizeL[dir]);

        const unsigned int *NTNow = NeighTable+cb*PadNeigh;

#ifdef MT_EXECUTION
// this code path is to taken if FF outer threading is implemented (similar to GF or CG).
// this isn't implemented or invoked because FF global threading isn't implemented.
	int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;
        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
#else
// prepare the send buffers in parallel using OpenMP parallel constructs
#pragma omp parallel for
        for(int pkt = 0 ; pkt < npkts; pkt+=1) {
#endif

                int coords[4];
                int tmp = pkt;
                for(int j = 0; j < 4; j++) {
                        int tmp1 = tmp / lens[j];
                        coords[j] = tmp - tmp1*lens[j];
                        tmp = tmp1;
                }

// fb = 0 means forward direction, = 1 means backward. Set coords for that direction accordingly.
                coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
                int x = coords[0];
                int yblock = coords[1];
                int z = coords[2];
                int t = coords[3];
                int y = yblock * nyg_pack;
// the xodd logic is used to determine if the first element in X direction is Even or Odd. Lattice at row 0 starts with Even.
// specifically for dir = 0, which is X, the row parameter is used as an offset to correctly index into the lattice.
                int xodd = (z + t + cb) & 1;
                if(dir == 0) {
                        int row = (fb == 1 ? 1-xodd : xodd); 
                        y += row;
                }
		QSU3M *giBase = &gi[t*Pxyz + z*Pxy + y*Vxh + x];
                fptype * __restrict rBuf = &res[18*pktsize*pkt];
		ff_pack_face_dir_vec_noinline(giBase, rBuf, dir*2+fb);
	}

}

void ff_unpack_face( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3M *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
	int pktsize = VECLEN;
        if((1 << dir) < VECLEN) pktsize = VECLEN/2;
// shift is used as an offset into the lattice index to correctly update the output lattice.
	int shifts[4] = {Vxh-1, Vxh*(Vy-1), Vxh*Vy*(Vz-1), Vxh*Vy*Vz*(Vt-1)}, shift;
// shift is +ve if forward direction data lookup and -ve if backward.
	shift = shifts[dir] * (fb == 0 ? 1 : -1);

#ifdef MT_EXECUTION
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
#else
#pragma omp parallel for
        for(int pkt = 0 ; pkt < npkts; pkt+=1) {
#endif
                int coords[4];
                int tmp = pkt;
                for(int j = 0; j < 4; j++) {
                        int tmp1 = tmp / lens[j];
                        coords[j] = tmp - tmp1*lens[j];
                        tmp = tmp1;
                }

                coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
                int x = coords[0];
                int yblock = coords[1];
                int z = coords[2];
                int t = coords[3];
                int y = yblock * nyg_pack;
                int xodd = (z + t + cb) & 1;

                if(dir == 0) {
                        int row = (fb == 1 ? 1-xodd : xodd);
                        y += row;
                }

		QSU3M *goBase = &res[t*Pxyz + z*Pxy + y*Vxh + x];
		fptype * __restrict rBuf = &gi[18*pktsize*pkt];

// note the use of shift into the lattice index when invoking the unpack routine.
		ff_unpack_face_dir_vec_noinline(goBase+shift, rBuf, dir*2+fb);
	}
}

/* cb==1 : EVEN */
// that is cb = 1 implies that we are operating on Even data and gathering (pack) is of Off data.
// the value of cb is supplied from the shift Matrix, Vector function is layout.cpp file.
// d is direction. 0 is X, 1 is Y, 2 is Z, 3 is T. Opposite direction is indicated by fb parameter.
// tid is always 0 because we are doing threading in the pack/ unpack routines and not global FF threading.
void ff_pack_and_send_boundaries(int tid, QSU3M *gi, int cb, int fb, int d)
{
// local_dir array is set or unset depending on whether the data is onsite or offsite. This is part of QPhiX Init code.
    if(!local_dir[d]) {
	if(tid == 0) {
		MYASSERT(MPI_Irecv((void *)fcommsBuf[4*d+fb], bufSizeL[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+(fb^1)], fb, MPI_COMM_THISJOB, &reqRecvsFF[2*d+fb]) == MPI_SUCCESS);
	}
	ff_pack_face(tid, NCores, n_threads_per_core, gi, (fptype *)fcommsBuf[2+4*d+fb], cb, d, fb);
	if(tid == 0) {
		MYASSERT(MPI_Isend((void *)fcommsBuf[2+4*d+fb], bufSizeL[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], fb, MPI_COMM_THISJOB, &reqSendsFF[2*d+fb]) == MPI_SUCCESS);
	}
    }
}

void ff_recv_and_unpack_boundaries(int tid, QSU3M *gio, int cb, int fb, int d)
{
   if(!local_dir[d]) {
	if(tid==0) {
		MYASSERT(MPI_Wait(&reqRecvsFF[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
	}
	ff_unpack_face(tid, NCores, n_threads_per_core, (fptype *)fcommsBuf[4*d+fb], gio, cb, d, fb);
   }
   if(!local_dir[d]) {
	if(tid == 0) {
		MYASSERT(MPI_Wait(&reqSendsFF[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
	}
   }
}

// destroy the comm buffers created.
void ff_destroy_comms_D()
{
    if(fcommsBuf[0]!=0x00)
      FREE(fcommsBuf[0]);
    for(int i=0; i<16; ++i)
    {
      fcommsBuf[i]=0x00;
    }
    if(vcommsBuf[0]!=0x00)
      FREE(vcommsBuf[0]);
    for(int i=0; i<16; ++i)
    {
      vcommsBuf[i]=0x00;
    }
    if(v3commsBuf[0]!=0x00)
      FREE(v3commsBuf[0]);
    for(int i=0; i<16; ++i)
    {
      v3commsBuf[i]=0x00;
    }
}

// this is code for vector pack. much similar to matrix pack above...
void ff_pack_face_v(int tid, int num_cores, int threads_per_core, QSU3V *gi, fptype *res, int cb, int dir, int fb)
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;

        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pktsize = VECLEN;
	if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	MYASSERT(6*pktsize*npkts == bufSizeV[dir]);

        const unsigned int *NTNow = NeighTable+cb*PadNeigh;

#ifdef MT_EXECUTION
	int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;
        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
#else
#pragma omp parallel for
        for(int pkt = 0 ; pkt < npkts; pkt+=1) {
#endif

                int coords[4];
                int tmp = pkt;
                for(int j = 0; j < 4; j++) {
                        int tmp1 = tmp / lens[j];
                        coords[j] = tmp - tmp1*lens[j];
                        tmp = tmp1;
                }

                coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
                int x = coords[0];
                int yblock = coords[1];
                int z = coords[2];
                int t = coords[3];
                int y = yblock * nyg_pack;
                int xodd = (z + t + cb) & 1;
		int ytmp = y;
                if(dir == 0) {
                        int row = (fb == 1 ? 1-xodd : xodd); 
                        y += row;
			ytmp += 1-row;
                }

		QSU3V *giBase = &gi[t*Pxyz + z*Pxy + y*Vxh + x];
                fptype * __restrict rBuf = &res[6*pktsize*pkt];
		ff_pack_face_v_dir_vec_noinline(giBase, rBuf, dir*2+fb);
	}

}

void ff_unpack_face_v( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3V *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
	int pktsize = VECLEN;
        if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	int shifts[4] = {Vxh-1, Vxh*(Vy-1), Vxh*Vy*(Vz-1), Vxh*Vy*Vz*(Vt-1)}, shift;
	shift = shifts[dir] * (fb == 0 ? 1 : -1);

#ifdef MT_EXECUTION
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
#else
#if 1
#pragma omp parallel for
#endif
        for(int pkt = 0 ; pkt < npkts; pkt+=1) {
#endif
                int coords[4];
                int tmp = pkt;
                for(int j = 0; j < 4; j++) {
                        int tmp1 = tmp / lens[j];
                        coords[j] = tmp - tmp1*lens[j];
                        tmp = tmp1;
                }

                coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
                int x = coords[0];
                int yblock = coords[1];
                int z = coords[2];
                int t = coords[3];
                int y = yblock * nyg_pack;
                int xodd = (z + t + cb) & 1;

                if(dir == 0) {
                        int row = (fb == 1 ? 1-xodd : xodd);
                        y += row;
                }

		QSU3V *goBase = &res[t*Pxyz + z*Pxy + y*Vxh + x];
		fptype * __restrict rBuf = &gi[6*pktsize*pkt];
		ff_unpack_face_v_dir_vec_noinline(goBase+shift, rBuf, dir*2+fb);
	}
}

void ff_pack_and_send_boundaries_v(int tid, QSU3V *gi, int cb, int fb, int d)
{
//    MYASSERT(Vxh%2 == 0);
//    MYASSERT(Vx%4 == 0);
//    MYASSERT(Vy%4 == 0);
//    MYASSERT(Vz%4 == 0);
//    MYASSERT(Vt%4 == 0);
    if(!local_dir[d]) {
	if(tid == 0) {
		MYASSERT(MPI_Irecv((void *)vcommsBuf[4*d+fb], bufSizeV[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+(fb^1)], fb, MPI_COMM_THISJOB, &reqRecvsV[2*d+fb]) == MPI_SUCCESS);
	}
	ff_pack_face_v(tid, NCores, n_threads_per_core, gi, (fptype *)vcommsBuf[2+4*d+fb], cb, d, fb);
	if(tid == 0) {
		MYASSERT(MPI_Isend((void *)vcommsBuf[2+4*d+fb], bufSizeV[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], fb, MPI_COMM_THISJOB, &reqSendsV[2*d+fb]) == MPI_SUCCESS);
	}
    }
}

void ff_recv_and_unpack_boundaries_v(int tid, QSU3V *gio, int cb, int fb, int d)
{
    if(!local_dir[d]) {
	if(tid==0) {
		MYASSERT(MPI_Wait(&reqRecvsV[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
	}
	ff_unpack_face_v(tid, NCores, n_threads_per_core, (fptype *)vcommsBuf[4*d+fb], gio, cb, d, fb);
   }
   if(!local_dir[d]) {
	if(tid == 0) {
		MYASSERT(MPI_Wait(&reqSendsV[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
	}
   }
}

// face pack for 3rd neighbor. this code is much similar to 1st neighbor code above except that we gather more data.
// this is only invokved for vector data and not matrix. this is kind-of specialization of 1st neighbor code.
// validated for forward only since Fermion Force needs only forward comm/ 3rd neighbors.
void ff_pack_face_v3(int tid, int num_cores, int threads_per_core, QSU3V *gi, fptype *res, int cb, int dir, int fb)
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	int shifts[4] = {Vxh-1, Vxh*(Vy-1), Vxh*Vy*(Vz-1), Vxh*Vy*Vz*(Vt-1)}, shift;
        lens[dir] = 1;

        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pktsize = VECLEN;
	if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	MYASSERT(18*pktsize*npkts == bufSizeV3[dir]);
	shift = shifts[dir] * (fb == 0 ? 1 : -1);

#ifdef MT_EXECUTION
	int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;
        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
#else
#pragma omp parallel for
        for(int pkt = 0 ; pkt < npkts; pkt+=1) {
#endif
		int coords[4];
		int tmp = pkt;
		for(int j = 0; j < 4; j++) {
			int tmp1 = tmp / lens[j];
			coords[j] = tmp - tmp1*lens[j];
			tmp = tmp1;
		}

		coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
		int x = coords[0];
		int yblock = coords[1];
		int z = coords[2];
		int t = coords[3];
		int y = yblock * nyg_pack;
// xshift is an additional parameter used to offset into the lattice index correctly.
		int xshift=0;
		int xodd = (z + t + 1-cb) & 1;
	 	int row;
		if(dir == 0) {
			row = xodd; // picking y with single pack site first
			y += row; 
			xshift = (row == 0 ? Vxh : -Vxh);
		}

		QSU3V *giBase = &gi[t*Pxyz + z*Pxy + y*Vxh + x];
		fptype * __restrict rBuf = &res[18*pktsize*pkt];

// handle packing in each direction separately. this is needed for X direction specifically.
// Y, Z, T pack code can be combined but for brevity they are left separate.
		if(dir == 0) {
			ff_pack_face_v3_dir_vec_noinline(giBase, rBuf, dir*2+fb);
			ff_pack_face_v3_dir_vec_noinline(giBase+xshift, rBuf+6*pktsize, dir*2+fb);
			ff_pack_face_v3_dir_vec_noinline(giBase+xshift+1, rBuf+12*pktsize, dir*2+fb);
//			ff_pack_face_v3_dir_vec_noinline(giBase+xshift+shift, rBuf+12*pktsize, dir*2+fb);

		}
		else if(dir == 1){
			ff_pack_face_v_dir_vec_noinline(giBase, rBuf, dir*2+fb);
			ff_pack_face_v_dir_vec_noinline(giBase+Vxh, rBuf+6*pktsize, dir*2+fb);
			ff_pack_face_v_dir_vec_noinline(giBase+2*Vxh, rBuf+12*pktsize, dir*2+fb);
		}
		else if(dir == 2){
			ff_pack_face_v_dir_vec_noinline(giBase, rBuf, dir*2+fb);
			ff_pack_face_v_dir_vec_noinline(giBase+Vxh*Vy, rBuf+6*pktsize, dir*2+fb);
			ff_pack_face_v_dir_vec_noinline(giBase+2*Vxh*Vy, rBuf+12*pktsize, dir*2+fb);
		}
		else if(dir == 3){
			ff_pack_face_v_dir_vec_noinline(giBase, rBuf, dir*2+fb);
			ff_pack_face_v_dir_vec_noinline(giBase+Vxh*Vy*Vz, rBuf+6*pktsize, dir*2+fb);
			ff_pack_face_v_dir_vec_noinline(giBase+2*Vxh*Vy*Vz, rBuf+12*pktsize, dir*2+fb);
		}
	}
}

void ff_unpack_face_v3( int tid, int num_cores, int threads_per_core, fptype *gi, QSU3V *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
	int pktsize = VECLEN;
	int shifts[4] = {Vxh-1, Vxh*(Vy-1), Vxh*Vy*(Vz-1), Vxh*Vy*Vz*(Vt-1)}, shift;
        if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	shift = shifts[dir] * (fb == 0 ? 1 : -1);

#ifdef MT_EXECUTION
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
#else
#pragma omp parallel for
        for(int pkt = 0 ; pkt < npkts; pkt+=1) {
#endif
		int coords[4];
		int tmp = pkt;
		for(int j = 0; j < 4; j++) {
			int tmp1 = tmp / lens[j];
			coords[j] = tmp - tmp1*lens[j];
			tmp = tmp1;
		}

		coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
		int x = coords[0];
		int yblock = coords[1];
		int z = coords[2];
		int t = coords[3];
		int y = yblock * nyg_pack;
		int xshift=0;
		int xodd = (z + t + 1-cb) & 1;
		if(dir == 0) {
			int row = xodd; // picking y with single pack site first
			y += row; 
			xshift = (row == 0 ? Vxh : -Vxh); 
		}

		QSU3V *goBase = &res[t*Pxyz + z*Pxy + y*Vxh + x];
		fptype * __restrict rBuf = &gi[18*pktsize*pkt];

		if(dir == 0) {
			ff_unpack_face_v3_dir_vec_noinline(goBase+shift, rBuf, dir*2+fb);
			ff_unpack_face_v3_dir_vec_noinline(goBase+xshift+shift-1, rBuf+6*pktsize, dir*2+fb);
			ff_unpack_face_v3_dir_vec_noinline(goBase+xshift+shift, rBuf+12*pktsize, dir*2+fb);
//			ff_unpack_face_v3_dir_vec_noinline(goBase+xshift, rBuf+6*pktsize, dir*2+fb);
//			ff_unpack_face_v3_dir_vec_noinline(goBase+xshift+shift, rBuf+12*pktsize, dir*2+fb);
		}
		else if(dir == 1){
			ff_unpack_face_v_dir_vec_noinline(goBase+shift-2*Vxh, rBuf, dir*2+fb);
			ff_unpack_face_v_dir_vec_noinline(goBase+shift-Vxh, rBuf+6*pktsize, dir*2+fb);
			ff_unpack_face_v_dir_vec_noinline(goBase+shift, rBuf+12*pktsize, dir*2+fb);
		}
		else if(dir == 2){
			ff_unpack_face_v_dir_vec_noinline(goBase+shift-2*Vxh*Vy, rBuf, dir*2+fb);
			ff_unpack_face_v_dir_vec_noinline(goBase+shift-Vxh*Vy, rBuf+6*pktsize, dir*2+fb);
			ff_unpack_face_v_dir_vec_noinline(goBase+shift, rBuf+12*pktsize, dir*2+fb);
		}
		else if(dir == 3){
			ff_unpack_face_v_dir_vec_noinline(goBase+shift-2*Vxh*Vy*Vz, rBuf, dir*2+fb);
			ff_unpack_face_v_dir_vec_noinline(goBase+shift-Vxh*Vy*Vz, rBuf+6*pktsize, dir*2+fb);
			ff_unpack_face_v_dir_vec_noinline(goBase+shift, rBuf+12*pktsize, dir*2+fb);
		}
	}
}

void ff_pack_and_send_boundaries_v3(int tid, QSU3V *gi, int cb, int fb, int d)
{
    if(!local_dir[d]) {
	if(tid == 0) {
		MYASSERT(MPI_Irecv((void *)v3commsBuf[4*d+fb], bufSizeV3[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+(fb^1)], fb, MPI_COMM_THISJOB, &reqRecvsV[2*d+fb]) == MPI_SUCCESS);
	}
	ff_pack_face_v3(tid, NCores, n_threads_per_core, gi, (fptype *)v3commsBuf[2+4*d+fb], cb, d, fb);
	if(tid == 0) {
		MYASSERT(MPI_Isend((void *)v3commsBuf[2+4*d+fb], bufSizeV3[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], fb, MPI_COMM_THISJOB, &reqSendsV[2*d+fb]) == MPI_SUCCESS);
	}
    }
}

void ff_recv_and_unpack_boundaries_v3(int tid, QSU3V *gio, int cb, int fb, int d)
{
   if(!local_dir[d]) {
        if(tid==0) {
		MYASSERT(MPI_Wait(&reqRecvsV[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
        }
        ff_unpack_face_v3(tid, NCores, n_threads_per_core, (fptype *)v3commsBuf[4*d+fb], gio, cb, d, fb);
   }
   if(!local_dir[d]) {
	if(tid == 0) {
		MYASSERT(MPI_Wait(&reqSendsV[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
	}
   }
}

#endif /* ENABLE_MPI */

