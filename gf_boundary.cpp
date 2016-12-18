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

#include "ks_config.h"
#include "ks_globals.h"        /* All global values and macros */
#include "gf_boundary.h"
#include "qphix_internal.h"
#include "gauge_force_imp_complete_specialization.h"
#include "misc.h"

#if QPHIX_PrecisionInt == 1
#undef gf_setup_comms
#define gf_setup_comms gf_setup_comms_F
#undef gf_destroy_comms
#define gf_destroy_comms gf_destroy_comms_F
#elif QPHIX_PrecisionInt == 2
#undef gf_setup_comms
#define gf_setup_comms gf_setup_comms_D
#undef gf_destroy_comms
#define gf_destroy_comms gf_destroy_comms_D
#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif

extern int PadBound;
extern int PadNeigh;
extern char * BoundTable;
extern unsigned int * NeighTable;

#ifndef ENABLE_MPI
/* gauge & force remappings */
void pack_gauge_face_dir(int tid, Gauge *gi, int dir, int cb){}
void pack_gauge_face_dir(int tid, fptype *gi, int dir, int cb){}
void pack_momentum_face_dir(int tid, Hermit *hi, int dir, int cb){}
void pack_momentum_face_dir(int tid, fptype *hi, int dir, int cb){}
void unpack_gauge_face_dir(int tid, Gauge *go, int dir, int cb){}
void unpack_gauge_face_dir(int tid, fptype *go, int dir, int cb){}
void unpack_momentum_face_dir(int tid, Hermit *ho, int dir, int cb){}
void unpack_momentum_face_dir(int tid, fptype *ho, int dir, int cb){}

void gf_pack_and_send_boundaries_u(int tid, Gauge *gi, int cb){}
void gf_recv_and_unpack_boundaries_h(int tid, Hermit *hio, HermitHelperYZT *htmp, Gauge *gio, int cb){}
void gf_recv_and_unpack_u_and_send_boundaries_h(int tid, Hermit *hio, HermitHelperYZT *htmp, Gauge *gio, int cb){}
#if QPHIX_PrecisionInt == 1
void gf_print_boundary_timings(unsigned long long t_gf){}
#endif
void gf_setup_comms(){}
void gf_destroy_comms(){}
#else // ENABLE_MPI is defined
#ifndef ASSUME_MULTINODE
#include <mpi.h>
#endif

#define GF_U_MPI_SEND_TAG(x) (102+(x))
#define GF_U_MPI_RECV_TAG(x) (103-(x))
#define GF_H_MPI_SEND_TAG(x) (104+(x))
#define GF_H_MPI_RECV_TAG(x) (105-(x))

const int nGX = (VECLEN < 2 ? 1 : 2);
const int nGY = (VECLEN < 4 ? 1 : 2);
const int nGZ = (VECLEN < 8 ? 1 : 2);
const int nGT = (VECLEN < 16 ? 1 : 2);

extern int neigh_ranks[8];
extern int myCoord[4];
extern int Ls[4];
static int bufSizeG[4], bufSizeG2[4];
static int bufSizeH[4], bufSizeH2[4], bufSizeHtmp[3];
#if QPHIX_PrecisionInt == 2
fptype *gcommsBuf[16], *hcommsBuf[16];
fptype *gcommsBuf2[16], *hcommsBuf2[16];
fptype *hlocalBuf[6], *glocalBuf[8];
#else
#if QPHIX_PrecisionInt == 1
extern fptype *gcommsBuf[16], *hcommsBuf[16];
extern fptype *gcommsBuf2[16], *hcommsBuf2[16];
extern fptype *hlocalBuf[6], *glocalBuf[8];
#endif
#endif

static unsigned int PackMask[2][4] = {{0,0,0,0},{0,0,0,0}};

#ifndef ASSUME_MULTINODE
static MPI_Request reqSends[16];
static MPI_Request reqRecvs[16];
static MPI_Request reqSendsU[8];
static MPI_Request reqSendsH[8];
static MPI_Request reqRecvsU[8];
static MPI_Request reqRecvsH[8];
#endif

static int nCommsDirs = 0;
static unsigned long send_lat[8];
static unsigned long recv_lat[8];
static unsigned int send_order[8];
static unsigned int recv_order[8];

static unsigned long long tt_pack[255][8];
static unsigned long long tt_unpack[255][8];
static unsigned long long t_boundary[6][8];

void gf_setup_comms()
{
        int gDims[4] = {Gx, Gy, Gz, Gt};
        int lDims[4] = {Gx, Gy, Gz, Gt};
	for(int i = 0; i < 4; i++) {
	    lDims[i] = gDims[i] / geometry[i];
	    for(int j=0; j<((1<<i)<VECLEN ? VECLEN/2 : VECLEN); ++j)
	    {
		PackMask[0][i] ^= (1<<(j%(1<<i)+(j/(1<<i))*((2<<i))));
	    }
	    PackMask[1][i] = ~PackMask[0][i];
	    if((PackMask[1][i]<<(sizeof(unsigned int)*8-VECLEN))==0) PackMask[1][i] = PackMask[0][i];
	}

        int numProcs = geometry[0] * geometry[1] * geometry[2] * geometry[3];
        MYASSERT(numProcs == nRanks);
        for(int i = 0; i < 4; i++) {
                bufSizeG[i] = 1; bufSizeH[i] = 1;
		bufSizeG2[i] = 1; bufSizeH2[i] = 1;
		if(i<3) bufSizeHtmp[i] = 1;
                for(int j = 0; j < 4; j++) if(j != i) {
		    bufSizeG[i] *= lDims[j]; 
		    bufSizeH[i] *= lDims[j]; 
		    bufSizeG2[i] *= lDims[j]; 
		    bufSizeH2[i] *= lDims[j]; 
		    if(i!=0) bufSizeHtmp[i-1] *= lDims[j]; 
		}
		bufSizeG[i] *= 7*12/2 /* Gauge */; 
		bufSizeH[i] *= 7*8/2 /* Momentum */;
		bufSizeG2[i] *= 12/2 /* SU(3) matrix */; 
		bufSizeH2[i] *= 8/2 /* Hermitian */;
		if(i<3)
		{
		    bufSizeHtmp[i] *= 7*8 /* Temporary momentum storage in + & - x direction */;
		    if((1<<(i+1))<VECLEN) bufSizeHtmp[i] *= 2;
		}
	}

        int alignedBufSizeG[4] = {0}, alignedBufSizeG2[4] = {0};
	int alignedBufSizeH[4] = {0}, alignedBufSizeH2[4] = {0};
	int alignedBufSizeHtmp[3] = {0};
        int totalBufSizeG = 0, totalBufSizeG2 = 0; 
	int totalBufSizeH = 0, totalBufSizeH2 = 0;
	int totalBufSizeHtmp = 0;
        int pgSize = 4096 / sizeof(fptype);
        for(int i = 0; i < 4; i++) {
                alignedBufSizeG[i] = (bufSizeG[i] + (pgSize-1)) & ~(pgSize-1);
		alignedBufSizeH[i] = (bufSizeH[i] + (pgSize-1)) & ~(pgSize-1);
		alignedBufSizeG2[i] = (bufSizeG2[i] + (pgSize-1)) & ~(pgSize-1);
		alignedBufSizeH2[i] = (bufSizeH2[i] + (pgSize-1)) & ~(pgSize-1);
		if(i>0) alignedBufSizeHtmp[i-1] = (bufSizeHtmp[i-1] + (pgSize-1)) & ~(pgSize-1);
                if(!local_dir[i]) {
                        totalBufSizeG += alignedBufSizeG[i];
			totalBufSizeH += alignedBufSizeH[i];
			totalBufSizeG2 += alignedBufSizeG2[i];
			totalBufSizeH2 += alignedBufSizeH2[i];
			if(i>0) totalBufSizeHtmp += alignedBufSizeHtmp[i-1];
                }
        }
#if QPHIX_PrecisionInt == 2
        fptype *htmpBuf = NULL, *gtmpBuf = NULL;
	fptype *htmpBuf2 = NULL, *gtmpBuf2 = NULL;
        if(totalBufSizeG > 0) {
		gtmpBuf = (fptype*)MALLOC(6 * totalBufSizeG * sizeof(fptype), 4096);
		MYASSERT(gtmpBuf != NULL);
		gtmpBuf2 = (fptype*)MALLOC(4 * totalBufSizeG2 * sizeof(fptype), 4096);
		MYASSERT(gtmpBuf2 != NULL);
	}
	if(totalBufSizeH > 0) {
#ifndef USE_WAITANY
                htmpBuf = (fptype*)MALLOC((4 * totalBufSizeH + totalBufSizeHtmp * 2) * sizeof(fptype), 4096);
#else
		/* H & G shares the same commsBuf */
		htmpBuf = (fptype*)MALLOC((4 * totalBufSizeH + totalBufSizeHtmp * 2) * sizeof(fptype), 4096);
#endif
                MYASSERT(htmpBuf != NULL);
		htmpBuf2 = (fptype*)MALLOC(4 * totalBufSizeH2 * sizeof(fptype), 4096);
		MYASSERT(htmpBuf2 != NULL);
        }
        for(int i = 0; i < 4; i++) {
                if(!local_dir[i]) {
                        nCommsDirs += 2;
#ifndef ASSUME_MULTINODE
                        for(int j = 0; j < 4; j++) {
                                gcommsBuf[4*i+j] = gtmpBuf;
				gcommsBuf2[2*i+(j&1)+(j/2)*8] = gtmpBuf2;
				hcommsBuf2[2*i+(j&1)+(j/2)*8] = htmpBuf2;
#ifdef USE_WAITANY
				hcommsBuf[4*i+j] = gtmpBuf;
#else
				hcommsBuf[4*i+j] = htmpBuf;
				htmpBuf += alignedBufSizeH[i];
#endif
				gtmpBuf += alignedBufSizeG[i];
				gtmpBuf2 += alignedBufSizeG2[i];
				htmpBuf2 += alignedBufSizeH2[i];
                        }
			if(i>0)	
			{
                            hlocalBuf[2*i-2] = htmpBuf;
                            htmpBuf += alignedBufSizeHtmp[i-1];
                            hlocalBuf[2*i-1] = htmpBuf;
                            htmpBuf += alignedBufSizeHtmp[i-1];
			}
			glocalBuf[2*i+0] = gtmpBuf;
			gtmpBuf += alignedBufSizeG[i];
			glocalBuf[2*i+1] = gtmpBuf;
			gtmpBuf += alignedBufSizeG[i];
#else
                        for(int j = 0; j < 2; j++) {
				gcommsBuf[4*i+j] = gtmpBuf;
				gcommsBuf[4*i+3-j] = gtmpBuf;
#ifdef USE_WAITANY
				hcommsBuf[4*i+j] = gtmpBuf;
				hcommsBuf[4*i+3-j] = gtmpBuf;
#else
				hcommsBuf[4*i+j] = htmpBuf;
				hcommsBuf[4*i+3-j] = htmpBuf;
				htmpBuf += 2*alignedBufSizeH[i];
#endif
				gtmpBuf += 2*alignedBufSizeG[i];
                        }
#endif
                }
                else {
                        gcommsBuf[4*i] = 0x00;
                        gcommsBuf[4*i+1] = 0x00;
                        gcommsBuf[4*i+2] = 0x00;
                        gcommsBuf[4*i+3] = 0x00;
			hcommsBuf[4*i] = 0x00;
			hcommsBuf[4*i+1] = 0x00;
			hcommsBuf[4*i+2] = 0x00;
			hcommsBuf[4*i+3] = 0x00;
                        if(i>0) 
			{
			    hlocalBuf[2*i-2] = 0x00;
                            hlocalBuf[2*i-1] = 0x00;
                        }
			glocalBuf[2*i+0] = 0x00;
                        glocalBuf[2*i+1] = 0x00;
			gcommsBuf2[2*i] = 0x00;
			gcommsBuf2[2*i+1] = 0x00;
			gcommsBuf2[2*i+8] = 0x00;
			gcommsBuf2[2*i+9] = 0x00;
			hcommsBuf2[2*i] = 0x00;
			hcommsBuf2[2*i+1] = 0x00;
			hcommsBuf2[2*i+8] = 0x00;
			hcommsBuf2[2*i+9] = 0x00;
                }
        }
#endif
#ifndef ASSUME_MULTINODE
        for(int d = 0; d < 16; d++) {
                reqSends[d] = reqRecvs[d] = MPI_REQUEST_NULL;
        }
	for(int d = 0; d < 8; d++) {
		reqSendsU[d] = reqSendsH[d] = reqRecvsU[d] = reqRecvsH[d] = MPI_REQUEST_NULL;
	}
#endif
        for(int d = 0; d < 8; d++) {
                for(int t = 0; t < 255; t++) {
                        tt_pack[t][d] = tt_unpack[t][d] = 0;
                }
                for(int t = 0; t < 6; t++) {
                        t_boundary[t][d] = 0;
                }
        }
        if(nParallelSends > nCommsDirs) nParallelSends = nCommsDirs;
        if(nParallelSends < 1) nParallelSends = 1;
        while(nCommsDirs % nParallelSends != 0) nParallelSends--;

        if(nParallelRecvs < 1) nParallelRecvs = 1;
        else if(nParallelRecvs > 2) nParallelRecvs = 2;

#ifdef USE_CML_TG
        setup_cml_thread_groups();
#endif
}

#if QPHIX_PrecisionInt == 1

void gf_setup_boundaries(void)
{
}

#endif

static void gf_pack_face_u(int tid, int num_cores, int threads_per_core, Gauge *gi, fptype *res, int cb, int dir, int fb)
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;

        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pktsize = VECLEN;
	if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	/* buffer size : 7*12*pktsize */
	MYASSERT(84*pktsize*npkts == bufSizeG[dir]);

        const unsigned int *NTNow = NeighTable+cb*PadNeigh;

	int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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
                int xodd = (z + t + 1-cb) & 1;
		int ytmp = y;
                if(dir == 0) {
                        int row = (fb == 1 ? 1-xodd : xodd); 
                        y += row;
			ytmp += 1-row;
                }
#if 0
#pragma omp critical
		{
		printf("fb = %d , thread %d : pkt = %d (%d, %d, %d, %d) ind = %d\n", fb, tid, pkt, x, y, z, t, t*Pxyz + z*Pxy + y*Vxh + x);
		}
#endif
                int pkt_next = pkt + threads_per_core < high_pkt ? pkt + threads_per_core : low_pkt + smtid;

                tmp = pkt_next;
                for(int j = 0; j < 4; j++) {
                        int tmp1 = tmp / lens[j];
                        coords[j] = tmp - tmp1*lens[j];
                        tmp = tmp1;
                }

                coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
                int x_next = coords[0];
                int yblock_next = coords[1];
                int z_next = coords[2];
                int t_next = coords[3];
                int y_next = yblock_next * nyg_pack;
                int xodd_next = (z_next + t_next + 1-cb) & 1;

                if(dir == 0) y_next +=(fb == 1 ? 1-xodd_next : xodd_next);
		Gauge *giBase = &gi[t*Pxyz + z*Pxy + y*Vxh + x];
		//MYASSERT(t*Pxyz + z*Pxy + y*Vxh + x==NTNow[16*(t*Pxyz + z*Pxy + ytmp*Vxh + x)+fb+2*dir]);
		int off_next = (t_next-t)*Pxyz+(z_next-z)*Pxy+(y_next-y)*Vxh+(x_next-x);

                fptype * __restrict rBuf = &res[84*pktsize*pkt];
                fptype * __restrict lBuf = &glocalBuf[2*dir+1-fb][84*pktsize*pkt];
//#pragma omp critical
		//if(fb==1) {
		//printf("thread %d : dir = %d\n", tid, dir*2+fb);
		gauge_pack_face_dir_vec_noinline(giBase, lBuf, rBuf, dir*2+fb);
		//}
	}

}

/* pack gauge links along one direction */
static void gf_pack_face_u_dir( int tid, int num_cores, int threads_per_core, Gauge *gi, fptype *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;

        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pktsize = VECLEN;
        if((1 << dir) < VECLEN) pktsize = VECLEN/2;
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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

                int pkt_next = pkt + threads_per_core < high_pkt ? pkt + threads_per_core : low_pkt + smtid;

                tmp = pkt_next;
                for(int j = 0; j < 4; j++) {
                        int tmp1 = tmp / lens[j];
                        coords[j] = tmp - tmp1*lens[j];
                        tmp = tmp1;
                }

                coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
                int x_next = coords[0];
                int yblock_next = coords[1];
                int z_next = coords[2];
                int t_next = coords[3];
                int y_next = yblock_next * nyg_pack;
                int xodd_next = (z_next + t_next + cb) & 1;

                if(dir == 0) y_next +=(fb == 1 ? 1-xodd_next : xodd_next);
                Gauge *giBase = &gi[t*Pxyz + z*Pxy + y*Vxh + x];
                int off_next = (t_next-t)*Pxyz+(z_next-z)*Pxy+(y_next-y)*Vxh+(x_next-x);

                fptype * __restrict rBuf = &res[12*pktsize*pkt];

		gauge_pack_face_dir_vec_noinline(giBase, rBuf, dir*2+fb);
	}
}

static void gf_pack_face_u_dir( int tid, int num_cores, int threads_per_core, fptype *gi, fptype *res, int cb, int dir, int fb )
{
	int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
	int lens[4] = {Nxh, Ny/nyg_pack, Nz, Nt};
	int lens1[4] = {Nxh, Ny/nyg_pack, Nz, Nt};
	lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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
		for(int i=0; i<12; ++i)
			res[12*pkt+i] = gi[(((t*Nz+z)*Ny+y)*Nx+x)*72+18*dir+i];
	}
}

/* pack gauge momentum along one direction */
static void gf_pack_face_h_dir( int tid, int num_cores, int threads_per_core, Hermit *hi, fptype *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;

        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pktsize = VECLEN;
        if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	MYASSERT(8*pktsize*npkts == bufSizeH2[dir]);

        int pkts_per_core = (npkts + num_cores - 1) / num_cores;
        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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

		int pkt_next = pkt + threads_per_core < high_pkt ? pkt + threads_per_core : low_pkt + smtid;
                tmp = pkt_next;
                for(int j = 0; j < 4; j++) {
                        int tmp1 = tmp / lens[j];
                        coords[j] = tmp - tmp1*lens[j];
                        tmp = tmp1;
                }

                int x_next = coords[0];
                int yblock_next = coords[1];
                int z_next = coords[2];
                int t_next = coords[3];
                int y_next = yblock_next * nyg_pack;
                int xodd_next = (z_next + t_next + cb) & 1;

                if(dir == 0) y_next +=(fb == 1 ? 1-xodd_next : xodd_next);
                Hermit *hiBase = &hi[t*Pxyz + z*Pxy + y*Vxh + x];
                int off_next = (t_next-t)*Pxyz+(z_next-z)*Pxy+(y_next-y)*Vxh+(x_next-x);

                fptype * __restrict rBuf = &res[8*pktsize*pkt];
#if 0
		for(int k=0; k<8; ++k) for(int v=0; v<pktsize; ++v)
		{
			printf("thread %d : rBuf[%d][%d][%d] = %f\n", tid, pkt, k, v, res[8*pktsize*pkt+k*pktsize+v]);
			fflush(stdout);
		}		
		printf("thread %d : x= %d y= %d z= %d t= %d pkt= %d\n", tid, x, y, z, t, pkt);
		fflush(stdout);
#endif
                momentum_pack_face_dir_vec_noinline(hiBase, rBuf, dir*2+fb);
        }
}

static void gf_pack_face_h_dir( int tid, int num_cores, int threads_per_core, fptype *hi, fptype *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens[4] = {Nxh, Ny/nyg_pack, Nz, Nt};
        int lens1[4] = {Nxh, Ny/nyg_pack, Nz, Nt};
        lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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
                for(int i=0; i<8; ++i)
		{
			int c1 = i/4;
			int c2=2;
			int ir = i&1;
			if(i<2 || i==7) c2=1;
			else if(i==6) { c1=c2=0; ir=1; }
                        res[8*pkt+i] = hi[(((t*Nz+z)*Ny+y)*Nx+x)*72+18*dir+6*c1+2*c2+ir];
		}
        }
}

/* Unpack guage field and update boundary dH in buffer */
/* cb == 1 : EVEN */
static void gf_unpack_face_u(int tid, int num_cores, int threads_per_core, 
				fptype *gibuf, fptype *gtmpbuf[8], 
				Gauge *gio, Hermit *hio, HermitHelperYZT *htmp, 
				fptype kappaS, fptype kappaR, fptype kappaB,
				fptype epsilonH,
				char accumulate, int cb, int dir, int fb)
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;
	int npkts = lens[0] * lens[1] * lens[2] * lens[3];

#if 0
	int lens0[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	lens0[dir] = 1;
	for(int i=0; i<d; ++i)
	    if((accumulate & (1<<(2*i+fb)))==0) 
		lens0[i]--;
        int npkts0 = lens0[0] * lens0[1] * lens0[2] * lens0[3];
	/* No data to work on */
	if(npkts0==0) return;
#endif
        int pktsize = VECLEN;
	if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	int hpktsize[4] = {8*VECLEN, 8*VECLEN, 8*VECLEN, 8*VECLEN};
	for(int i=0; i<4; ++i) if((1<<i) < VECLEN) hpktsize[i] = 4*VECLEN;

	MYASSERT(84*pktsize*npkts == bufSizeG[dir]);
	MYASSERT(56*pktsize*npkts == bufSizeH[dir]);
	if(dir!=0) MYASSERT(112*npkts*VECLEN == bufSizeHtmp[dir-1]);

        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;
	const unsigned int *NTNow = NeighTable+cb*PadNeigh;
	const char *BTNow = BoundTable+cb*PadBound;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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

		int pkth[2] = {0};
		if(dir!=0) for(int i=0; i<2; ++i)
		{
		    int coordsh[4] = {x, yblock, z, t};
		    int xoddh = (xodd + y) & 1;
		    coordsh[0] = (coordsh[0] + (i==0 ? (-1+xoddh) : xoddh) + Vxh)%Vxh;
		    for(int j=3; j>=0; --j) if(j!=dir) 
			pkth[i] = pkth[i]*lens[j] + coordsh[j];
		}
		
                int pkt_next = pkt + threads_per_core < high_pkt ? pkt + threads_per_core : low_pkt + smtid;

                tmp = pkt_next;
                for(int j = 0; j < 4; j++) {
                        int tmp1 = tmp / lens[j];
                        coords[j] = tmp - tmp1*lens[j];
                        tmp = tmp1;
                }

                coords[dir] = (fb == 0 ? 0 : lens1[dir] - 1);
                int x_next = coords[0];
                int yblock_next = coords[1];
                int z_next = coords[2];
                int t_next = coords[3];
                int y_next = yblock_next * nyg_pack;
                int xodd_next = (z_next + t_next + cb) & 1;
                if(dir == 0) y_next +=(fb == 1 ? 1-xodd_next : xodd_next);

		int ind = t*Pxyz + z*Pxy + y*Vxh + x;
		const char isBoundary = BTNow[4*ind+2];
		int flag = 1;
		for(int i=0; i<8; ++i) if(i!=dir*2+fb && ((isBoundary>>i) & 1) && !((accumulate>>i) & 1)) { flag = 0; break; }
		if(flag==0) continue;

		Gauge *giBase[7];
                fptype *rBuf[8];
                fptype *lBuf[8];
                fptype *hBuf[8];
		fptype *hlBuf[6];
		Hermit *hoBase[8];
		int len = 0;
		for(int i=0, n=0; i<8; ++i) {
			/* Buffers to read in and write to */
			if(i!=2*dir+fb) {
				if(!((BTNow[4*ind]>>i)&1)) { /* Read from buffers */
				    int nygtmp=1;
				    if(i<2) nygtmp=2;
				    int lentmp[4] = {Vxh, Vy/nygtmp, Vz, Vt};
				    int coordtmp[4] = {x, y/nygtmp, z, t};
				    int ii=0;
				    for(int itmp=3; itmp>=0; --itmp) if(itmp!=(i>>1))
					ii = coordtmp[itmp] + lentmp[itmp]*ii;
				    rBuf[i] = &gtmpbuf[i][84*hpktsize[i>>1]/8*ii];
			 	    lBuf[i] = &glocalBuf[i][84*hpktsize[i>>1]/8*ii];
				    hBuf[i] = &hcommsBuf[2+2*(i>>1)+i][7*hpktsize[i>>1]*ii];
				}
				else {
				    rBuf[i] = 0x00;
				    lBuf[i] = 0x00;
				    hBuf[i] = 0x00;
				}
				giBase[n] = &gio[NTNow[16*ind+i]];
				n++;
			}
			else {
			    rBuf[i] = &gibuf[84*pktsize*pkt];
			    lBuf[i] = &glocalBuf[2*dir+fb][84*pktsize*pkt];
			    hBuf[i] = &hcommsBuf[2+4*dir+fb][56*pktsize*pkt];
			    len++;
			}
			if(i<2)
			{
			    	if(dir==0) hoBase[i] = &hio[NTNow[16*ind+i]];
			    	else hoBase[i] = (Hermit *)&hlocalBuf[2*(dir-1)+fb][(112*pkth[i]+56*(1-i))*VECLEN];
			}
			else
			    	hoBase[i] = (Hermit *)&htmp[NTNow[16*ind+i]][(i&1)? i-3 : i-1];
		}
		//fptype * __restrict rBuf = &gibuf[84*pktsize*pkt];
		//fptype * __restrict lBuf = &glocalBuf[2*dir+fb][84*pktsize*pkt];
		//fptype * __restrict hBuf = &hcommsBuf[2+4*dir+fb][56*pktsize*pkt];
		//rBuf = {rBuf0, rBuf1, rBuf2, rBuf3, rBuf4, rBuf5, rBuf6, rBuf7};
		//lBuf = { lBuf0, lBuf1, lBuf2, lBuf3, lBuf4, lBuf5, lBuf6, lBuf7 };
		//hBuf = { hBuf0, hBuf1, hBuf2, hBuf3, hBuf4, hBuf5, hBuf6, hBuf7 };
//#pragma omp critical
		{
//		printf("node %d thread %d : calling gauge_unpack_face_pack_gforce_dir_vec_noinline\n", myRank, tid);
		gauge_unpack_face_pack_gforce_dir_vec_noinline( 
			hBuf, /*hlBuf,*/ rBuf, lBuf,
			hoBase,
			giBase,
			dir*2+fb, hpktsize, /* 12*pktsize, 8*pktsize, */
			kappaS, kappaR, kappaB, 
			epsilonH,
			isBoundary
		);
		}
	}
}

static void gf_unpack_face_h(int tid, int num_cores, int threads_per_core, fptype *hbuf, Hermit *hio, HermitHelperYZT *htmp, int cb, int dir, int fb, char accx)
{
	int nyg_pack = 1;
	if(dir == 0) nyg_pack = 2;
	int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	lens[dir] = 1;

	int npkts = lens[0] * lens[1] * lens[2] * lens[3];
	int pktsize = VECLEN;
	if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	int pkts_per_core = (npkts + num_cores - 1) / num_cores;
	MYASSERT(56*pktsize*npkts == bufSizeH[dir]);
	if(dir>0) MYASSERT(112*npkts*VECLEN == bufSizeHtmp[dir-1]);

	int cid = tid/threads_per_core;
	int smtid = tid - threads_per_core * cid;
	int low_pkt = cid*pkts_per_core;
	int high_pkt = (cid+1)*pkts_per_core;
	if ( high_pkt > npkts ) high_pkt = npkts;
	const char *BTNow = BoundTable + (1-cb)*PadBound;

	for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
		int coords[4], tmp = pkt;
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
                int xodd = (z + t + 1-cb) & 1;

                if(dir == 0) {
                        int row = (fb == 1 ? 1-xodd : xodd);
                        y += row;
                }

		int ind = t*Pxyz+z*Pxy+y*Vxh+x;

		fptype * __restrict rBuf = &hbuf[56*pkt*pktsize];
		Hermit *hiBase[2] = {0x00, 0x00};
		const char accumulate = BTNow[4*ind];
		if(dir!=0)
		{
		    hiBase[0] = (Hermit*)&hlocalBuf[2*(dir-1)+fb][112*pkt*VECLEN];
		    hiBase[1] = (Hermit*)&hlocalBuf[2*(dir-1)+fb][(112*pkt+56)*VECLEN];
		    for(int j=2; j<8; ++j)
		    {
			if((j>>1)!=dir && !(accumulate&(1<<j)) && !(accx&(1<<j)))
			{
			    hiBase[0] = 0x00;
			    hiBase[1] = 0x00;
			    break;
			}
		    }
		}
		gforce_unpack_face_dir_vec_noinline( rBuf, &hio[ind], hiBase, &htmp[ind], 2*dir+fb, 8*pktsize, !(accumulate&1), !(accumulate&2) );
	}
}

static void gf_unpack_face_u_dir( int tid, int num_cores, int threads_per_core, fptype *gi, Gauge *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
	int pktsize = VECLEN;
        if((1 << dir) < VECLEN) pktsize = VECLEN/2;
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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
                int y0 = y;
                int xodd = (z + t + 1-cb) & 1;

                if(dir == 0) {
                        int row = (fb == 1 ? 1-xodd : xodd);
                        y += row;
                        y0 += 1-row;
                }

		Gauge *goBase = &res[t*Pxyz + z*Pxy + y0*Vxh + x];
		fptype * __restrict rBuf = &gi[12*pktsize*pkt];

		gauge_unpack_face_dir_vec_noinline(goBase, rBuf, dir*2+fb);
	}
}

static void gf_unpack_face_u_dir( int tid, int num_cores, int threads_per_core, fptype *gi, fptype *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
	int lens1[4] = {Nxh, Ny/nyg_pack, Nz, Nt};
	int lens[4] = {Nxh, Ny/nyg_pack, Nz, Nt};
	lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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

		for(int i=0; i<12; ++i)
			res[(((t*Nz+z)*Ny+y)*Nx+x)*72+18*dir+i] = gi[12*pkt+i];
	}
}

static void gf_unpack_face_h_dir( int tid, int num_cores, int threads_per_core, fptype *hi, Hermit *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
        int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
        int pktsize = VECLEN;
        if((1 << dir) < VECLEN) pktsize = VECLEN/2;
	MYASSERT(8*pktsize*npkts == bufSizeH2[dir]);
        int pkts_per_core = (npkts + num_cores - 1) / num_cores;

        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) {
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
                int y0 = y;
                int xodd = (z + t + 1-cb) & 1;

                if(dir == 0) {
                        int row = (fb == 1 ? 1-xodd : xodd);
                        y += row;
                        y0 += 1-row;
                }

                Hermit *hoBase = &res[t*Pxyz + z*Pxy + y0*Vxh + x];
                fptype * __restrict rBuf = &hi[8*pktsize*pkt];

                momentum_unpack_face_dir_vec_noinline(hoBase, rBuf, dir*2+fb);
        }
}

static void gf_unpack_face_h_dir( int tid, int num_cores, int threads_per_core, fptype *hi, fptype *res, int cb, int dir, int fb )
{
        int nyg_pack = 1;
        if(dir == 0) nyg_pack = 2;
	int pktsize = VECLEN;
	if((1<<dir) < VECLEN) pktsize = VECLEN/2; 

	int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
        lens[dir] = 1;
        int npkts = lens[0] * lens[1] * lens[2] * lens[3];
	MYASSERT(8*pktsize*npkts == bufSizeH2[dir]);

        int pkts_per_core = (npkts + num_cores - 1) / num_cores;
        int cid = tid/threads_per_core;
        int smtid = tid - threads_per_core * cid;

        int low_pkt = cid*pkts_per_core;
        int high_pkt = (cid+1)*pkts_per_core;
        if ( high_pkt > npkts ) high_pkt = npkts;

        for(int pkt = low_pkt + smtid; pkt < high_pkt; pkt+=threads_per_core) 
	{
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
		int indv = 0;
		for(int v=0; v<VECLEN; ++v) if((PackMask[fb][dir]>>v)&1)
		{
		    int y0 = v/nGX;
		    int x0 = v - y0*nGX;
		    int z0 = y0/nGY;
		    y0 -= z0*nGY;
		    int t0 = z0/nGZ;
		    z0 -= t0*nGZ;
		    x0 = Vxh*x0+x;
		    y0 = Vy*y0+y;
		    z0 = Vz*z0+z;
		    t0 = Vt*t0+t;
		    int ind0 = ((t0*Nz+z0)*Ny+y0)*Nxh+x0;

                    for(int i=0; i<8; ++i)
		    {
                        int c1 = i/4;
                        int c2=2;
                        int ir = i&1;
                        if(i<2 || i==7) {
			    c2=1;
			    if(i==7) ir = 0;
			}
                        else if(i==6) { c1=c2=0; }
			ir = 1-ir;
                        res[ind0*72+18*dir+6*c1+2*c2+ir] = (-1.+2.*ir)*hi[(8*pkt+i)*pktsize+indv];
			res[ind0*72+18*dir+6*c2+2*c1+ir] = hi[(8*pkt+i)*pktsize+indv];
#if 0
#pragma omp critical
			{
			if(res[ind0*72+18*dir+6*c2+2*c1+ir]!=0.)
				printf("rank %d thread %d : res[%d][%d][%d] = %f\n", myRank, tid, ind0, dir, 6*c2+2*c1+ir, res[ind0*72+18*dir+6*c2+2*c1+ir]);
			}
#endif
		    }
		    for(int i=0; i<3; ++i)
			res[ind0*72+18*dir+i*8] = 0.0;
		    res[ind0*72+18*dir+17] = -res[ind0*72+18*dir+1] - res[ind0*72+18*dir+9];
		    indv++;
		}
        }
}

#ifndef USE_CML_TG

void pack_gauge_face_dir(int tid, Gauge *gi, int dir, int cb)
{
        unsigned long long t0, t1;
	int nv[4] = {Nx, Ny, Nz, Nt};
	long bufsize = 6*Nx*Ny*Nz*Nt/nv[dir>>1];
        if(tid == 0) {
                t0 = __rdtsc();
                MYASSERT(MPI_Irecv((void *)gcommsBuf2[QPHIX_OPP_DIR(dir)], bufsize*sizeof(fptype), MPI_BYTE, neigh_ranks[QPHIX_OPP_DIR(dir)], GF_U_MPI_RECV_TAG((dir&1)^1), MPI_COMM_WORLD, &reqRecvsU[QPHIX_OPP_DIR(dir)]) == MPI_SUCCESS);
                t1 = __rdtsc();
                t_boundary[3][QPHIX_OPP_DIR(dir)] += (t1 - t0);
        }
        t0 = __rdtsc();
        gf_pack_face_u_dir(tid, NCores, n_threads_per_core, gi, gcommsBuf2[8+dir], cb, dir>>1, dir&1);
        t1 = __rdtsc();
        tt_pack[tid][dir] += (t1 - t0);
        gBar->wait(tid);
//#pragma omp barrier
        if(tid == 0) {
		t1 = __rdtsc();
                t_boundary[0][dir] += (t1 - t0);
#ifndef ASSUME_MULTINODE
                MYASSERT(MPI_Isend((void *)gcommsBuf2[8+dir], bufsize*sizeof(fptype), MPI_BYTE, neigh_ranks[dir], GF_U_MPI_SEND_TAG(dir&1), MPI_COMM_WORLD, &reqSendsU[dir]) == MPI_SUCCESS);
#endif
                t0 = __rdtsc();
                t_boundary[2][dir] += (t0 - t1);
        }

}

void pack_gauge_face_dir(int tid, fptype *gi, int dir, int cb)
{
        unsigned long long t0, t1;
        int nv[4] = {Nx, Ny, Nz, Nt};
        long bufsize = 6*Nx*Ny*Nz*Nt/nv[dir>>1];
        if(tid == 0) {
                t0 = __rdtsc();
                MYASSERT(MPI_Irecv((void *)gcommsBuf2[QPHIX_OPP_DIR(dir)], bufsize*sizeof(fptype), MPI_BYTE, neigh_ranks[QPHIX_OPP_DIR(dir)], GF_U_MPI_RECV_TAG((dir&1)^1), MPI_COMM_WORLD, &reqRecvsU[QPHIX_OPP_DIR(dir)]) == MPI_SUCCESS);
                t1 = __rdtsc();
                t_boundary[3][QPHIX_OPP_DIR(dir)] += (t1 - t0);
        }
        t0 = __rdtsc();
        gf_pack_face_u_dir(tid, NCores, n_threads_per_core, gi, gcommsBuf2[8+dir], cb, dir>>1, dir&1);
        t1 = __rdtsc();
        tt_pack[tid][dir] += (t1 - t0);
        gBar->wait(tid);
//#pragma omp barrier
        if(tid == 0) {
                t1 = __rdtsc();
                t_boundary[0][dir] += (t1 - t0);
#ifndef ASSUME_MULTINODE
                MYASSERT(MPI_Isend((void *)gcommsBuf2[8+dir], bufsize*sizeof(fptype), MPI_BYTE, neigh_ranks[dir], GF_U_MPI_SEND_TAG(dir&1), MPI_COMM_WORLD, &reqSendsU[dir]) == MPI_SUCCESS);
#endif
                t0 = __rdtsc();
                t_boundary[2][dir] += (t0 - t1);
        }
}

void pack_momentum_face_dir(int tid, Hermit *hi, int dir, int cb)
{
        unsigned long long t0, t1;
        int nv[4] = {Nx, Ny, Nz, Nt};
        long bufsize = 4*Nx*Ny*Nz*Nt/nv[dir>>1];
        if(tid == 0) {
                t0 = __rdtsc();
		//printf("rank %d : start receiving momentum from node %d... d = %d\n", myRank, neigh_ranks[QPHIX_OPP_DIR(dir)], QPHIX_OPP_DIR(dir));
		MYASSERT(MPI_Irecv((void *)hcommsBuf2[QPHIX_OPP_DIR(dir)], bufsize*sizeof(fptype), MPI_BYTE, neigh_ranks[QPHIX_OPP_DIR(dir)], GF_H_MPI_RECV_TAG((dir&1)^1), MPI_COMM_WORLD, &reqRecvsH[QPHIX_OPP_DIR(dir)]) == MPI_SUCCESS);
		//printf("rankd %d : done start receiving momentum from node %d... d = %d\n", myRank, neigh_ranks[QPHIX_OPP_DIR(dir)], QPHIX_OPP_DIR(dir));
		t1 = __rdtsc();
                t_boundary[3][QPHIX_OPP_DIR(dir)] += (t1 - t0);
        }
        t0 = __rdtsc();
        gf_pack_face_h_dir(tid, NCores, n_threads_per_core, hi, hcommsBuf2[8+dir], cb, dir>>1, dir&1);
        t1 = __rdtsc();
        tt_pack[tid][dir] += (t1 - t0);
        gBar->wait(tid);
//#pragma omp barrier
        if(tid == 0) {
                t1 = __rdtsc();
                t_boundary[0][dir] += (t1 - t0);
#ifndef ASSUME_MULTINODE
                MYASSERT(MPI_Isend((void *)hcommsBuf2[8+dir], bufsize*sizeof(fptype), MPI_BYTE, neigh_ranks[dir], GF_H_MPI_SEND_TAG(dir&1), MPI_COMM_WORLD, &reqSendsH[dir]) == MPI_SUCCESS);
#endif
                t0 = __rdtsc();
                t_boundary[2][dir] += (t0 - t1);
        }
}

void pack_momentum_face_dir(int tid, fptype *hi, int dir, int cb)
{
        unsigned long long t0, t1;
        int nv[4] = {Nx, Ny, Nz, Nt};
        long bufsize = 4*Nx*Ny*Nz*Nt/nv[dir>>1];
        if(tid == 0) {
                t0 = __rdtsc();
                MYASSERT(MPI_Irecv((void *)hcommsBuf2[QPHIX_OPP_DIR(dir)], bufsize*sizeof(fptype), MPI_BYTE, neigh_ranks[QPHIX_OPP_DIR(dir)], GF_H_MPI_RECV_TAG((dir&1)^1), MPI_COMM_WORLD, &reqRecvsH[QPHIX_OPP_DIR(dir)]) == MPI_SUCCESS);
                t1 = __rdtsc();
                t_boundary[3][QPHIX_OPP_DIR(dir)] += (t1 - t0);
        }
        t0 = __rdtsc();
        gf_pack_face_h_dir(tid, NCores, n_threads_per_core, hi, hcommsBuf2[8+dir], cb, dir>>1, dir&1);
        t1 = __rdtsc();
        tt_pack[tid][dir] += (t1 - t0);
        gBar->wait(tid);
//#pragma omp barrier
        if(tid == 0) {
                t1 = __rdtsc();
                t_boundary[0][dir] += (t1 - t0);
#ifndef ASSUME_MULTINODE
                MYASSERT(MPI_Isend((void *)hcommsBuf2[8+dir], bufsize*sizeof(fptype), MPI_BYTE, neigh_ranks[dir], GF_H_MPI_SEND_TAG(dir&1), MPI_COMM_WORLD, &reqSendsH[dir]) == MPI_SUCCESS);
#endif
                t0 = __rdtsc();
                t_boundary[2][dir] += (t0 - t1);
        }
}

void unpack_gauge_face_dir(int tid, Gauge *go, int dir, int cb)
{
    unsigned long long t0, t1, t2;
    if(tid == 0) {
	t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
        MYASSERT(MPI_Wait(&reqRecvsU[dir], MPI_STATUS_IGNORE) == MPI_SUCCESS);
#endif
        t1 = __rdtsc();
        t_boundary[5][dir] += (t1 - t0);
    }
    gBar->wait(tid);
//#pragma omp barrier
    t0 = __rdtsc();
#ifndef NO_COMPUTE
    gf_unpack_face_u_dir(tid, NCores, n_threads_per_core, gcommsBuf2[dir], go, cb, dir>>1, dir&1);
#endif
    t2 = __rdtsc();
    tt_unpack[tid][dir] += (t2 - t0);
}

void unpack_gauge_face_dir(int tid, fptype *go, int dir, int cb)
{
    unsigned long long t0, t1, t2;
    if(tid == 0) {
        t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
        MYASSERT(MPI_Wait(&reqRecvsU[dir], MPI_STATUS_IGNORE) == MPI_SUCCESS);
#endif
        t1 = __rdtsc();
        t_boundary[5][dir] += (t1 - t0);
    }
    gBar->wait(tid);
//#pragma omp barrier
    t0 = __rdtsc();
#ifndef NO_COMPUTE
    gf_unpack_face_u_dir(tid, NCores, n_threads_per_core, gcommsBuf2[dir], go, cb, dir>>1, dir&1);
#endif
    t2 = __rdtsc();
    tt_unpack[tid][dir] += (t2 - t0);
}

void unpack_momentum_face_dir(int tid, Hermit *ho, int dir, int cb)
{
    unsigned long long t0, t1, t2;
    if(tid == 0) {
        t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
        MYASSERT(MPI_Wait(&reqRecvsH[dir], MPI_STATUS_IGNORE) == MPI_SUCCESS);
#endif
        t1 = __rdtsc();
        t_boundary[5][dir] += (t1 - t0);
    }
    gBar->wait(tid);
//#pragma omp barrier
    t0 = __rdtsc();
#ifndef NO_COMPUTE
    gf_unpack_face_h_dir(tid, NCores, n_threads_per_core, hcommsBuf2[dir], ho, cb, dir>>1, dir&1);
#endif
    t2 = __rdtsc();
    tt_unpack[tid][dir] += (t2 - t0);
}

void unpack_momentum_face_dir(int tid, fptype *ho, int dir, int cb)
{
    unsigned long long t0, t1, t2;
    if(tid == 0) {
        t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
        MYASSERT(MPI_Wait(&reqRecvsH[dir], MPI_STATUS_IGNORE) == MPI_SUCCESS);
#endif
        t1 = __rdtsc();
        t_boundary[5][dir] += (t1 - t0);
    }
    gBar->wait(tid);
//#pragma omp barrier
    t0 = __rdtsc();
#ifndef NO_COMPUTE
    gf_unpack_face_h_dir(tid, NCores, n_threads_per_core, hcommsBuf2[dir], ho, cb, dir>>1, dir&1);
#endif
    t2 = __rdtsc();
    tt_unpack[tid][dir] += (t2 - t0);
}

/* Send gauge fields off-node */
/* Always send compressed gauges */
/* cb==1 : EVEN */
void gf_pack_and_send_boundaries_u(int tid, Gauge *gi, int cb)
{
  //for(int d = 3; d >= 0; d--) {
  for(int d=0; d<4; ++d) {
    if(!local_dir[d]) {
	unsigned long long t0, t1;
	if(tid == 0) {
	  //for(int fb = 1; fb >= 0; fb--) {
	  for(int fb=0; fb<2; ++fb) {
		t0 = __rdtsc();
		MYASSERT(MPI_Irecv((void *)gcommsBuf[4*d+fb], bufSizeG[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], GF_U_MPI_RECV_TAG(fb), MPI_COMM_WORLD, &reqRecvs[2*d+fb]) == MPI_SUCCESS);
		t1 = __rdtsc();
		t_boundary[3][2*d+fb] += (t1 - t0);
	  }
	}
	//for(int fb = 1; fb >=0; fb--) {
	for(int fb=0; fb<2; ++fb) {
		t0 = __rdtsc();
#ifndef NO_COMPUTE
		gf_pack_face_u(tid, NCores, n_threads_per_core, gi, gcommsBuf[2+4*d+fb], cb, d, fb);
#endif
		t1 = __rdtsc();
		tt_pack[tid][2*d+fb] += (t1 - t0);
		gBar->wait(tid);
//#pragma omp barrier
		if(tid == 0) {
			t1 = __rdtsc();
			t_boundary[0][2*d+fb] += (t1 - t0);
#ifndef ASSUME_MULTINODE
			MYASSERT(MPI_Isend((void *)gcommsBuf[2+4*d+fb], bufSizeG[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], GF_U_MPI_SEND_TAG(fb), MPI_COMM_WORLD, &reqSends[2*d+fb]) == MPI_SUCCESS);
#endif
			t0 = __rdtsc();
			t_boundary[2][2*d+fb] += (t0 - t1);
		}
	}
    }
  } // d loop
}

/* Receive dH off-node and update momentum in face */
/* cb == 1 : EVEN */
void gf_recv_and_unpack_boundaries_h(int tid, Hermit *hio, HermitHelperYZT *htmp, Gauge *gio, int cb)
{
   //bool accx=true;
   char accx = 0x00;
   for(int d = 0; d < 4; d++) if(!local_dir[d]) {
	unsigned long long t0, t1;
	for(int fb = 0; fb < 2; fb++) {
	    if(tid==0) {
		t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
		MYASSERT(MPI_Wait(&reqRecvs[8+2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
#endif
		t1 =  __rdtsc();
		t_boundary[5][2*d+fb] += (t1 - t0);
	    }
	    gBar->wait(tid);
//#pragma omp barrier
	    t0 = __rdtsc();
#ifndef NO_COMPUTE
	    gf_unpack_face_h(tid, NCores, n_threads_per_core, hcommsBuf[4*d+fb], hio, htmp, cb, d, fb, accx);
#endif
	    t1 =  __rdtsc();
	    tt_unpack[tid][2*d+fb] += (t1 - t0);
	    accx ^= (1<<(2*d+fb));
	}
   }
}

/* Receive gauge fields and calculate dH on boundary and send dH off-nodes */
/* cb == 1 : EVEN */
void gf_recv_and_unpack_u_and_send_boundaries_h(int tid, Hermit *hio, HermitHelperYZT *htmp, Gauge *gio, fptype kappaS, fptype kappaR, fptype kappaB, fptype epsilonH, int cb)
{
    fptype *gtmpbuf[8] = {0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
    char accumulate = 0xff;
    for(int d=0; d<4; ++d) if(!local_dir[d]) accumulate ^= (3 << (2*d));
    //unsigned long long t1[4];
#ifndef USE_WAITANY
    for(int d = 0; d < 4; ++d) {
	if(!local_dir[d]) {
	    unsigned long long t0, t1, t2;
	    for(int fb = 0; fb < 2; ++fb) {
		if(tid == 0) {
		    t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
		    MYASSERT(MPI_Wait(&reqRecvs[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
		    MYASSERT(MPI_Irecv((void *)hcommsBuf[4*d+fb], bufSizeH[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], GF_H_MPI_RECV_TAG(fb), MPI_COMM_WORLD, &reqRecvs[8+2*d+fb]) == MPI_SUCCESS);
#endif
		    t1 = __rdtsc();
		    t_boundary[5][2*d+fb] += (t1 - t0);
		}
		gBar->wait(tid);
//#pragma omp barrier
		accumulate ^= (1 << (2*d+fb));
		t0 = __rdtsc();
#ifndef NO_COMPUTE
/*
		if(tid==0) 
		{
		    printf("calling gf_unpack_face_u...\n");
		    fflush(stdout);
		}
*/
		gf_unpack_face_u(tid, NCores, n_threads_per_core, gcommsBuf[4*d+fb], gtmpbuf, gio, hio, htmp, kappaS, kappaR, kappaB, epsilonH, accumulate, cb, d, fb);
#endif
		gtmpbuf[2*d+fb] = gcommsBuf[4*d+fb];

		t2 = __rdtsc();
		tt_unpack[tid][2*d+fb] += (t2 - t0);
	    }
	}
    }
    gBar->wait(tid);
//#pragma omp barrier
    /* Send dH */
    if(tid == 0) for(int d = 0; d < 4; d++) 
	if(!local_dir[d]) {
	    unsigned long long t0, t2;
            for(int fb = 0; fb < 2; fb++) {
                    t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
		    MYASSERT(MPI_Isend((void *)hcommsBuf[2+4*d+fb], bufSizeH[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], GF_H_MPI_SEND_TAG(fb), MPI_COMM_WORLD, &reqSends[8+2*d+fb]) == MPI_SUCCESS);
#endif
		    t2 = __rdtsc();
		    t_boundary[1][2*d+fb] += (t2 - t0);	
	    }
	}
#else /* USE_WAITANY */
#error " GF USE_WAITANY NOT supported! "
#endif
    
}

#else // USE_CML_TG
void pack_gauge_face_dir(int tid, Gauge *gi, int dir, int cb)
{

}

void pack_gauge_face_dir(int tid, fptype *gi, int dir, int cb)
{

}

void pack_momentum_face_dir(int tid, Hermit *hi, int dir, int cb)
{

}

void pack_momentum_face_dir(int tid, fptype *hi, int dir, int cb)
{

}

void unpack_gauge_face_dir(int tid, Gauge *go, int dir, int cb)
{

}

void unpack_gauge_face_dir(int tid, fptype *go, int dir, int cb)
{

}

void unpack_momentum_face_dir(int tid, Hermit *ho, int dir, int cb)
{

}

void unpack_momentum_face_dir(int tid, fptype *ho, int dir, int cb)
{

}

void gf_pack_and_send_boundaries_u(int tid, Gauge *gi, int cb)
{
/* Not supported */
}

void gf_recv_and_unpack_boundaries_h(int tid, Hermit *hio, HermitHelperYZT *htmp, Gauge *gio, int cb)
{
/* Not supported */
}

void gf_recv_and_unpack_u_and_send_boundaries_h(int tid, Hermit *hio, HermitHelperYZT *htmp, Gauge *gio, int cb)
{
/* Not supported */
}

#endif

void gf_destroy_comms()
{
#if QPHIX_PrecisionInt == 2
    if(gcommsBuf[0]!=0x00)
	FREE(gcommsBuf[0]);
    if(gcommsBuf2[0]!=0x00)
	FREE(gcommsBuf2[0]);
    if(hcommsBuf[0]!=0x00)
	FREE(hcommsBuf[0]);
    if(hcommsBuf2[0]!=0x00)
	FREE(hcommsBuf2[0]);
    for(int i=0; i<16; ++i)
    {
	gcommsBuf[i]=0x00;
	gcommsBuf2[i]=0x00;
	hcommsBuf[i]=0x00;
	hcommsBuf2[i]=0x00;
    }
    for(int i=0; i<8; ++i)
    {
	glocalBuf[i]=0x00;
	hlocalBuf[i]=0x00;
    }
#endif
    for(int i=0; i<2; ++i) for(int j=0; j<4; ++j)
        PackMask[i][j] = 0;
}

#if QPHIX_PrecisonInt == 1

void gf_print_boundary_timings(unsigned long long t_gf){

}

#endif

#endif /* ENABLE_MPI */

