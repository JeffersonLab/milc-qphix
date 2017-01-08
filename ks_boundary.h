#ifndef _BOUNDARY_H_
#define _BOUNDARY_H_
#include <unistd.h>

extern int Gx, Gy, Gz, Gt;
extern int Nxh, Nx, Ny, Nz, Nt;
//extern int Vxh, Vx, Vy, Vz, Vt;
//extern int Lsxh, Lsx, Lsy, Lsz, Lst;
extern int n_threads_per_core;
//extern int Pxy, Pxyz;
extern int NCores;
extern Barrier* gBar;
extern int geometry[4];
extern int myRank;
extern int nRanks;
extern bool local_dir[4];
extern bool printAllRanks;
extern int nParallelRecvs, nParallelSends;
extern double Freq;

double getFreq();

void setup_mpi(int argc, char *argv[]);
void setup_boundaries(void);
void setup_comms_F(int myCoord_[4], int neigh_ranks_[8]);
void setup_comms_D(int myCoord_[4], int neigh_ranks_[8]);
void destroy_comms_F();
void destroy_comms_D();
void pack_and_send_boundaries(int tid, KS *si, int cb);
void recv_and_unpack_boundaries(int tid, KS *so, Gauge18 *u, Gauge *ull, const fptype beta, int cb);
void print_boundary_timings(int iters, unsigned long long t_dslash);
double mpi_allreduce(double sum);

#define MSEC(x) ((double)x/Freq*1000)
#define MYASSERT(cond) if(!(cond)) { printf("Rank %d: %s:%d Assertion failed\n\t%s evaluates false\n", myRank, __FILE__, __LINE__, #cond); exit(1);}
#ifndef ENABLE_MPI

/*
inline void pack_and_send_boundaries(int tid, KS *si, int cb) {}
inline void recv_and_unpack_boundaries(int tid, KS *so, Gauge18 *u, Gauge *ull, const fptype beta, int cb) {}
inline void print_boundary_timings(int iters, unsigned long long t_dslash) {}
inline double mpi_allreduce(double sum) { return sum; }
*/

#else // ENABLE_MPI is defined
#ifndef ASSUME_MULTINODE
#include <mpi.h>
#endif

#define DSLASH_MPI_SEND_TAG(x) (12+(x))
#define DSLASH_MPI_RECV_TAG(x) (12+(1-x))
#if 0
int face_npkts[4] = {0}, face_pktsize[4] = {0};
unsigned int face_mask_x[2][2] = {0};
unsigned int face_mask_y[2] = {0};
__declspec(align(64)) int face_offs[VECLEN] = {0};
__declspec(align(64)) int face_gOffs[VECLEN] = {0};
typedef struct {
	short dims[4];
} face_pkt_info;
face_pkt_info *pkt_info_buf[4];
#endif

void setup_cml_thread_groups(void);

void pack_face(int tid, int num_cores, int threads_per_core, KS *si, fptype *res, int cb, int dir, int fb);
void unpack_face(int tid, int num_cores, int threads_per_core, fptype *inbuf, KS *res, Gauge18 *u, Gauge * ull, const fptype beta, int cb, int dir, int fb);

#if 0
inline void pack_face(int tid, int num_cores, int threads_per_core, KS *si, fptype *res, int cb, int dir, int fb) {
	int nyg_pack = 1;
	if(dir == 0) nyg_pack = 2;
	int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	int shifts[4] = {1, Vxh, Pxy, Pxyz};

	lens[dir] = 1;

	int npkts = lens[0] * lens[1] * lens[2] * lens[3];

	int pktsize = VECLEN;
	int is1stAnd3rdFused = 0;
	if((1 << dir) < VECLEN) {
		pktsize = VECLEN/2;
		is1stAnd3rdFused = (dir == 0 ? (Vxh == 1) : (lens1[dir] == 2));
	}

	int shift = shifts[dir] * (fb == 0 ? 1 : -1);

	// printf("pkts = %d, pktsize=%d\n", npkts, pktsize);
	MYASSERT(18*pktsize*npkts == bufSize[dir]);

	int pkts_per_core = (npkts + num_cores - 1) / num_cores;

	// Assume thread indexing runs SMT-ID fastest
	// smtid + threads_per_core*cid
	int cid = tid/threads_per_core;
	int smtid = tid - threads_per_core * cid;

	int low_pkt = cid*pkts_per_core;
	int high_pkt = (cid+1)*pkts_per_core;
	if ( high_pkt > npkts ) high_pkt = npkts;

	// OK Each core can now work on its vectors:
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
		int xshift=0;

		int xodd = (z + t + cb) & 1;

		if(dir == 0) {
			int row = (fb == 1 ? 1-xodd : xodd); // picking y with single pack site first
			y += row; 
			xshift = (row == 0 ? Vxh : -Vxh);
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
		// yi is always going to be even for dir==0
		int xodd_next = (z_next + t_next + cb) & 1;

		if(dir == 0) y_next +=(fb == 1 ? 1-xodd_next : xodd_next);
		// printf("Pack %d %d %d %d\n", t, z, yi, xblock*SOALEN);

		// Now we have x,y,z coordinates, we need the base address of the KS
		KS *siBase = &si[t*Pxyz + z*Pxy + y*Vxh + x];
		// Offset to next KS. T-values are the same.
		int off_next = (t_next-t)*Pxyz+(z_next-z)*Pxy+(y_next-y)*Vxh+(x_next-x);

		// Convert to prefetch distance in floats
		//int si_offset = off_next*sizeof(KS)/sizeof(fptype);
		//int hsprefdist = (pkt_next - pkt)*sizeof(KS)/sizeof(fptype)/2;

		// We are streaming out in sequence
		fptype * __restrict rBuf = &res[18*pktsize*pkt];
		fptype * __restrict lBuf = &localBuf[2*dir+fb][18*pktsize*pkt];

		//printf("rank = %d, pkt = %d, outbuf=%p (%lld)\n", myRank, pkt, outbuf, outbuf-res);
		// OK: now we have xyBase, offs, and oubuf -- we should call the kernel.
		if(dir == 0) {
			ks_face_pack_dir(siBase, lBuf, rBuf, dir*2+fb);
			ks_face_pack_dir(siBase+xshift, lBuf+6*pktsize, rBuf+6*pktsize, dir*2+fb);
			if(is1stAnd3rdFused) {
				rBuf[12*pktsize:6*pktsize] = lBuf[6*pktsize:6*pktsize];
			}
			else {
				ks_face_pack_dir(siBase+xshift+shift, lBuf+12*pktsize, rBuf+12*pktsize, dir*2+fb);
			}
		}
		else {
			ks_face_pack_dir(siBase, lBuf, rBuf, dir*2+fb);
			ks_face_pack_dir(siBase+shift, lBuf+6*pktsize, rBuf+6*pktsize, dir*2+fb);
			if(is1stAnd3rdFused) {
				rBuf[12*pktsize:6*pktsize] = lBuf[0:6*pktsize];
			}
			else {
				ks_face_pack_dir(siBase+2*shift, lBuf+12*pktsize, rBuf+12*pktsize, dir*2+fb);
			}
		}
	}
}

inline void unpack_face(int tid, int num_cores, int threads_per_core, fptype *inbuf, KS *res, Gauge18 *u, Gauge * ull, const fptype beta, int cb, int dir, int fb) {
	int nyg_pack = 1;
	if(dir == 0) nyg_pack = 2;
	int lens1[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	int lens[4] = {Vxh, Vy/nyg_pack, Vz, Vt};
	int shifts[4] = {1, Vxh, Pxy, Pxyz};

	lens[dir] = 1;

	int npkts = lens[0] * lens[1] * lens[2] * lens[3];

	int pktsize = VECLEN;
	int is1stAnd3rdFused = 0;
	if((1 << dir) < VECLEN) {
		pktsize = VECLEN/2;
		is1stAnd3rdFused = (dir == 0 ? (Vxh == 1) : (lens1[dir] == 2));
	}

	int shift = shifts[dir] * (fb == 0 ? 1 : -1);

	// printf("pkts = %d, pktsize=%d\n", npkts, pktsize);
	MYASSERT(18*pktsize*npkts == bufSize[dir]);

	int pkts_per_core = (npkts + num_cores - 1) / num_cores;

	// Assume thread indexing runs SMT-ID fastest
	// smtid + threads_per_core*cid
	int cid = tid/threads_per_core;
	int smtid = tid - threads_per_core * cid;

	int low_pkt = cid*pkts_per_core;
	int high_pkt = (cid+1)*pkts_per_core;
	if ( high_pkt > npkts ) high_pkt = npkts;

	// OK Each core can now work on its vectors:
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
		int xshift=0;

		int xodd = (z + t + cb) & 1;

		if(dir == 0) {
			int row = (fb == 0 ? 1-xodd : xodd); // picking y with single pack site first
			y += row; 
			xshift = (row == 0 ? Vxh : -Vxh); 
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
		// yi is always going to be even for dir==0
		int xodd_next = (z_next + t_next + cb) & 1;

		if(dir == 0) y_next +=(fb == 1 ? xodd_next : 1-xodd_next);
		// printf("Unpack %d %d %d %d\n", t, z, yi, xblock*SOALEN);

		// Now we have x,y,z coordinates, we need the base address of the KS
		KS *oBase = &res[t*Pxyz + z*Pxy + y*Vxh + x];
		// Offset to next KS. T-values are the same.
		int off_next = (t_next-t)*Pxyz+(z_next-z)*Pxy+(y_next-y)*Vxh+(x_next-x);

		//int soprefdist = off_next*sizeof(KS)/sizeof(fptype);
		//int hsprefdist = (pkt_next - pkt)*sizeof(KS)/sizeof(fptype)/2;
		//int goff_next  = ((t_next-t)*Pxyz+(z_next-z)*Pxy+(yi_next-yi)*nvecs)/nyg + (xblock_next-xblock);
		//int gprefdist = goff_next*sizeof(Gauge)/sizeof(fptype);

		Gauge18 *gBase = &u[t*Pxyz+z*Pxy+y*Vxh+x];
		Gauge *gllBase = &ull[t*Pxyz+z*Pxy+y*Vxh+x];
		// We are streaming out in sequence
		fptype * __restrict rBuf = &inbuf[18*pktsize*pkt];
		fptype * __restrict lBuf = &localBuf[2*dir+1-fb][18*pktsize*pkt];

		// OK: now we have xyBase, offs, and oubuf -- we should call the kernel.

		if(dir == 0) {
			// single 3rd neighbor
			ks_dslash_face_unpack_dir(lBuf, rBuf, gllBase, oBase, beta, xodd, y, t, dir*2+fb);
			// process first neighbor
			ks_dslash_face_unpack18_dir(lBuf+6*pktsize, rBuf+6*pktsize, gBase+xshift, oBase+xshift, beta, xodd, y, t, dir*2+fb);
			if(is1stAnd3rdFused) {
				ks_dslash_face_unpack_dir(rBuf+6*pktsize, rBuf+12*pktsize, gllBase+xshift, oBase+xshift, beta, xodd, y, t, dir*2+fb);
			}
			else {
				// Now inner 3rd neighbour
				ks_dslash_face_unpack_dir(lBuf+6*pktsize, rBuf+6*pktsize, gllBase+xshift+shift, oBase+xshift+shift, beta, xodd, y, t, dir*2+fb);
				// outer 3rd neighbor
				ks_dslash_face_unpack_dir(lBuf+12*pktsize, rBuf+12*pktsize, gllBase+xshift, oBase+xshift, beta, xodd, y, t, dir*2+fb);
			}
		}
		else {
			// process first neighbor first
			ks_dslash_face_unpack18_dir(lBuf, rBuf, gBase, oBase, beta, xodd, y, t, dir*2+fb);
			//need to check for fused 1st and 3rd 3rd-neighbour
			if(is1stAnd3rdFused) {
				ks_dslash_face_unpack_dir(rBuf, rBuf+12*pktsize, gllBase, oBase, beta, xodd, y, t, dir*2+fb);
			} 
			else {
				// Now inner most 3rd neighbour
				ks_dslash_face_unpack_dir(lBuf, rBuf, gllBase+2*shift, oBase+2*shift, beta, xodd, y, t, dir*2+fb);
				// outer most 3rd neighbor
				ks_dslash_face_unpack_dir(lBuf+12*pktsize, rBuf+12*pktsize, gllBase, oBase, beta, xodd, y, t, dir*2+fb);
			}
			// middle 3rd neighbor
			ks_dslash_face_unpack_dir(lBuf+6*pktsize, rBuf+6*pktsize, gllBase+shift, oBase+shift, beta, xodd, y, t, dir*2+fb);
		}
	}
}


#ifndef USE_CML_TG
inline void pack_and_send_boundaries(int tid, KS *si, int cb) {
	for(int d = 3; d >= 0; d--) {
		if(!local_dir[d]) {
			unsigned long long t0, t1;
			if(tid == 0) {
				for(int fb = 0; fb < 2; fb++) {
					//printf("MPI_Irecv: src=%d, dst=%d, dir=%d, fb=%d, size=%lld, buf=%p\n", myRank, neigh_ranks[2*d+fb], d, fb, bufSize[d]*sizeof(fptype), commsBuf[4*d+fb]);
					t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
					MYASSERT(MPI_Irecv((void *)commsBuf[4*d+fb], bufSize[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], DSLASH_MPI_RECV_TAG(fb), MPI_COMM_WORLD, &reqRecvs[2*d+fb]) == MPI_SUCCESS);
#endif
					t1 = __rdtsc();
					t_boundary[3][2*d+fb] += (t1 - t0);
				}
			}
			for(int fb = 1; fb >=0; fb--) {
				//printf("pack buf size = %d s = %p, e = %p\n", bufSize[d], commsBuf[2+4*d+fb], commsBuf[2+4*d+fb] + bufSize[d]);
				t0 = __rdtsc();
#ifndef NO_COMPUTE
				pack_face(tid, NCores, n_threads_per_core, si, commsBuf[2+4*d+fb], cb, d, fb);
#endif
				t1 = __rdtsc();
				tt_pack[tid][2*d+fb] += (t1 - t0);
				gBar->wait(tid);
				if(tid == 0) {
					//printf("MPI_Isend: src=%d, dst=%d, dir=%d, fb=%d, size=%lld, buf=%p\n", myRank, neigh_ranks[2*d+fb], d, fb, bufSize[d]*sizeof(fptype), commsBuf[2+4*d+fb]);
					t1 = __rdtsc();
					t_boundary[0][2*d+fb] += (t1 - t0);
#ifndef ASSUME_MULTINODE
					MYASSERT(MPI_Isend((void *)commsBuf[2+4*d+fb], bufSize[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], DSLASH_MPI_SEND_TAG(fb), MPI_COMM_WORLD, &reqSends[2*d+fb]) == MPI_SUCCESS);
#endif
					t0 = __rdtsc();
					t_boundary[2][2*d+fb] += (t0 - t1);
				}
				//gBar->wait(tid);
			}
		}
	}
}

volatile int waitany_dir;
inline void recv_and_unpack_boundaries(int tid, KS *so, Gauge18 *u, Gauge *ull, const fptype beta, int cb) {
	//gBar->wait(tid);
#ifndef USE_WAITANY
	for(int d = 3; d >= 0; d--) {
		if(!local_dir[d]) {
			unsigned long long t0, t1, t2;
			for(int fb = 0; fb < 2; fb++) {
				if(tid == 0) {
					t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
					MYASSERT(MPI_Wait(&reqRecvs[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
#endif
					t1 = __rdtsc();
					t_boundary[5][2*d+fb] += (t1 - t0);
				}
				gBar->wait(tid);
				t0 = __rdtsc();
#ifndef NO_COMPUTE
				unpack_face(tid, NCores, n_threads_per_core, commsBuf[4*d+fb], so, u, ull, beta, cb, d, fb);
#endif
				t2 = __rdtsc();
				tt_unpack[tid][2*d+fb] += (t2 - t0);
				if(tid == 0) t_boundary[1][2*d+fb] += (t2 - t1);
			}
		}
	}
	for(int d = 3; d >= 0; d--) {
		if(!local_dir[d]) {
			unsigned long long t0, t1;
			if(tid == 0) {
				for(int fb = 1; fb >= 0; fb--) {
					t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
					MYASSERT(MPI_Wait(&reqSends[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
#endif
					t1 = __rdtsc();
					t_boundary[4][2*d+fb] += (t1 - t0);
				}
			}
		}
	}
#else // USE_WAITANY
	for(int i = 0; i < nCommsDirs; i++) {
		int d, fb;
		unsigned long long t0, t1, t2;
		if(tid == 0) {
			int dir;
			t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
			MYASSERT(MPI_Waitany(8, reqRecvs, &dir, MPI_STATUS_IGNORE) == MPI_SUCCESS);
#else
			dir = i;
#endif
			t1 = __rdtsc();
			t_boundary[5][dir] += (t1 - t0);
			waitany_dir = dir;
		}
		gBar->wait(tid);
		d = waitany_dir / 2;
		fb = waitany_dir & 1;
		t0 = __rdtsc();
#ifndef NO_COMPUTE
		unpack_face(tid, NCores, n_threads_per_core, commsBuf[4*d+fb], so, u, ull, beta, cb, d, fb);
#endif
		t2 = __rdtsc();
		tt_unpack[tid][2*d+fb] += (t2 - t0);
		if(tid == 0) t_boundary[1][waitany_dir] += (t2 - t1);
	}
	for(int d = 3; d >= 0; d--) {
		if(!local_dir[d]) {
			unsigned long long t0, t1;
			if(tid == 0) {
				for(int fb = 1; fb >= 0; fb--) {
					int dir;
					t0 = __rdtsc();
#ifndef ASSUME_MULTINODE
					MYASSERT(MPI_Waitany(8, reqSends, &dir, MPI_STATUS_IGNORE) == MPI_SUCCESS);
#else
					dir = 2*d+fb;
#endif
					t1 = __rdtsc();
					t_boundary[4][dir] += (t1 - t0);
				}
			}
		}
	}
#endif // USE_WAITANY
}

#else // USE_CML_TG
inline void pack_and_send_boundaries(int tid, KS *si, int cb) {
	//printf("Rank %d: Tid %d: Inside %s @%d\n", myRank, tid, __func__, __LINE__);
	int sendTgNum = 0, recvTgNum = 0;
	for(int d = 3; d >= 0; d--) {
		if(!local_dir[d]) {
			unsigned long long t0, t1;
			for(int fb = 0; fb < 2; fb++) {
				int tgtid = -1;
				if(CML_Threadgroup_member(tgRecvs[recvTgNum], tid, &tgtid) == MPI_SUCCESS) {
					t0 = __rdtsc();
					MYASSERT(CML_Thread_irecv(tgRecvs[recvTgNum], tid, (void *)commsBuf[4*d+fb], bufSize[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], DSLASH_MPI_RECV_TAG(fb), MPI_COMM_WORLD, &cmlReqRecvs[2*d+fb]) == MPI_SUCCESS);
					t1 = __rdtsc();
					if(tgtid==0) t_boundary[3][2*d+fb] += (t1 - t0);
				}
				recvTgNum++;
				if(recvTgNum == nParallelRecvs) recvTgNum = 0;
			}
			for(int fb = 1; fb >=0; fb--) {
				int tgtid = -1;
				if(CML_Threadgroup_member(tgSends[sendTgNum], tid, &tgtid) == MPI_SUCCESS) {
					t0 = __rdtsc();
#ifndef NO_COMPUTE
					pack_face(tgtid, numSendCores, numSendThPerCore, si, commsBuf[2+4*d+fb], cb, d, fb);
#endif
					t1 = __rdtsc();
					tt_pack[tgtid][2*d+fb] += (t1 - t0);
					//gBar->wait(tid);
					if(tgtid==0) t_boundary[0][2*d+fb] += (t1 - t0);
					MYASSERT(CML_Thread_isend(tgSends[sendTgNum], tid, CML_NO_AGGREGATE, (void *)commsBuf[2+4*d+fb], (void *)commsBuf[2+4*d+fb], bufSize[d]*sizeof(fptype), MPI_BYTE, neigh_ranks[2*d+fb], DSLASH_MPI_SEND_TAG(fb), MPI_COMM_WORLD, &cmlReqSends[2*d+fb]) == MPI_SUCCESS);
					t0 = __rdtsc();
					if(tgtid==0) t_boundary[2][2*d+fb] += (t0 - t1);
					//gBar->wait(tid);
				}
				sendTgNum++;
				if(sendTgNum == nParallelSends) sendTgNum = 0;
			}
		}
	}
}

volatile int waitany_dir;
inline void recv_and_unpack_boundaries(int tid, KS *so, Gauge18 *u, Gauge *ull, const fptype beta, int cb) {
//	gBar->wait(tid);
	int sendTgNum = 0, recvTgNum = 0;
#ifndef USE_WAITANY
	for(int d = 3; d >= 0; d--) {
		if(!local_dir[d]) {
			unsigned long long t0, t1, t2;
			for(int fb = 0; fb < 2; fb++) {
				int tgtid = -1;
				if(CML_Threadgroup_member(tgRecvs[recvTgNum], tid, &tgtid) == MPI_SUCCESS) {
					t0 = __rdtsc();
					MYASSERT(CML_Thread_wait(tgRecvs[recvTgNum], tid, &cmlReqRecvs[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
					t1 = __rdtsc();
					if(tgtid==0) t_boundary[5][2*d+fb] += (t1 - t0);
					gBar->wait(tid);
#ifndef NO_COMPUTE
					unpack_face(tgtid, numRecvCores, numRecvThPerCore, commsBuf[4*d+fb], so, u, beta, cb, d, fb);
#endif
					t2 = __rdtsc();
					tt_unpack[tgtid][2*d+fb] += (t2 - t1);
					if(tgtid==0) t_boundary[1][2*d+fb] += (t2 - t1);
				}
				recvTgNum++;
				if(recvTgNum == nParallelRecvs) recvTgNum = 0;
			}
		}
	}
	for(int d = 3; d >= 0; d--) {
		if(!local_dir[d]) {
			unsigned long long t0, t1;
			for(int fb = 1; fb >= 0; fb--) {
				int tgtid = -1;
				if(CML_Threadgroup_member(tgSends[sendTgNum], tid, &tgtid) == MPI_SUCCESS) {
					t0 = __rdtsc();
					MYASSERT(CML_Thread_wait(tgSends[sendTgNum], tid, &cmlReqSends[2*d+fb], MPI_STATUS_IGNORE) == MPI_SUCCESS);
					t1 = __rdtsc();
					if(tgtid==0) t_boundary[4][2*d+fb] += (t1 - t0);
				}
				sendTgNum++;
                                if(sendTgNum == nParallelSends) sendTgNum = 0;
			}
		}
	}
#else // USE_WAITANY
#error "Waitany version is not yet supported for CML_Threading"
#endif // USE_WAITANY
}
#endif

inline double mpi_allreduce(double sum) {
	double globalSum;
#ifndef ASSUME_MULTINODE
	MPI_Allreduce(&sum, &globalSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	globalSum = sum;
#endif
	return globalSum; 
}

inline void print_boundary_timings(int iters, unsigned long long t_dslash) {
	const char *dirs[] = {"X Back", "X Forw", "Y Back", "Y Forw", "Z Back", "Z Forw", "T Back", "T Forw" };
	for(int k = 0; k < nRanks; k++) {
		if(!printAllRanks && k > 0) break; 
#ifndef ASSUME_MULTINODE
		MPI_Barrier(MPI_COMM_WORLD);
#endif
		if(myRank == k) {
			unsigned long long total[8] = {0}, gTotal = 0, dtotal;
			char tmpBuf[4096];
			int pos = 0;
			pos+=snprintf(&tmpBuf[pos], 4095-pos, "Average DSlash time: %10.3f msec\n", MSEC(t_dslash/iters));
			pos+=snprintf(&tmpBuf[pos], 4095-pos, "Boundary Timings in msec:\n  Dir    - %10s %10s %10s %10s %10s %10s %10s\n", "Pack", "Unpack", "Send", "Recv", "SendWait", "RecvWait", "Total");
			for(int d = 0; d < 8; d++) {
				dtotal = 0;
				if(!local_dir[d/2]) {
					pos+=snprintf(&tmpBuf[pos], 4095-pos, "  %s - ", dirs[d]);
					for(int i = 0; i < 6; i++) {
						pos+=snprintf(&tmpBuf[pos], 4095-pos, "%10.3f ", MSEC(t_boundary[i][d]/iters));
						total[i] += t_boundary[i][d];
						dtotal += t_boundary[i][d];
						gTotal += t_boundary[i][d];
					}
					pos+=snprintf(&tmpBuf[pos], 4095-pos, "%10.3f %10.3f %10.3f\n", MSEC(dtotal/iters), MSEC(tt_pack[0][d]/iters), MSEC(tt_unpack[0][d]/iters));
					total[6] += tt_pack[0][d];
					total[7] += tt_unpack[0][d];
				}
			}
			if(gTotal > 0) {
				pos+=snprintf(&tmpBuf[pos], 4095-pos, "  %s - ", "Total ");
				for(int i = 0; i < 6; i++) {
					pos+=snprintf(&tmpBuf[pos], 4095-pos, "%10.3f ", MSEC(total[i]/iters));
				}
				pos+=snprintf(&tmpBuf[pos], 4095-pos, "%10.3f %10.3f %10.3f\n", MSEC(gTotal/iters), MSEC(total[6]/iters), MSEC(total[7]/iters));
			}
			tmpBuf[pos] = 0;
#ifdef USE_CML_TG
			for(int d = 0; d < numTgs; d++) {
				CML_Threadgroup_time_print(tgList[d], d==0);
			}
#endif
			printf("Rank %3d:\n%s\n", myRank, tmpBuf); 
		}
	}
}
#endif

#endif // ENABLE_MPI

#endif // _BOUNDARY_H_

