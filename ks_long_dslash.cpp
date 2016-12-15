#include <cassert>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <immintrin.h>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "ks_config.h"
#include "ks_globals.h"        /* All global values and macros */
#include "ks_long_dslash.h"
#include "ks_dslash_complete_specialization.h"
#include "ks_boundary.h"

#undef allocKS_fptype
#undef allocKS2_fptype
#undef allocGauge_fptype
#undef allocGauge18_fptype
#undef setKS_fptype
#undef getKS_fptype
#undef setFullGauge18_fptype
#undef setFullGauge_fptype
#undef setGauge_fptype
#undef getGauge_fptype
#undef setGauge18_fptype
#undef getGauge18_fptype
#undef qphix_ks_dslash_fptype

/* defining the globals */
#if QPHIX_PrecisionInt==1
#define allocKS_fptype allocKS_F
#define allocKS2_fptype allocKS2_F
#define allocGauge_fptype allocGauge_F
#define allocGauge18_fptype allocGauge18_F
#define setKS_fptype setKS_F
#define getKS_fptype getKS_F
#define setFullGauge18_fptype setFullGauge18_F
#define setFullGauge_fptype setFullGauge_F
#define setGauge_fptype setGauge_F
#define getGauge_fptype getGauge_F
#define setGauge18_fptype setGauge18_F
#define getGauge18_fptype getGauge18_F
#define qphix_ks_dslash_fptype qphix_ks_dslash_F

#elif QPHIX_PrecisionInt==2
#define allocKS_fptype allocKS_D
#define allocKS2_fptype allocKS2_D
#define allocGauge_fptype allocGauge_D
#define allocGauge18_fptype allocGauge18_D
#define setKS_fptype setKS_D
#define getKS_fptype getKS_D
#define setFullGauge18_fptype setFullGauge18_D
#define setFullGauge_fptype setFullGauge_D
#define setGauge_fptype setGauge_D
#define getGaug_fptypee getGauge_D
#define setGauge18_fptype setGauge18_D
#define getGauge18_fptype getGauge18_D
#define qphix_ks_dslash_fptype qphix_ks_dslash_D

#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif

#if 0
EXTERN int Gxh, Gx, Gy, Gz, Gt;
EXTERN int Nxh, Nx, Ny, Nz, Nt;
EXTERN int Vxh, Vx, Vy, Vz, Vt;
EXTERN int Lsxh, Lsx, Lsy, Lsz, Lst;

EXTERN int Pxy, Pxyz;
EXTERN int NCores;
EXTERN int num_floats_in_ks_array;
EXTERN int n_threads_per_core;
EXTERN int nThreads;
EXTERN int geometry[4];
EXTERN int myRank;
EXTERN int nRanks;
EXTERN bool printAllRanks;
EXTERN int nParallelRecvs, nParallelSends;

EXTERN Phaser *phaser;
EXTERN Phaser *phaser_dslash;
EXTERN Barrier* gBar;
EXTERN bool local_dir[4];
#endif

#if QPHIX_PrecisionInt==1
int Gxh, Gx, Gy, Gz, Gt;
int Nxh, Nx, Ny, Nz, Nt;
int Vxh, Vx, Vy, Vz, Vt;
int Lsxh, Lsx, Lsy, Lsz, Lst;

int Pxy, Pxyz;
int NCores;
int num_floats_in_ks_array;
int n_threads_per_core;
int nThreads;
int geometry[4];
int myRank;
int nRanks;
bool printAllRanks;
int nParallelRecvs, nParallelSends;

Phaser *phaser;
//Phaser *phaser_dslash;
Barrier* gBar;
bool local_dir[4];
#endif

const int nGX = (VECLEN < 2 ? 1 : 2);
const int nGY = (VECLEN < 4 ? 1 : 2);
const int nGZ = (VECLEN < 8 ? 1 : 2);
const int nGT = (VECLEN < 16 ? 1 : 2);
 
#if QPHIX_PrecisionInt==1
const int MILC2MBENCH[8] = { 1, 3, 5, 7, 6, 4, 2, 0 };
const int MBENCH2MILC[8] = { 7, 0, 6, 1, 5, 2, 4, 3 };
#endif

#if QPHIX_PrecisionInt==1
double Freq = getFreq();
char * BoundTable = 0x00;
unsigned int * NeighTable = 0x00;
//unsigned int * Neigh3Table = 0x00;
int BLENGTH = 0;
int PadBound = 0;
int PadNeigh = 0;
#elif QPHIX_PrecisionInt==2
extern double Freq;
extern char * BoundTable;
extern unsigned int * NeighTable;
//unsigned int * Neigh3Table = 0x00;
extern int BLENGTH;
extern int PadBound;
extern int PadNeigh;
#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif

static
void dump_ks_spinor(void *si_arg)
{
    KS* si = (KS*)si_arg;
    //dump_spinor(si);

}

// r = (s1*s2-s3*s4)'
//r[0] = (s1[0]*s2[0])-(s1[1]*s2[1])-(s3[0]*s4[0])+(s3[1]*s4[1])
//r[1] = (s3[0]*s4[1])+(s3[1]*s4[0])-(s1[0]*s2[1])-(s1[1]*s2[0])
void 
Conj_CrossProd(fptype *r, fptype *s1, fptype *s2, fptype *s3, fptype *s4)
{
    r[0] =  s1[0] *  s2[0];
    r[0] =  r[0] -  s1[1] *  s2[1];
    r[0] =  r[0] -  s3[0] *  s4[0];
    r[0] =  r[0] +  s3[1] *  s4[1];

    r[1] =  s3[0] *  s4[1];
    r[1] =  r[1] +  s3[1] *  s4[0];
    r[1] =  r[1] -  s1[0] *  s2[1];
    r[1] =  r[1] -  s1[1] *  s2[0];
}

#if 0
void 
setup_mbench( int *lattice, int *geom, int *my_coord, int *my_neigh_ranks, int my_id
            , int NCores_, int ThreadsPerCore, int minCt_)
{
    int sy =1, sz = 1;
    Gx = lattice[0];
    Gxh = Gx/2; 
    Gy = lattice[1];
    Gz = lattice[2]; 
    Gt = lattice[3]; 
    NCores = NCores_; 
    n_threads_per_core = ThreadsPerCore;
    myRank = my_id;
    nRanks = 1;
    for(int i = 0; i < 4; i++)
    {
        geometry[i] = geom[i];
        nRanks *= geom[i];
    }
#ifndef ENABLE_MPI
		if(nRanks != 1) {
			printf("MPI not enabled in QPhiX, please compile with ENABLE_MPI=1\n"); exit(1);
		}
#else
		//printf("Rank: %d (pid %d) waiting for gdb to connect...\n", myRank, getpid());
		//volatile int dummy = 0;
		//while(dummy == 0) ;
		MPI_Barrier(MPI_COMM_WORLD);
#endif
    setup_comms(my_coord, my_neigh_ranks);

		//checkParams();
		MYASSERT(Nxh % 2 == 0);
		MYASSERT(Ny % 4 == 0);
		MYASSERT(Nz % 4 == 0); 
		MYASSERT(Nt % 4 == 0);

		Vx = (VECLEN > 1 ? Nx/2 : Nx);
		Vxh = (VECLEN > 1 ? Nxh/2 : Nxh);
		Vy = (VECLEN > 2 ? Ny/2 : Ny);
		Vz = (VECLEN > 4 ? Nz/2 : Nz);
		Vt = (VECLEN > 8 ? Nt/2 : Nt);
		MYASSERT(Vy % BY == 0);
		MYASSERT(Vz % BZ == 0);
		Pxy = (Vxh*Vy+XY_PAD);
		Pxyz = (Pxy*Vz+XYZ_PAD);

    //FIXME: need better way to set sy and sz
    if(ThreadsPerCore == 1) sz = 1;
    else if(ThreadsPerCore == 2) sz = 2;
    else if(ThreadsPerCore == 4) {
        sy = 2;
        sz = 2;
    }
    else {
        printf("Threads per core must be set to 1, 2 or 4. Value is %d\n", ThreadsPerCore);
        exit(1);
    }

    num_floats_in_ks_array = (Pxyz*Vt)* sizeof(KS)/sizeof(fptype);

		gBar = new Barrier(NCores,ThreadsPerCore);
		phaser = new Phaser(Vxh, Vy, Vz, Vt, Vxh
                      , BY, BZ, NCores, sy, sz, minCt_);

    nThreads = NCores * ThreadsPerCore;

#ifdef _OPENMP
#pragma omp parallel default (shared) num_threads(nThreads)
    {
        int tid = omp_get_thread_num();
				gBar->init(tid);
        phaser->init(tid);
    }
#else
    if(nThreads != 1) {printf("No. of threads set to %d while not using OpenMP\n", nThreads); }
    phaser->init(0);  
#endif    
 printf("Done setting up mbench\n");
}
#endif

static void 
*initBuf(void *b, size_t s)
{
	fptype *buf = (fptype*)b;
#pragma omp parallel for num_threads (nThreads)
	for(int i = 0; i < s/sizeof(fptype); i++)
		buf[i] = 0.0f;
	return b;
}

void 
* allocKS_fptype()
{
    return initBuf(_mm_malloc((Pxyz*Vt+1)*sizeof(KS), 64), (Pxyz*Vt+1)*sizeof(KS));
}

/* For both parities */
void 
* allocKS2_fptype()
{
    return initBuf(_mm_malloc((Pxyz*Vt+1)*2*sizeof(KS), 64), (Pxyz*Vt+1)*2*sizeof(KS));
}

#if QPHIX_PrecisionInt==1
void 
freeKS(void *p)
{
    _mm_free(p);
}
#endif

void 
* allocGauge_fptype()
{
        return initBuf(_mm_malloc((Pxyz*Vt)*sizeof(Gauge), 64)
                      , (Pxyz*Vt)*sizeof(Gauge));
}

void 
* allocGauge18_fptype()
{
        return initBuf(_mm_malloc((Pxyz*Vt)*sizeof(Gauge18), 64)
                       , (Pxyz*Vt)*sizeof(Gauge18));
}

#if QPHIX_PrecisionInt==1
void 
freeGauge(void *p)
{
    _mm_free(p);
}

void 
freeGauge18(void *p)
{
    _mm_free(p);
}
#endif

void 
setKS_fptype(void *in, void *su3_v, int x, int y, int z, int t)
{
    KS *s_in = (KS*)in;
    fptype *su3_vec = (fptype*)su3_v;

	int x1, x2, y1, y2, z1, z2, t1, t2;
	x1 = (x-Lsxh) / Vxh; x2 = x % Vxh;
	y1 = (y-Lsy) / Vy; y2 = y % Vy;
	z1 = (z-Lsz) / Vz; z2 = z % Vz;
	t1 = (t-Lst) / Vt; t2 = t % Vt;

    int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
	int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

    for(int c = 0; c < 3; c++) {
        s_in[ind][c][0][v] = su3_vec[2*c];
        s_in[ind][c][1][v] = su3_vec[2*c+1];
    }
}

void 
getKS_fptype(void *out, void *su3_v, int x, int y, int z, int t)
{
    KS *s_out = (KS*)out;
    fptype *su3_vec = (fptype*)su3_v;
	int x1, x2, y1, y2, z1, z2, t1, t2;
	x1 = (x-Lsxh) / Vxh; x2 = x % Vxh;
	y1 = (y-Lsy) / Vy; y2 = y % Vy;
	z1 = (z-Lsz) / Vz; z2 = z % Vz;
	t1 = (t-Lst) / Vt; t2 = t % Vt;

    int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
	int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

    for(int c = 0; c < 3; c++) {
        su3_vec[2*c] = s_out[ind][c][0][v];
        su3_vec[2*c+1] = s_out[ind][c][1][v];
    }
}

void 
setFullGauge18_fptype(void *in, void *su3_m, int len)
{
    Gauge18 *g_in = (Gauge18*)in;
    fptype *su3_mat = (fptype*)su3_m;
    const int nrows = 3;
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads) collapse (4)
#endif
	for(int t1 = 0; t1 < Vt; t1++) {
		for(int z1 = 0; z1 < Vz; z1++) {
			for(int y1 = 0; y1 < Vy; y1++) {
				for(int x1 = 0; x1 < Vxh; x1++) {
					int ind = t1*Pxyz+z1*Pxy+y1*Vxh+x1;
#pragma unroll
					for(int d = 0; d < 8; d++) {
#pragma unroll
						for(int c1 = 0; c1 < 3; c1++) {
#pragma unroll
							for(int c2 = 0; c2 < 3; c2++) {
								for(int t2 = 0; t2 < nGT; t2++) {
									for(int z2 = 0; z2 < nGZ; z2++) {
										for(int y2 = 0; y2 < nGY; y2++) {
											for(int x2 = 0; x2 < nGX ; x2++) {
												int v = ((t2*nGZ+z2)*nGY+y2)*nGX+x2;
			                                    int ind2 = ((((t2*Vt+t1)*Nz+(z2*Vz+z1))*Ny+(y2*Vy+y1))*Nxh+(x2*Vxh+x1))*8+d;

												g_in[ind][d][c1][c2][0][v] = su3_mat[ind2*18+2*(3*c1+c2)+0];
												g_in[ind][d][c1][c2][1][v] = su3_mat[ind2*18+2*(3*c1+c2)+1];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void 
setFullGauge_fptype(void *in, void *su3_m, int len)
{
    Gauge *g_in = (Gauge*)in;
    fptype *su3_mat = (fptype*)su3_m;
#ifdef COMPRESSED_12
    const int nrows = 2;
#else
    const int nrows = 3;
#endif
#ifdef _OPENMP
#pragma omp parallel for default(shared) num_threads(nThreads) collapse (4)
#endif
	for(int t1 = 0; t1 < Vt; t1++) {
		for(int z1 = 0; z1 < Vz; z1++) {
			for(int y1 = 0; y1 < Vy; y1++) {
				for(int x1 = 0; x1 < Vxh; x1++) {
					int ind = t1*Pxyz+z1*Pxy+y1*Vxh+x1;
#pragma unroll
					for(int d = 0; d < 8; d++) {
#pragma unroll
						for(int c1 = 0; c1 < nrows; c1++) {
#pragma unroll
							for(int c2 = 0; c2 < 3; c2++) {
								for(int t2 = 0; t2 < nGT; t2++) {
									for(int z2 = 0; z2 < nGZ; z2++) {
										for(int y2 = 0; y2 < nGY; y2++) {
											for(int x2 = 0; x2 < nGX ; x2++) {
												int v = ((t2*nGZ+z2)*nGY+y2)*nGX+x2;
			                                    int ind2 = ((((t2*Vt+t1)*Nz+(z2*Vz+z1))*Ny+(y2*Vy+y1))*Nxh+(x2*Vxh+x1))*8+d;

												g_in[ind][d][c1][c2][0][v] = su3_mat[ind2*18+2*(3*c1+c2)+0];
												g_in[ind][d][c1][c2][1][v] = su3_mat[ind2*18+2*(3*c1+c2)+1];
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void 
setGauge_fptype(void *in, void *su3_m, int dir, int x, int y, int z, int t)
{
    Gauge *g_in = (Gauge*)in;
    fptype *su3_mat = (fptype*)su3_m;
	int x1, x2, y1, y2, z1, z2, t1, t2;
	x1 = (x-Lsxh) / Vxh; x2 = x % Vxh;
	y1 = (y-Lsy) / Vy; y2 = y % Vy;
	z1 = (z-Lsz) / Vz; z2 = z % Vz;
	t1 = (t-Lst) / Vt; t2 = t % Vt;

    int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
	int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

#ifdef COMPRESSED_12
    const int nrows = 2;
#else
    const int nrows = 3;
#endif

    for(int c1 = 0; c1 < nrows; c1++) {
        for(int c2 = 0; c2 < 3; c2++) {
            g_in[ind][dir][c1][c2][0][v] = su3_mat[2*(3*c1+c2)+0];
            g_in[ind][dir][c1][c2][1][v] = su3_mat[2*(3*c1+c2)+1];
        }
    }
}

void getGauge_fptype(void *out, void *su3_m, int dir, int x, int y, int z, int t)
{
    Gauge *g_out = (Gauge*)out;
    fptype *su3_mat = (fptype*)su3_m;
	int x1, x2, y1, y2, z1, z2, t1, t2;
	x1 = (x-Lsxh) / Vxh; x2 = x % Vxh;
	y1 = (y-Lsy) / Vy; y2 = y % Vy;
	z1 = (z-Lsz) / Vz; z2 = z % Vz;
	t1 = (t-Lst) / Vt; t2 = t % Vt;

    int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
	int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

#ifdef COMPRESSED_12
    const int nrows = 2;
#else
    const int nrows = 3;
#endif
    for(int c1 = 0; c1 < nrows; c1++) {
        for(int c2 = 0; c2 < 3; c2++) {
            su3_mat[2*(3*c1+c2)+0] = g_out[ind][dir][c1][c2][0][v];
            su3_mat[2*(3*c1+c2)+1] = g_out[ind][dir][c1][c2][1][v];
        }
    }
    if(nrows == 2)
        for(int c =0 ; c < 3; c++) {
            Conj_CrossProd(&su3_mat[2*(3*2+c)], 
                    &su3_mat[2*(3*0+((c+1)%3))], &su3_mat[2*(3*1+((c+2)%3))], 
                    &su3_mat[2*(3*0+((c+2)%3))], &su3_mat[2*(3*1+((c+1)%3))]);
        }
}

void setGauge18_fptype(void *in, void *su3_m, int dir, int x, int y, int z, int t)
{
    Gauge18 *g_in = (Gauge18*)in;
    fptype *su3_mat = (fptype*)su3_m;
	int x1, x2, y1, y2, z1, z2, t1, t2;
	x1 = (x-Lsxh) / Vxh; x2 = x % Vxh;
	y1 = (y-Lsy) / Vy; y2 = y % Vy;
	z1 = (z-Lsz) / Vz; z2 = z % Vz;
	t1 = (t-Lst) / Vt; t2 = t % Vt;

    int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
	int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

    const int nrows = 3;

    //printf("\n");
    for(int c1 = 0; c1 < 3; c1++) {
        for(int c2 = 0; c2 < nrows; c2++) {
            g_in[ind][dir][c1][c2][0][v] = su3_mat[2*(3*c1+c2)+0];
            g_in[ind][dir][c1][c2][1][v] = su3_mat[2*(3*c1+c2)+1];
        /*    printf("XXX dir=%d c1=%d c2=%d {%g,%g}\n"
                    , dir, c1, c2, g_in[ind][dir][c1][c2][0][x2]
                    , g_in[ind][dir][c1][c2][1][x2]);
        */
        }
    }
}

void 
getGauge18_fptype(void *out, void *su3_m, int dir, int x, int y, int z, int t)
{
    Gauge18 *g_out = (Gauge18*)out;
    fptype *su3_mat = (fptype*)su3_m;
	int x1, x2, y1, y2, z1, z2, t1, t2;
	x1 = (x-Lsxh) / Vxh; x2 = x % Vxh;
	y1 = (y-Lsy) / Vy; y2 = y % Vy;
	z1 = (z-Lsz) / Vz; z2 = z % Vz;
	t1 = (t-Lst) / Vt; t2 = t % Vt;

    int ind = t2*Pxyz+z2*Pxy+y2*Vxh+x2;
	int v = ((t1*nGZ+z1)*nGY+y1)*nGX+x1;

    const int nrows = 3;

    for(int c1 = 0; c1 < 3; c1++) {
        for(int c2 = 0; c2 < nrows; c2++) {
            su3_mat[2*(3*c1+c2)+0] = g_out[ind][dir][c1][c2][0][v];
            su3_mat[2*(3*c1+c2)+1] = g_out[ind][dir][c1][c2][1][v];
        }
    }
}

void 
qphix_ks_dslash_fptype(void *s_in, void *gll_, void *gfl_, void *s_out_, int cb)
{
/*
    printf("qphix_ks_dslash_fptype: VECLEN = %d\n", VECLEN);
    printf("num_floats_in_ks_array = %d\n", num_floats_in_ks_array);
    printf("Vxh = %d Vy = %d Vz = %d Vt = %d Pxy = %d Pxyz = %d\n", Vxh, Vy, Vz, Vt, Pxy, Pxyz);
#if QPHIX_PrecisionInt==1
                printf("QPHIX_PrecisionInt==1\n");
#else
                printf("QPHIX_PrecisionInt==2\n");
#endif
*/
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
        KS *si = (KS*)s_in;
        KS *so = (KS*)s_out_;
        Gauge *gll = (Gauge*)gll_;
        Gauge18 *gfl = (Gauge18*)gfl_;

        // Get my thread ID
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        // accumulate[] is a flag per direction
		//unsigned int accumulate[8], accumulate3[8], isBoundary[8], isBoundary3[8];
		/* = { 1, 1, 1, 1, 1, 1, 1, 1 };
		unsigned int accumulate3[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
		unsigned int isBoundary[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
		unsigned int isBoundary3[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
		*/
		unsigned long long int t0, t1;
		pack_and_send_boundaries(tid, si, cb);
		t0 = __rdtsc();
		Phaser::PhaserState ps;
		int x, y, z, t;
		int x_next, y_next, z_next, t_next;
		//bool loop_ret = phaser_dslash->start(ps, tid/n_threads_per_core, x, y, z, t);
		bool loop_ret = phaser->start(ps, tid, x, y, z, t);
		const char *BTNow = BoundTable + cb*PadBound;
		const unsigned int *NTNow = NeighTable+cb*PadNeigh;
		int Bbytes=8*4/8;
		while(loop_ret) {
			//my_printf("XXX Tid %3d: xyzt = %2d %2d %2d %2d\n", tid, x, y, z, t);
			//loop_ret = phaser_dslash->next(ps, x_next, y_next, z_next, t_next);
			loop_ret = phaser->next(ps, x_next, y_next, z_next, t_next);
			const int xodd = (y + z + t + cb) & 1;
			KS *neighs[8];
			KS *neighs3[8];
			int ind = t*Pxyz+z*Pxy+y*Vxh+x;
			#pragma noprefetch
			for(int jj=0/*, nbit=1*/; jj<8; ++jj)
			{
			    neighs[jj] = &si[NTNow[16*ind+jj]];
			    neighs3[jj] = &si[NTNow[16*ind+8+jj]];
/*
			    accumulate[jj] = BTNow[Bbytes*ind] & nbit;
			    accumulate3[jj] = BTNow[Bbytes*ind+1] & nbit;
			    isBoundary[jj] = BTNow[Bbytes*ind+2] & nbit;
			    isBoundary3[jj] = BTNow[Bbytes*ind+3] & nbit;
			    nbit *= 2;
*/
			}
			
			//for(int ii=1;ii<8;ii++) accumulate[ii] = 0;
			int ind_next = (t_next - t)*Pxyz+(z_next-z)*Pxy+(y_next-y)*Vxh+(x_next-x);

			_mm_prefetch( (char *)&NTNow[16*(ind+ind_next)], _MM_HINT_T0 );
			_mm_prefetch( BTNow+Bbytes*(ind_next+ind), _MM_HINT_T0 );
			//if(tid == 0) my_printf("XXXXX %d %d %d %d\n", t, z, y, x);
			ks_long_dslash_vec_noinline(
			n_threads_per_core,
			&so[ind],
			&gfl[ind],
			&gll[ind],
			neighs,
			neighs3,
			BTNow[Bbytes*ind]/*accumulate*/,
			BTNow[Bbytes*ind+1]/*accumulate3*/,
			BTNow[Bbytes*ind+2]/*isBoundary*/,
			BTNow[Bbytes*ind+3]/*isBoundary3*/,
			ind_next,
			xodd, y, t);

#ifdef ENABLE_BARRIER
			//if(t != t_next &&  ps.ct % BARRIER_FREQ_T == 0 ) phaser->sync(ps); //barriers[ph1][cid_t[ph1]]->wait(group_tid[ph1]);
#endif
			x = x_next;
			y = y_next;
			z = z_next;
			t = t_next;
/*
			_mm_prefetch( BTNow+Bbytes*(ind_next+ind), _MM_HINT_T0 );
			for(int jj=0; jj<8; ++jj) {
			    _mm_prefetch( (char *)&si[NTNow[16*(ind+ind_next)+jj]], _MM_HINT_T1 );
			}
*/
		}; // while(loop_ret == true);
		recv_and_unpack_boundaries(tid, so, gfl, gll, -1.0, cb);
    } // End of OMP Parallel BEGIN SCOPE
}
