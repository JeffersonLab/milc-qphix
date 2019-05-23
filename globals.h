#ifndef _GLOBALS_H_
#define _GLOBALS_H_

// These sets the default values for lattice

#ifndef NX
#define NX                32
#endif

#ifndef NY
#define NY                32
#endif

#ifndef NZ
#define NZ                32
#endif

#ifndef NT
#define NT				  96
#endif

#ifndef BY
//#define BY                4
#define BY		  2
#endif

#ifndef BZ
//#define BZ                4
#define BZ		  1
#endif

// 2 D Threading
#ifndef NCORES
#define NCORES               60
#endif

#ifndef MINCT
#define MINCT              1
#endif

// SMT Threading
#ifndef SY
#define SY      1
#endif

#ifndef SZ
#define SZ      4
#endif

#ifndef XY_PAD
#define XY_PAD (0)
#endif

#ifndef XYZ_PAD
#define XYZ_PAD (0)
#endif
// Nrepeat is the number of iterations to run for a single timing
#ifndef NREPEAT
#define NREPEAT           4
#endif

// Bench Repeat is to repeat the benchmark BENCH_REPEAT times
// Sometimes the second - 3rd etc run has higher perf than the
// first as the caches get 'warmed up'
#ifndef BENCH_REPEAT
#define BENCH_REPEAT	  4
#endif

#ifndef ENABLE_BARRIER
// No barrier
#undef BARRIER_FREQ_T
#else //ENABLE_BARRIER
// Use barrier
// If frequency is not set, set it to something really big so it
// is not called
#ifndef BARRIER_FREQ_T
#define BARRIER_FREQ_T   1024
#endif

#endif // ENABLE_BARRIER

#define my_printf(x...) if(myRank == 0) printf(x)

int Gx=NX, Gy=NY, Gz=NZ, Gt=NT;
int Nxh, Nx, Ny, Nz, Nt;
int Vxh, Vx, Vy, Vz, Vt;
int Vxhf, Vxf, Vyf, Vzf, Vtf;
int Vxhd, Vxd, Vyd, Vzd, Vtd;
int Lsxh = 0, Lsx = 0, Lsy = 0, Lsz = 0, Lst = 0;
int Bx = 0, By = BY, Bz=BZ;
int Sy = SY, Sz = SZ, n_threads_per_core;
int Pxy, Pxyz;
int Pxyf, Pxyzf;
int Pxyd, Pxyzd;
int NCores = NCORES;
int iters = NREPEAT;
int bench_repeats = BENCH_REPEAT;
int minCt = MINCT;
Barrier* gBar;

int geometry[4] = {1, 1, 1, 1};
int myRank = 0;
int nRanks = 1;
bool local_dir[4] = {true, true, true, true};

int xy_pad = XY_PAD;
int xyz_pad = XYZ_PAD;

bool dwait = false, dump = false;
int waitrank=-1;
bool printAllRanks=false;
int nParallelRecvs = 1, nParallelSends = 1;

#endif // _GLOBALS_H_
