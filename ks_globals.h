#ifndef _KS_GLOABLS_H_
#define _KS_GLOABLS_H_

#define my_printf(x...) if(myRank == 0) printf(x)

/* Needed header files */
#include "ks_config.h"
#include "qcd_data_types.h"
#include "Barrier.h"
#include "phaser.h"
//#include "ubench_helper.h"

/* Global variables for the lattice dimensions and number of threads etc. */
extern int Gxh, Gx, Gy, Gz, Gt;
extern int Nxh, Nx, Ny, Nz, Nt;
extern int Lsxh, Lsx, Lsy, Lsz, Lst;
extern int Vxh, Vx, Vy, Vz, Vt;
extern int Pxy, Pxyz;
extern int NCores;
extern int n_threads_per_core;
extern int nThreads;
extern int nThreads_dslash;

extern int num_floats_in_ks_array;
extern int num_floats_in_gauge_array;
extern int num_floats_in_hermit_array;
/* Global variables for the offset calculations */
extern Phaser *phaser;
//extern Phaser *phaser_dslash;
extern Barrier* gBar;
extern bool local_dir[4];
extern int geometry[4];

extern int myRank;
extern int nRanks;

#endif
