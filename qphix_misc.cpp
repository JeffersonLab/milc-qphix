/******************************************/
/* QPHIX initialize and finalize routines */
/******************************************/

#include <string.h>
#include <iostream>
#include <iomanip>
#include <cassert>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ks_config.h"
#include "ks_globals.h"        /* All global values and macros */
#include "ks_long_dslash.h"
#include "ks_boundary.h"
#include "gf_boundary.h"
#include "misc.h"

#include "qphix_internal.h"

#if QPHIX_PrecisionInt == 1
int qphix_sites_on_node;
int qphix_even_sites_on_node;
//int minCt;
#else
extern int qphix_sites_on_node;
extern int qphix_even_sites_on_node;
#endif

extern char * BoundTableF, * BoundTableD;
extern unsigned int * NeighTableF, * NeighTableD;
//unsigned int * Neigh3Table = 0x00;
extern int BLENGTHF, BLENGTHD;
extern int PadBoundF, PadBoundD;
extern int PadNeighF, PadNeighD;
//int NLENGTH = 0;

/* Gauge force */
#ifdef QPHIX_GF
extern char * GFBoundTable;
extern unsigned int * GFNeighTable;
extern int GFBLENGTH;
#endif

/* Set up ks dslash boundary&neighbor tables */
static void init_ks_dslash_imp()
{
    int Bbytes=8*4/8;
    int Lsize=Pxyz*Vt*2;
    int Nbytes=8*2*sizeof(unsigned int);
    int Vbytes=VECLEN*sizeof(fptype);
    MYASSERT((BLENGTH = Vbytes/Bbytes)>0);
    if(BoundTable==0x00)
	MYASSERT((BoundTable = (char *)_mm_malloc(Lsize*Bbytes,64))!=0x00);
    if(NeighTable==0x00)
	MYASSERT((NeighTable = (unsigned int *)_mm_malloc(Lsize*Nbytes,64))!=0x00);
/*    if(Neigh3Table==0x00)
	MYASSERT((Neigh3Table = _mm_malloc(Lsize*Nbytes))!=0x00);
*/
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        // accumulate[] is a flag per direction
        Phaser::PhaserState ps;
	int x, y, z, t, x_next, y_next, z_next, t_next;
	PadBound=Lsize/2*Bbytes;
	PadNeigh=Lsize/2*8*2;
	for(int cb=0; cb<2; ++cb)
	{
	  bool loop_ret = phaser->start(ps, tid, x, y, z, t);
	  int PBnow = cb*PadBound;
	  int PNnow = cb*PadNeigh;
	  while(loop_ret) {
	    loop_ret = phaser->next(ps, x_next, y_next, z_next, t_next);
	    const int xodd = (y + z + t + cb) & 1;
	    int ind = t*Pxyz+z*Pxy+y*Vxh+x;
	    /* accumulate */
	    BoundTable[PBnow+Bbytes*ind] = ~0;
	    /* accumulate3 */
	    BoundTable[PBnow+Bbytes*ind+1] = ~0;
	    /* isBoundary */
	    BoundTable[PBnow+Bbytes*ind+2] = 0;
	    /* isBoundary3 */
	    BoundTable[PBnow+Bbytes*ind+3] = 0;
	    NeighTable[PNnow+16*ind] = ind+xodd-1;
	    if(x+xodd-1 < 0) {
		NeighTable[PNnow+16*ind] += Vxh;
		if(!local_dir[0]) BoundTable[PBnow+Bbytes*ind] ^= 1;
		BoundTable[PBnow+Bbytes*ind+2] ^= 1;
	    }

	    NeighTable[PNnow+16*ind+1] = ind+xodd;
	    if(x+xodd >= Vxh) {
		NeighTable[PNnow+16*ind+1] -= Vxh;
		if(!local_dir[0]) BoundTable[PBnow+Bbytes*ind] ^= 2;
		BoundTable[PBnow+Bbytes*ind+2] ^= 2;
	    }

	    NeighTable[PNnow+16*ind+2] = ind-Vxh;
	    if(y == 0) {
		NeighTable[PNnow+16*ind+2] += Vy*Vxh;
		if(!local_dir[1]) BoundTable[PBnow+Bbytes*ind] ^= 4;
		BoundTable[PBnow+Bbytes*ind+2] ^= 4;
	    }

	    NeighTable[PNnow+16*ind+3] = ind+Vxh;
	    if(y == Vy-1) {
		NeighTable[PNnow+16*ind+3] -= Vy*Vxh;
		if(!local_dir[1]) BoundTable[PBnow+Bbytes*ind] ^= 8;
		BoundTable[PBnow+Bbytes*ind+2] ^= 8;
	    }

	    NeighTable[PNnow+16*ind+4] = ind-Pxy;
	    if(z == 0) {
		NeighTable[PNnow+16*ind+4] += Vz*Pxy;
		if(!local_dir[2]) BoundTable[PBnow+Bbytes*ind] ^= 16;
		BoundTable[PBnow+Bbytes*ind+2] ^= 16;
	    }

	    NeighTable[PNnow+16*ind+5] = ind+Pxy;
	    if(z == Vz-1) {
		NeighTable[PNnow+16*ind+5] -= Vz*Pxy;
		if(!local_dir[2]) BoundTable[PBnow+Bbytes*ind] ^= 32;
		BoundTable[PBnow+Bbytes*ind+2] ^= 32;
	    }
	    NeighTable[PNnow+16*ind+6] = ind-Pxyz;
	    if(t == 0) {
		NeighTable[PNnow+16*ind+6] += Vt*Pxyz;
		if(!local_dir[3]) BoundTable[PBnow+Bbytes*ind] ^= 64;
		BoundTable[PBnow+Bbytes*ind+2] ^= 64;
	    }

	    NeighTable[PNnow+16*ind+7] = ind+Pxyz;
	    if(t == Vt-1) {
		NeighTable[PNnow+16*ind+7] -= Vt*Pxyz;
		if(!local_dir[3]) BoundTable[PBnow+Bbytes*ind] ^= 128;
		BoundTable[PBnow+Bbytes*ind+2] ^= 128;
	    }

	    NeighTable[PNnow+16*ind+8] = ind+xodd-2;
	    if(x+xodd-2 < 0) {
		NeighTable[PNnow+16*ind+8] += Vxh;
		if(!local_dir[0]) BoundTable[PBnow+Bbytes*ind+1] ^= 1;
		if(x+xodd-2+Vxh < 0) 
		    NeighTable[PNnow+16*ind+8] += Vxh;
		else
		    BoundTable[PBnow+Bbytes*ind+3] ^= 1;
	    }

	    NeighTable[PNnow+16*ind+9] = ind+xodd+1;
	    if(x+xodd+1 >= Vxh) {
		NeighTable[PNnow+16*ind+9] -= Vxh;
		if(!local_dir[0]) BoundTable[PBnow+Bbytes*ind+1] ^= 2;
		if(x+xodd+1-Vxh >= Vxh) 
		    NeighTable[PNnow+16*ind+9] -= Vxh;
		else
		    BoundTable[PBnow+Bbytes*ind+3] ^= 2;
	    }

	    NeighTable[PNnow+16*ind+10] = ind-3*Vxh;
	    if(y-3 < 0) {
		NeighTable[PNnow+16*ind+10] += Vy*Vxh;
		if(!local_dir[1]) BoundTable[PBnow+Bbytes*ind+1] ^= 4;
		if(y-3+Vy < 0) 
		    NeighTable[PNnow+16*ind+10] += Vy*Vxh;
		else
		    BoundTable[PBnow+Bbytes*ind+3] ^= 4;
	    }

	    NeighTable[PNnow+16*ind+11] = ind+3*Vxh;
	    if(y+3 >= Vy) {
		NeighTable[PNnow+16*ind+11] -= Vy*Vxh;
		if(!local_dir[1]) BoundTable[PBnow+Bbytes*ind+1] ^= 8;
		if(y+3-Vy >= Vy) 
		    NeighTable[PNnow+16*ind+11] -= Vy*Vxh;
		else
		    BoundTable[PBnow+Bbytes*ind+3] ^= 8;
	    }

	    NeighTable[PNnow+16*ind+12] = ind-3*Pxy;
	    if(z-3 < 0) {
		NeighTable[PNnow+16*ind+12] += Vz*Pxy;
		if(!local_dir[2]) BoundTable[PBnow+Bbytes*ind+1] ^= 16;
		if(z-3+Vz < 0) 
		    NeighTable[PNnow+16*ind+12] += Vz*Pxy;
		else
		    BoundTable[PBnow+Bbytes*ind+3] ^= 16;
	    }

	    NeighTable[PNnow+16*ind+13] = ind+3*Pxy;
	    if(z+3 >= Vz) {
		NeighTable[PNnow+16*ind+13] -= Vz*Pxy;
		if(!local_dir[2]) BoundTable[PBnow+Bbytes*ind+1] ^= 32;
		if(z+3-Vz >= Vz) 
		    NeighTable[PNnow+16*ind+13] -= Vz*Pxy;
		else
		    BoundTable[PBnow+Bbytes*ind+3] ^= 32;
	    }

	    NeighTable[PNnow+16*ind+14] = ind-3*Pxyz;
	    if(t-3 < 0) {
		NeighTable[PNnow+16*ind+14]+= Vt*Pxyz;
		if(!local_dir[3]) BoundTable[PBnow+Bbytes*ind+1] ^= 64;
		if(t-3+Vt < 0) 
		    NeighTable[PNnow+16*ind+14]+= Vt*Pxyz;
		else
		    BoundTable[PBnow+Bbytes*ind+3] ^= 64;
	    }

	    NeighTable[PNnow+16*ind+15] = ind+3*Pxyz;
	    if(t+3 >= Vt) {
		NeighTable[PNnow+16*ind+15] -= Vt*Pxyz;
		if(!local_dir[3]) BoundTable[PBnow+Bbytes*ind+1] ^= 128;
		if(t+3-Vt >= Vt) 
		    NeighTable[PNnow+16*ind+15] -= Vt*Pxyz;
		else
		    BoundTable[PBnow+Bbytes*ind+3] ^= 128;
	    }

            x=x_next;
            y=y_next;
            z=z_next;
            t=t_next;
	}
      }
    } // end omp loop
}

#if 0
#ifdef QPHIX_GF
static void init_gauge_force( )
{
    int Bbytes=8*2/8;
    int Lsize=Pxyz*Vt;
    int Nbytes=8*2*sizeof(unsigned int);
    int Vbytes=VECLEN*sizeof(fptype);
    MYASSERT((GFBLENGTH = Vbytes/Bbytes)>0);
    if(GFBoundTable==0x00)
        MYASSERT((GFBoundTable = (char *)_mm_malloc(Lsize*Bbytes,64))!=0x00);
    if(GFNeighTable==0x00)
        MYASSERT((GFNeighTable = NeighTable)!=0x00);
    int otherparity=1;
    if(GF_PARITY==QPHIX_ODD) otherparity=0;
    else MYASSERT(GF_PARITY==QPHIX_EVEN || GF_PARITY==QPHIX_EVENODD);
#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        Phaser::PhaserState ps;
        int x, y, z, t, x_next, y_next, z_next, t_next;
        cb = 0;
        {
          bool loop_ret = phaser->start(ps, tid, x, y, z, t);
          int PBnow = 0; //cb*PadBound;
          while(loop_ret) {
            loop_ret = phaser->next(ps, x_next, y_next, z_next, t_next);
            const int xodd = (y + z + t + otherparity) & 1;
            int ind = t*Pxyz+z*Pxy+y*Vxh+x;
            /* accumulate */
            GFBoundTable[PBnow+Bbytes*ind] = ~0;
            /* isBoundary */
            GFBoundTable[PBnow+Bbytes*ind+1] = 0;
            if(!local_dir[0]) {
              if(x+xodd-2 < 0) {
                GFBoundTable[PBnow+Bbytes*ind] ^= 1;
		GFBoundTable[PBnow+Bbytes*ind+1] ^= 1;
              }
              if(x+xodd+1 >= Vxh) {
                GFBoundTable[PBnow+Bbytes*ind] ^= 2;
		GFBoundTable[PBnow+Bbytes*ind+1] ^= 2;
              }
	    }
	    else {
	      if(x+xodd-1 < 0) GFBoundTable[PBnow+Bbytes*ind+1] ^= 1;
	      if(x+xodd >= Vxh) GFBoundTable[PBnow+Bbytes*ind+1] ^= 2;
	    }

	    if(!local_dir[1]) {
              if(y <= 2) {
                GFBoundTable[PBnow+Bbytes*ind] ^= 4;
                GFBoundTable[PBnow+Bbytes*ind+1] ^= 4;
              }
              if(y >= Vy-3) {
                GFBoundTable[PBnow+Bbytes*ind] ^= 8;
                GFBoundTable[PBnow+Bbytes*ind+1] ^= 8;
              }
	    }
	    else {
	      if(y == 0) GFBoundTable[PBnow+Bbytes*ind+1] ^= 4;
	      if(y == Vy-1) GFBoundTable[PBnow+Bbytes*ind+1] ^= 8;
	    }

            if(!local_dir[2]) {
              if(z <= 2) {
                GFBoundTable[PBnow+Bbytes*ind] ^= 16;
                GFBoundTable[PBnow+Bbytes*ind+1] ^= 16;
              }
	      if(z >= Vz-3) {
                GFBoundTable[PBnow+Bbytes*ind] ^= 32;
                GFBoundTable[PBnow+Bbytes*ind+1] ^= 32;
              }
	    }
	    else {
	      if(z == 0) GFBoundTable[PBnow+Bbytes*ind+1] ^= 16;
	      if(z == Vz-1) GFBoundTable[PBnow+Bbytes*ind+1] ^= 32;
	    }

            if(!local_dir[3]) {
              if(t <= 2) {
                GFBoundTable[PBnow+Bbytes*ind] ^= 64;
                GFBoundTable[PBnow+Bbytes*ind+1] ^= 64;
              }
              if(t >= Vt-3) {
                GFBoundTable[PBnow+Bbytes*ind] ^= 128;
                GFBoundTable[PBnow+Bbytes*ind+1] ^= 128;
              }
	    }
	    else {
	      if(t == 0) GFBoundTable[PBnow+Bbytes*ind+1] ^= 64;
	      if(t == Vt-1) GFBoundTable[PBnow+Bbytes*ind+1] ^= 128;
	    }

            x=x_next;
            y=y_next;
            z=z_next;
            t=t_next;
        }
      }
    } // end omp loop

}

#endif
#endif

/* Convert rank to coordinates */
static void lex_coords(int coords[], const int dim, const int size[], const size_t rank)
{
  int d;
  size_t r = rank;

  for(d = 0; d < dim; d++){
    MYASSERT(size[d]>0);
    coords[d] = r % size[d];
    r /= size[d];
  }
}

static void neigh_lex_coords(int coords[], const int latdim, const int size[], const size_t rank, int (*node_number)(const int coord[]))
{
  int d;
  const int orig[4] = {0,0,0,0};
  const int xend[4] = {Gx-1,0,0,0};
  const int yend[4] = {Gx-1,Gy-1,0,0};
  const int zend[4] = {Gx-1,Gy-1,Gz-1,0};
  const int tend[4] = {Gx-1,Gy-1,Gz-1,Gt-1};
  int nblock[4];
  nblock[0] = 1 + node_number(xend) - node_number(orig);
  MYASSERT(nblock[0]>0);
  nblock[1] = 1 + node_number(yend) - node_number(orig);
  MYASSERT(nblock[1]>0 && nblock[1]%nblock[0]==0);
  nblock[2] = 1 + node_number(zend) - node_number(orig);
  MYASSERT(nblock[2]>0 && nblock[2]%nblock[1]==0);
  nblock[3] = 1 + node_number(tend) - node_number(orig);
  MYASSERT(nblock[3]>0 && nblock[3]%nblock[2]==0);
  nblock[3] /= nblock[2];
  nblock[2] /= nblock[1];
  nblock[1] /= nblock[0];

  int inc = (rank%nblock[0] ? 0 : nblock[0]);
  coords[0] = rank - 1 + inc;
  inc = ((rank+1)%nblock[0] ? 0 : nblock[0]);
  coords[1] = rank + 1 - inc;
  int bsize = 1;
  MYASSERT(latdim<=4);
  int rankn = rank;
  for(d=1; d<latdim; ++d) {
    bsize *= nblock[d-1];
    rankn /= nblock[d-1];
    inc = (rankn%nblock[d] ? 0 : nblock[d]);
    coords[2*d] = rank - bsize * (1 - inc);
    inc = ((rankn+1)%nblock[d] ? 0 : nblock[d]);
    coords[2*d+1] = rank + bsize * (1 - inc);
  }

}

#if QPHIX_PrecisionInt==1
#define QPHIX_init_fptype QPHIX_init_F
#define setup_comms setup_comms_F
#define gf_setup_comms gf_setup_comms_F
#define QPHIX_finalize_fptype QPHIX_finalize_F
#define destroy_comms destroy_comms_F
#define gf_destroy_comms gf_destroy_comms_F
#elif QPHIX_PrecisionInt==2
#define QPHIX_init_fptype QPHIX_init_D
#define setup_comms setup_comms_D
#define gf_setup_comms gf_setup_comms_D
#define QPHIX_finalize_fptype QPHIX_finalize_D
#define destroy_comms destroy_comms_D
#define gf_destroy_comms gf_destroy_comms_D
#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif

QPHIX_status_t 
QPHIX_init_fptype(QPHIX_layout_t *layout)
{
    int sy =1, sz = 1;
    Gx = layout->latsize[0];
    Gxh = Gx/2;
    Gy = layout->latsize[1];
    Gz = layout->latsize[2];
    Gt = layout->latsize[3];

    Lsxh = Lsx = Lsy = Lsz = Lst = 0;
    qphix_even_sites_on_node = layout->even_sites_on_node;
    //printf("qphix_even_sites_on_node = %d\n", qphix_even_sites_on_node);
    qphix_sites_on_node = layout->sites_on_node;

    int minCt = 1;
    if(getenv("MINCT")) minCt = atoi(getenv("MINCT"));
    n_threads_per_core = 1;
#ifdef _OPENMP
    if(getenv("THREADS_PER_CORE")) n_threads_per_core = atoi(getenv("THREADS_PER_CORE"));
    int numThreads = omp_get_max_threads();
    NCores = numThreads / n_threads_per_core;
#endif
    myRank = layout->this_node;
    printAllRanks = false;
    nRanks = 1;
    for(int i=0; i<4; ++i)
	local_dir[i] = true;
    if(layout->machdim != layout->latdim) 
    {
	if(myRank==0) cout << "QPHIX ERROR: Right now only takes Machine dimension same as Lattice dimension" << endl;
	return QPHIX_FAIL;
    }
    for(int i = 0; i < layout->machdim; i++)	
    {
	geometry[i] = layout->machsize[i];
	nRanks *= geometry[i];
    }
#ifndef ENABLE_MPI
    if(nRanks != 1) {
	printf("MPI not enabled in QPhiX, please compile with ENABLE_MPI=1\n"); 
	return QPHIX_FAIL;
    }
#else
    MPI_Barrier(MPI_COMM_WORLD);    
#endif

    int m_coord[4], n_coord[8];
    const int mdim = layout->machdim;
    //lex_coords(m_coord, mdim, geometry, (const size_t)myRank);
    lex_coords(m_coord, mdim, geometry, myRank);
    //neigh_lex_coords(n_coord, (const int)layout->latdim, geometry, (const size_t)myRank, layout->node_number);
    neigh_lex_coords(n_coord, layout->latdim, geometry, myRank, layout->node_number);

    setup_comms( m_coord, n_coord );
    gf_setup_comms();

                MYASSERT(Nxh % 2 == 0);
                MYASSERT(Ny % 4 == 0);
                MYASSERT(Nz % 4 == 0);
                MYASSERT(Nt % 4 == 0);

	//printf("QPHIX_init: VECLEN = %d\n", VECLEN);

                Vx = (VECLEN > 1 ? Nx/2 : Nx);
                Vxh = (VECLEN > 1 ? Nxh/2 : Nxh);
                Vy = (VECLEN > 2 ? Ny/2 : Ny);
                Vz = (VECLEN > 4 ? Nz/2 : Nz);
                Vt = (VECLEN > 8 ? Nt/2 : Nt);
                MYASSERT(Vy % BY == 0);
                MYASSERT(Vz % BZ == 0);
                Pxy = (Vxh*Vy+XY_PAD);
                Pxyz = (Pxy*Vz+XYZ_PAD);

	//printf("Vxh = %d Vy = %d Vz = %d Vt = %d\n", Vxh, Vy, Vz, Vt);

    if(n_threads_per_core == 1) sz = 1;
    else if(n_threads_per_core == 2) sz = 2;
    else if(n_threads_per_core == 4) {
        sy = 2;
        sz = 2;
    }
    else {
        printf("Threads per core must be set to 1, 2 or 4. Value is %d\n", n_threads_per_core);
	return QPHIX_FAIL;
    }

    num_floats_in_ks_array = (Pxyz*Vt)* sizeof(KS)/sizeof(fptype);
/*
    num_floats_in_gauge_array = (Pxyz*Vt)* sizeof(Gauge)/sizeof(fptype);
    num_floats_in_hermit_array = (Pxyz*Vt)* sizeof(Hermit)/sizeof(fptype);
*/
                gBar = new Barrier(NCores,n_threads_per_core);
		if(gBar==0x00) return QPHIX_MEM_ERROR;
                phaser = new Phaser(Vxh, Vy, Vz, Vt, Vxh
                      , BY, BZ, NCores, sy, sz, minCt);
/*
		if(phaser==0x00) return QPHIX_MEM_ERROR;
		phaser_dslash = new Phaser(Vxh, Vy, Vz, Vt, Vxh
                      , BY, BZ, NCores, 1, 1, minCt);
		if(phaser_dslash==0x00) return QPHIX_MEM_ERROR;
*/
    nThreads = NCores * n_threads_per_core;

#ifdef _OPENMP
#pragma omp parallel default (shared) num_threads(nThreads)
    {
        int tid = omp_get_thread_num();
        gBar->init(tid);
        phaser->init(tid);
	//if(tid%n_threads_per_core==0) phaser->init(tid);
    }
#else
    if(nThreads != 1) {printf("No. of threads set to %d while not using OpenMP\n", nThreads); }
    phaser->init(0);
    //phaser_dslash->init(0);
#endif

    init_ks_dslash_imp();
/*
#ifdef QPHIX_GF
    init_gauge_force();
#endif
*/
    //printf("Done initializing QPHIX\n");

    return QPHIX_SUCCESS;

}

QPHIX_status_t
QPHIX_finalize_fptype(void)
{
  delete gBar;
  delete phaser;
  gBar=0x00;
  phaser=0x00;
  _mm_free((void *)BoundTable);
  _mm_free((void *)NeighTable);
  BoundTable = 0x00;
  NeighTable = 0x00;
  destroy_comms();
  gf_destroy_comms();
  BLENGTH = 0;
  PadBound = 0;
  PadNeigh = 0;
  Gxh = Gx = Gy = Gz = Gt = 0;
  Nxh = Nx = Ny = Nz = Nt = 0;
  Vxh = Vx = Vy = Vz = Vt = 0;
  Lsxh = Lsx = Lsy = Lsz = Lst = 0;
  Pxy = Pxyz = 0;
  NCores = 1;
  n_threads_per_core = 1;
  nThreads = 1;
  geometry[0] = geometry[1] = geometry[2] = geometry[3] = 1;
  qphix_sites_on_node = qphix_even_sites_on_node = 0;
  num_floats_in_ks_array = 0;

  return QPHIX_SUCCESS;
}

#if QPHIX_PrecisionInt == 1
QPHIX_status_t
//QPHIX_init(QPHIX_layout_t *layout, int prec)
QPHIX_init(QPHIX_layout_t *layout)
{
  return (QPHIX_status_t)(QPHIX_init_F(layout) | QPHIX_init_D(layout));
}

QPHIX_status_t
QPHIX_finalize()
{
  return (QPHIX_status_t)(QPHIX_finalize_F() | QPHIX_finalize_D());
}

#endif

