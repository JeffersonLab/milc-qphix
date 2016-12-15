/***************** gauge_force_imp.cpp *****************/
/*
 *
 */ 
/*******************************************************/

#include "qphix_internal.h"
#include "gauge_force_imp.h"

#define TIME_GF 1
#define CHECK 1

#include <sys/time.h>
#if TIME_GF
#include "hrtimer.h"   /* High resolution timer */
#endif

/*
#if VTUNE
#include <ittnotify.h>
#endif
*/

//#undef QPHIX_fptype_init_gauge_force_imp
//#undef QPHIX_fptype_update_h_u_gauge_force_imp
#undef allocHermit_fptype
#undef allocHermitHelper_fptype
#undef allocHermitHelperYZT_fptype
#undef update_h_u_gforce_imp_fptype
#undef QPHIX_fptype_Force
#undef QPHIX_fptype_GaugeField

#if QPHIX_PrecisionInt==1
#define allocHermit_fptype allocHermit_F
#define allocHermitHelper_fptype allocHermitHelper_F
#define allocHermitHelperYZT_fptype allocHermitHelperYZT_F
#define QPHIX_symanzik_1loop_gauge_force_fptype QPHIX_F3_symanzik_1loop_gauge_force
//#define QPHIX_fptype_init_gauge_force_imp QPHIX_F3_init_gauge_force_imp
//#define QPHIX_fptype_update_h_u_gauge_force_imp QPHIX_F3_update_h_u_gauge_force_imp
#define update_h_u_gforce_imp_fptype update_h_u_gforce_imp_F
#define QPHIX_fptype_Force QPHIX_F3_Force
#define QPHIX_fptype_GaugeField QPHIX_F3_GaugeField
#elif QPHIX_PrecisionInt==2
#define allocHermit_fptype allocHermit_D
#define allocHermitHelper_fptype allocHermitHelper_D
#define allocHermitHelperYZT_fptype allocHermitHelperYZT_D
#define QPHIX_symanzik_1loop_gauge_force_fptype QPHIX_D3_symanzik_1loop_gauge_force
//#define QPHIX_fptype_init_gauge_force_imp QPHIX_D3_init_gauge_force_imp
//#define QPHIX_fptype_update_h_u_gauge_force_imp QPHIX_D3_update_h_u_gauge_force_imp
#define update_h_u_gforce_imp_fptype update_h_u_gforce_imp_D
#define QPHIX_fptype_Force QPHIX_D3_Force
#define QPHIX_fptype_GaugeField QPHIX_D3_GaugeField
#else
#error "QPHIX_PrecisionInt Not defined"
#endif

#if QPHIX_PrecisionInt==1
int num_floats_in_hermit_array;
int num_floats_in_gauge_array;
#endif

const int nGX = (VECLEN < 2 ? 1 : 2);
const int nGY = (VECLEN < 4 ? 1 : 2);
const int nGZ = (VECLEN < 8 ? 1 : 2);
const int nGT = (VECLEN < 16 ? 1 : 2);

#if QPHIX_PrecisionInt==1
const int MILC2MBENCH[8] = { 1, 3, 5, 7, 6, 4, 2, 0 };
const int MBENCH2MILC[8] = { 7, 0, 6, 1, 5, 2, 4, 3 };
#endif

extern double Freq;
extern int PadBound;
extern int PadNeigh;
extern char * BoundTable;
extern unsigned int * NeighTable;
extern int qphix_even_sites_on_node;

static void
*initBuf(void *b, size_t s)
{
        fptype *buf = (fptype*)b;
#pragma omp parallel for num_threads (nThreads)
        for(int i = 0; i < s/sizeof(fptype); i++)
                buf[i] = 0.0f;
        return b;
}

void *allocHermit_fptype()
{
    return initBuf(_mm_malloc((Pxyz*Vt+1)*sizeof(Hermit), 64), (Pxyz*Vt+1)*sizeof(Hermit));
}

void *allocHermitHelper_fptype()
{
    return initBuf(_mm_malloc((Pxyz*Vt+1)*sizeof(HermitHelper), 64), (Pxyz*Vt+1)*sizeof(HermitHelper));
}

void *allocHermitHelperYZT_fptype()
{
    return initBuf(_mm_malloc((Pxyz*Vt+1)*sizeof(HermitHelperYZT), 64), (Pxyz*Vt+1)*sizeof(HermitHelperYZT));
}

/*
void QPHIX_fptype_init_gauge_force_imp( QPHIX_fptype_Force *momout )
{
  momout->parity = GF_PARITY;
  MYASSERT(momout->parity==QPHIX_EVEN || momout->parity==QPHIX_ODD || momout->parity->QPHIX_EVENODD);
  momout->heo = _mm_malloc(Pxyz*Vt*sizeof(Hermit), 64);
  momout->htmp = _mm_malloc(Pxyz*Vt*sizeof(HermitHelperYZT), 64);
  MYASSERT(momout->heo!=0x00);
}

void QPHIX_fptype_finalize_gauge_force_imp( QPHIX_fptype_Force *momout )
{
  if(momout->heo!=0x00)
    delete momout->heo;
  momout->heo=0x00;
  momout->parity = 0;
}
*/
/*
#if QPHIX_PrecisionInt==1
static int MASK[8] = { 0xAAAA, 0x5555, 0xCCCC, 0x3333, 0xF0F0, 0x0F0F, 0xFF00, 0x00FF };
#elif QPHIX_PrecisionInt==2
static int MASK[8] = { 0xAA, 0x55, 0xCC, 0x33, 0xF0, 0x0F, 0x00, 0x00 };
#endif
*/

/*
#if VTUNE
static int count=0;
#endif
*/

void QPHIX_symanzik_1loop_gauge_force_fptype(QPHIX_info_t *info,
                                         QPHIX_fptype_GaugeField *gauge,
                                         QPHIX_fptype_Force *force,
                                         QPHIX_gauge_coeffs_t *coeffs,
                                         fptype eps)
{
  MYASSERT(force->parity == gauge->parity);
  struct timeval tv;
  gettimeofday( &tv, NULL );
  info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec;
#if TIME_GF
    hrtimer_t t = hrtimer_get();
#endif

#ifdef CHECK
#if 0
  printf("node %d : \n", myRank);
  fflush(stdout);
  char fn[10];
  sprintf(fn, "node%d", myRank);
  //Gauge *gi = (Gauge*)gauge->geo;
  Hermit *hi = (Hermit*)force->heo;
  FILE *fp = fopen(fn, "w");
  //for(int i=0; i<Pxyz*Vt; ++i) for(int d=0; d<8; ++d) for(int r=0; r<3; ++r) for(int c=0; c<3; ++c) for(int ir=0; ir<2; ++ir) for(int v=0; v<VECLEN; ++v)
  for(int i=0; i<Pxyz*Vt; ++i) for(int d=0; d<8; ++d) for(int r=0; r<8; ++r) for(int v=0; v<VECLEN; ++v)
    if(hi[i][d][r][v]!=0.) 
    {
	fprintf(fp, "hermit[%d][%d][%d][%d] = %f\n", i+myRank*Pxyz*Vt, d, r, v, hi[i][d][r][v]);
	//printf("gauge[%d][%d][%d][%d][%d][%d] = %f\n", i+myRank*Pxyz*Vt, d, r, c, ir, v, gi[i][d][r][c][ir][v]);
    }
  fflush(fp);
  fclose(fp);
    exit(0);
#endif  
#endif
/*
#if VTUNE
  if(count<10)
  {
    printf("Starting VTune\n");
    __itt_resume();
  }
#endif
*/
  update_h_u_gforce_imp(force->heo, force->htmp, gauge->geo, (gauge->parity==QPHIX_EVEN ? 1 : 0), eps, eps, coeffs->plaquette, coeffs->rectangle, coeffs->parallelogram);
/*
#if VTUNE
  if(count<10)
  {
    printf("Pausing VTune\n");
    __itt_pause();
    count++;
  }
#endif
*/
#if TIME_GF
  //cout << "QPHIX_symanzik_1loop_gauge_force time (sec) : " << (double)(hrtimer_get()-t)/1.e+9 << endl;
#endif  
  gettimeofday( &tv, NULL );
  info->status = QPHIX_SUCCESS;
  info->final_sec = tv.tv_sec + 1.e-6 * tv.tv_usec - info->final_sec;
  info->final_flop = 61092.0; //153004.0;

}

#if 0
/* Momentum and gauge field update from gauge force wrapper */
void QPHIX_fptype_update_h_u_gauge_force_imp( QPHIX_fptype_Force *mom,
                            	   QPHIX_fptype_GaugeField *gauge,
				   double epsilonU, double epsilonH )
{
  MYASSERT(mom->parity & gauge->parity);
  int parity = mom->parity & QPHIX_EVEN; /* parity==0 : ODD */
  if(parity==0) {
    update_h_u_gforce_imp(mom->heo, (gauge->geo), 0, epsilonH, epsilonU, mom->kappaS, mom->kappaR, mom->kappaB); 
  }
  else {
    update_h_u_gforce_imp(mom->heo, (gauge->geo), 1, epsilonH, epsilonU, mom->kappaS, mom->kappaR, mom->kappaB);
  }
}
#endif

/* Calculate improved gauge forces */
/* Update Momentum */
/* Then update Gauge field */
void update_h_u_gforce_imp_fptype( void *h_io, void *h_tmp, void *g_io, int cb, double epsilonH, double epsilonU, double kappaS, double kappaR, double kappaB ) 
{

#ifdef _OPENMP
#pragma omp parallel default(shared) num_threads(nThreads)
#endif
    {
        Hermit *hio = (Hermit*)h_io;
	HermitHelperYZT *htmp = (HermitHelperYZT*)h_tmp;
        Gauge *gio = (Gauge*)g_io;

#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
	/* Send boundary gauge fields */
/*
	if(tid==0) printf("sending u...\n");
	fflush(stdout);
*/
	gf_pack_and_send_boundaries_u(tid, gio, cb);
/*
	if(tid==0) printf("done sending u\n");
	fflush(stdout);
*/
                unsigned long long int t0, t1;
                t0 = __rdtsc();
                Phaser::PhaserState ps;
                int x, y, z, t;
                int x_next, y_next, z_next, t_next;
                bool loop_ret = phaser->start(ps, tid, x, y, z, t);

                const char *BTNow = BoundTable+cb*PadBound;
                const unsigned int *NTNow = NeighTable+cb*PadNeigh;
                int Bbytes=8*4/8;

                while(loop_ret) {
                        loop_ret = phaser->next(ps, x_next, y_next, z_next, t_next);
                        Gauge *g_in_neighs[8];
                        Hermit *h_out_neighs[8];
                        int ind = t*Pxyz+z*Pxy+y*Vxh+x;
                        #pragma noprefetch
                        for(int jj=0/*, nbit=1*/; jj<8; ++jj)
                        {
                            g_in_neighs[jj] = &gio[NTNow[16*ind+jj]];
                            if(jj<2) /* X direction */
			      h_out_neighs[jj] = &(hio[NTNow[16*ind+jj]]);
			    else
                              h_out_neighs[jj] = (Hermit*)&(htmp[NTNow[16*ind+jj]][(jj&1)?jj-3:jj-1][0][0][0]);
                        }

                        int ind_next = (t_next - t)*Pxyz+(z_next-z)*Pxy+(y_next-y)*Vxh+(x_next-x);
                        //const char accumulate = BTNow[Bbytes*ind];
                        const char isBoundary = BTNow[Bbytes*ind+2];
                        const long hprefdist = ind_next * num_floats_in_hermit_array;
                        const long gprefdist = ind_next * num_floats_in_gauge_array;

                        /* Momentum update body part dH from U */
			if((~BTNow[Bbytes*ind])==0) //if( mask != 0 )
			  //if(tid==0&&x+y+z+t==0)
                          gforce_imp_body_vec_noinline(
                                epsilonH,
                                //&hio[ind],
                                h_out_neighs,
                                g_in_neighs,
                                //accumulate,
                                //mask,
				isBoundary,
                                g_in_neighs[1],
                                g_in_neighs[3],
                                g_in_neighs[5],
                                g_in_neighs[7],
				kappaS, kappaR, kappaB,
                                gprefdist, gprefdist, gprefdist, gprefdist
                        );

                        _mm_prefetch( (char *)&NTNow[16*(ind+ind_next)], _MM_HINT_T0 );
                        _mm_prefetch( BTNow+Bbytes*(ind_next+ind), _MM_HINT_T0 );

                        x = x_next;
                        y = y_next;
                        z = z_next;
                        t = t_next;

                }

		/* Update momentum in face and send dH */
/*
		if(tid==0) printf("receiving u and sending h...\n");
		fflush(stdout);
*/
		gf_recv_and_unpack_u_and_send_boundaries_h(tid, hio, htmp, gio, kappaS, kappaR, kappaB, epsilonH, cb);
/*
		if(tid==0) printf("done receiving u and sending h\n");	
		fflush(stdout);
*/
		loop_ret = phaser->start(ps, tid, x, y, z, t);
		const char *BTNow2 = BoundTable+(1-cb)*PadBound;

		while(loop_ret) {
			loop_ret = phaser->next(ps, x_next, y_next, z_next, t_next);
			int ind = t*Pxyz+z*Pxy+y*Vxh+x;
			int ind_next = (t_next - t)*Pxyz+(z_next-z)*Pxy+(y_next-y)*Vxh+(x_next-x);
			int gprefdist = ind_next*sizeof(Gauge);
			char acctmp = BTNow2[Bbytes*ind];
			if(~(acctmp|3)!=0) acctmp = ((acctmp>>2)<<2); //acctmp &= (~3);
                        const char accumulate = acctmp;
			/* Update momentum and gauge fields body part */
			//if( (~BTNow[Bbytes*ind]) == 0 )
			if(accumulate != 0)
//#pragma omp critical
			  //if(tid==0&&x+y+z+t==0)
			  update_mg_body_vec_noinline(
                                epsilonU,
                                &gio[ind],
                                &hio[ind],
				&htmp[ind],
				accumulate, 
				gprefdist, gprefdist, gprefdist, gprefdist
                          );

			_mm_prefetch( BTNow2+Bbytes*(ind_next+ind), _MM_HINT_T0 );

                        x = x_next;
                        y = y_next;
                        z = z_next;
                        t = t_next;

		}
		/* Now update momentum and gauge fields requiring off-node dH */
/*
		if(tid==0) printf("receiving h...\n");
		fflush(stdout);
*/
		gf_recv_and_unpack_boundaries_h(tid, hio, htmp, gio, cb);
/*
		if(tid==0) printf("done receiving h\n");
		fflush(stdout);
*/
    } /* OMP loop end */
#if 0
    printf("calling gf_destroy_comms...\n");
    fflush(stdout);
    gf_destroy_comms_D();
    printf("done gf_destroy_comms\n");
    fflush(stdout);
#endif
}

