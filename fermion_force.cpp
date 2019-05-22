#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "qphix_internal.h"
#include "fermion_force.h"
#include "ks_globals.h"
#include "layout.h"
#include "qphix_su3_algebra.h"

#define set_temps() if(!setcount++) set_temps0()
#define free_temps() if(!--setcount) free_temps0()


#define NMTMP 4
#define NFTMP 3
#define NBTMP 3

extern int qphix_sites_on_node;
extern int qphix_even_sites_on_node;
extern site *lattice;
extern int qphix_fused_sites_on_node;
void dumpVal(QPHIX_ColorMatrix *p1, QPHIX_ColorMatrix *p2, QPHIX_ColorMatrix *p3, QPHIX_ColorMatrix *p4, int nu);
void dumpVal_1(QPHIX_ColorMatrix *p1, int nu, char *c1);
void dumpVal_1_s(QPHIX_ColorMatrix *p1, int nu, char *c1);
void dumpVal_V(QPHIX_ColorVector *p1, int nu, char *c1);

static QPHIX_ColorMatrix *mtmp[NMTMP], *ftmp0[NFTMP], *ftmp[NFTMP][4];
static QPHIX_ColorMatrix *btmp0[NBTMP], *btmp[NBTMP][4];
static int setcount=0;

#ifdef FF_PROFILE
static double shift_total = 0.0;
#endif

#if QPHIX_Colors == 'N'
static int gnc;
#define NC gnc
#define SETNC(x) gnc = x
#else
#define SETNC(x) (void)0
#endif

#define SETNCF(x) SETNC(QPHIX_get_nc(x))

QPHIX_hisq_links_t QPHIX_hisq_links =
{
  .inited = 0,
  .want_deps = 0,
  .want_aux=1,
  .reunit_allow_svd = 1,
  .reunit_svd_only = 0,
  .reunit_svd_rel_error = 1e-8,
  .reunit_svd_abs_error = 1e-8,
  .svd_values_info = 0,
  .use_fat7_lepage = 0,
};

QPHIX_hisq_force_t QPHIX_hisq_ff =
{
  .inited = 0,
  .fnmat_src_min = 1,
  .veclength = 4,
  .force_filter = 5e-5,
};

static double QPHIX_time(void)
{
  struct timeval tp;
  gettimeofday(&tp,NULL);
  return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6);
}

static void set_temps0(void)
{
  for(int i=0; i<NMTMP; i++) QPHIX_create_M(&mtmp[i], QPHIX_EVENODD);

  for(int i=0; i<NFTMP; i++)
  {
    QPHIX_create_M(&ftmp0[i], QPHIX_EVENODD);
    for(int j=0; j<4; j++) QPHIX_create_M(&ftmp[i][j], QPHIX_EVENODD);
  }

  for(int i=0; i<NBTMP; i++)
  {
    QPHIX_create_M(&btmp0[i], QPHIX_EVENODD);
    for(int j=0; j<4; j++) QPHIX_create_M(&btmp[i][j], QPHIX_EVENODD);
  }
}

static void free_temps0(void)
{
  for(int i=0; i<NMTMP; i++) QPHIX_destroy_M(mtmp[i]);

  for(int i=0; i<NFTMP; i++)
  {
    QPHIX_destroy_M(ftmp0[i]);
    for(int j=0; j<4; j++) QPHIX_destroy_M(ftmp[i][j]);
  }

  for(int i=0; i<NBTMP; i++)
  {
    QPHIX_destroy_M(btmp0[i]);
    for(int j=0; j<4; j++) QPHIX_destroy_M(btmp[i][j]);
  }
}

void QPHIX_hisq_force_multi(QPHIX_info_t *info, QPHIX_FermionLinksHisq *flh,
                            QPHIX_ColorMatrix *force[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            QPHIX_Real *residues, QPHIX_ColorVector *x[],
                            int *n_orders_naik)
{
  QPHIX_ColorMatrix *deriv[4];
  QPHIX_ColorMatrix *mtmp;
  double dtime = QPHIX_time();
  
  for(int mu=0; mu<4; mu++)
  {
    QPHIX_create_M(&deriv[mu], QPHIX_EVENODD);
    QPHIX_M_eq_zero(deriv[mu]);
  }
  
  QPHIX_hisq_deriv_multi(info, flh, deriv, hisq_coeff, residues, x, n_orders_naik);
  
  QPHIX_create_M(&mtmp, QPHIX_EVENODD);
  
  for(int dir=0; dir<4; dir++)
  {
    QPHIX_M_eq_M_times_Ma(mtmp, flh->U_links[dir], deriv[dir]);
    QPHIX_M_eq_antiherm_M(deriv[dir], mtmp);
    QPHIX_M_peq_M(force[dir], deriv[dir]); 
  }
  
  QPHIX_destroy_M(mtmp);
  for(int mu=0; mu<4; mu++)
  {
    QPHIX_destroy_M(deriv[mu]);
  }
  
  info->final_flop += (4.*(198+24+18))*qphix_sites_on_node;
  info->final_sec = QPHIX_time() - dtime;
}

void QPHIX_hisq_deriv_multi(QPHIX_info_t *info, QPHIX_FermionLinksHisq *flh,
                            QPHIX_ColorMatrix *deriv[],
                            QPHIX_hisq_coeffs_t *hisq_coeff,
                            QPHIX_Real *residues, QPHIX_ColorVector *x[],
                            int *n_orders_naik)
{
  struct timeval qphix_start, qphix_end;
  double qphix_total;
  double dtime = QPHIX_time();
  double totalflops = 0;
  int siteflops = 0;
  QPHIX_info_t tinfo;
  
  QPHIX_ColorMatrix *Ugf[4], *Vgf[4], *Wgf[4];
  QPHIX_ColorMatrix *force_accum_0[4];
  QPHIX_ColorMatrix *force_accum_0_naik[4];
  QPHIX_ColorMatrix *force_accum_1[4];
  QPHIX_ColorMatrix *force_accum_1u[4];
  QPHIX_ColorMatrix *force_accum_2[4];
  QPHIX_ColorMatrix *force_final[4];
  QPHIX_ColorMatrix *tmat;
  
  int n_naiks, nterms, n_naik_shift;
  
  for(int i=0; i<4; i++)
  {
    Ugf[i] = flh->U_links[i];
    Vgf[i] = flh->V_links[i];
    Wgf[i] = flh->W_unitlinks[i];
  }
  
  QPHIX_create_M(&tmat, QPHIX_EVENODD);
  for(int i=0; i<4; i++)
  {
    QPHIX_create_M(&force_accum_0[i], QPHIX_EVENODD);
    QPHIX_create_M(&force_accum_0_naik[i], QPHIX_EVENODD);
    QPHIX_create_M(&force_accum_1[i], QPHIX_EVENODD);
    QPHIX_create_M(&force_accum_1u[i], QPHIX_EVENODD);
    QPHIX_create_M(&force_accum_2[i], QPHIX_EVENODD);
    QPHIX_create_M(&force_final[i], QPHIX_EVENODD);
    QPHIX_M_eq_zero(force_accum_2[i]);
  }
  
  nterms = 0;
  n_naik_shift = 0;
  n_naiks = hisq_coeff->n_naiks;
  
  for(int inaik = 0; inaik < n_naiks; inaik++)
  {
    nterms += n_orders_naik[inaik];
  }
  
  for(int inaik=0; inaik<n_naiks; inaik++)
  {
    int n_orders_naik_current;
    if( inaik==0 )
    {
      n_orders_naik_current = nterms;
    }
    else
    {
      n_orders_naik_current = n_orders_naik[inaik];
    }

    /* get_mid */
#ifdef FF_PROFILE
    gettimeofday(&qphix_start, NULL);
#endif    
    QPHIX_get_mid(&tinfo, force_accum_0, x+n_naik_shift, residues+n_naik_shift,
                  1, n_orders_naik_current, 1);
#ifdef FF_PROFILE
    gettimeofday(&qphix_end, NULL);
    qphix_total = (qphix_end.tv_sec + (double) qphix_end.tv_usec/1000000) -
                  (qphix_start.tv_sec + (double) qphix_start.tv_usec/1000000);
    printf("qphix_get_mid time = %f\n", qphix_total);
#endif
    totalflops += tinfo.final_flop;
#ifdef FF_DEBUG
    //QPHIX_checksum_M(force_accum_0[0], 0);
#endif
    QPHIX_get_mid(&tinfo, force_accum_0_naik, x+n_naik_shift, residues+n_naik_shift,
                  1, n_orders_naik_current, 3);
    totalflops += tinfo.final_flop;
    
    for(int dir=0; dir<4; dir++)
    {
      QPHIX_M_eqm_M(force_accum_0[dir], force_accum_0[dir], QPHIX_ODD);
      QPHIX_M_eqm_M(force_accum_0_naik[dir], force_accum_0_naik[dir], QPHIX_ODD);
    }
    
    // Smearing 
    for(int i=0; i<4; i++) QPHIX_M_eq_zero(force_accum_1[i]);
    if(inaik==0) {
      QPHIX_asqtad_coeffs_t acoef;
      acoef.one_link = hisq_coeff->asqtad_one_link;
      acoef.three_staple = hisq_coeff->asqtad_three_staple;
      acoef.five_staple = hisq_coeff->asqtad_five_staple;
      acoef.seven_staple = hisq_coeff->asqtad_seven_staple;
      acoef.lepage = hisq_coeff->asqtad_lepage;
      acoef.naik = hisq_coeff->asqtad_naik;
      QPHIX_asqtad_deriv(&tinfo, Wgf, force_accum_1, &acoef,
          force_accum_0, force_accum_0_naik);
      totalflops += tinfo.final_flop;
    } else {
      QPHIX_asqtad_coeffs_t acoef;
      acoef.one_link = hisq_coeff->difference_one_link;
      acoef.three_staple = 0;
      acoef.five_staple = 0;
      acoef.seven_staple = 0;
      acoef.lepage = 0;
      acoef.naik = hisq_coeff->difference_naik;
      QPHIX_asqtad_deriv(&tinfo, Wgf, force_accum_1, &acoef,
          force_accum_0, force_accum_0_naik);
      totalflops += tinfo.final_flop;
    }

    QPHIX_Real coeff_mult;
    if( inaik==0 ) {
      coeff_mult = 1.0;
    } else {
      coeff_mult = hisq_coeff->eps_naik[inaik];
    }
    for(int dir=0; dir<4; dir++) {
      QPHIX_M_peq_r_times_M(force_accum_2[dir], coeff_mult, force_accum_1[dir]);
    }
    siteflops += 4*36;

    n_naik_shift += n_orders_naik[inaik];
  }

  QPHIX_asqtad_coeffs_t acoef;
  acoef.one_link = hisq_coeff->fat7_one_link;
  acoef.three_staple = hisq_coeff->fat7_three_staple;
  acoef.five_staple = hisq_coeff->fat7_five_staple;
  acoef.seven_staple = hisq_coeff->fat7_seven_staple;
  acoef.lepage = 0;
  acoef.naik = 0;
  if(QPHIX_hisq_links.use_fat7_lepage) {
    acoef.lepage = hisq_coeff->fat7_lepage;
  }

  /* QPHIX_hisq_unitarize_method_t umethod = hisq_coeff->umethod; */
  QPHIX_hisq_unitarize_method_t umethod = QPHIX_UNITARIZE_RATIONAL; 
  if ( umethod==QPHIX_UNITARIZE_NONE ){

    for(int dir=0; dir<4; dir++)
      QPHIX_M_eq_zero(force_accum_1[dir]);
#ifdef FF_DEBUG      
    printf("umethod==QPHIX_UNITARIZE_NONE\n");
#endif      
    QPHIX_asqtad_deriv(&tinfo, Ugf, force_accum_1, &acoef, force_accum_2, NULL);
    totalflops += tinfo.final_flop;

  } else if ( umethod==QPHIX_UNITARIZE_RATIONAL ) {

    for(int mu=0; mu<4; mu++) QPHIX_M_eq_Ma(force_accum_1u[mu], force_accum_2[mu]);

    QPHIX_hisq_force_multi_reunit(&tinfo, Vgf, force_accum_2, force_accum_1u);
    for(int mu=0; mu<4; mu++) QPHIX_M_eq_Ma(force_accum_1u[mu], force_accum_2[mu]);
    totalflops += tinfo.final_flop;

    for(int dir=0; dir<4; dir++) QPHIX_M_eq_zero(force_accum_1[dir]);
    QPHIX_asqtad_deriv(&tinfo, Ugf, force_accum_1, &acoef,
        force_accum_1u, NULL);
totalflops += tinfo.final_flop;

} else {
  printf("Unknown or unsupported unitarization method\n");
  exit(1);
}

for(int dir=0; dir<4; dir++) {
  QPHIX_Real treal = 2;
  QPHIX_M_peq_r_times_M(deriv[dir], treal, force_accum_1[dir]);
}
siteflops += 4*36;

for(int i=0; i<4; i++) {
  QPHIX_destroy_M( force_accum_0[i] );
  QPHIX_destroy_M( force_accum_0_naik[i] );
  QPHIX_destroy_M( force_accum_1[i] );
  QPHIX_destroy_M( force_accum_1u[i] );
  QPHIX_destroy_M( force_accum_2[i] );
  QPHIX_destroy_M( force_final[i] );
}
QPHIX_destroy_M( tmat );

totalflops += ((double)siteflops)*qphix_sites_on_node;
info->final_sec = QPHIX_time() - dtime;
info->final_flop = totalflops;
info->status = QPHIX_SUCCESS;
}


  void QPHIX_get_mid(QPHIX_info_t *info, QPHIX_ColorMatrix **mid,
      QPHIX_ColorVector *multix[], QPHIX_Real eps[],
      QPHIX_Real scale, int nterms, int ndist)
  {
#ifdef FF_PROFILE
    struct timeval sv_start, sv_end;
    double sv_total = 0.0;
#endif
    register int dim, dir, term, s, i, j, v, k;
    register QPHIX_Real r;
    QSU3V *v1e, *v1o, *v2e, *v2o;
    QSU3M *m1e, *m1o;
    QPHIX_ColorVector *shiftx[2];
    QPHIX_ColorMatrix *temp;

    QPHIX_create_V(&shiftx[0], QPHIX_EVENODD);
    QPHIX_create_V(&shiftx[1], QPHIX_EVENODD);
    QPHIX_create_M(&temp, QPHIX_EVENODD);

    //Algo using only primitives

    for(dim=XUP; dim<=TUP; dim++)
    {
      dir = (ndist == 3) ? dim+8 : dim;
      QPHIX_M_eq_zero(mid[dim]);

      for(int i=-1; i<nterms; i++) {
       if(i+1<nterms) {
        int k = (i+1)%2;
        QPHIX_V_eq_sV(shiftx[k], multix[i+1], dir, FORWARDS);
       }
       if(i>=0) {
        int k = i%2;
        r = eps[i] * scale;
        QPHIX_M_eq_V_times_Va(temp, multix[i], shiftx[k]);
        QPHIX_M_peq_r_times_M(mid[dim], r, temp);
       }
      }
    }

QPHIX_destroy_V(shiftx[0]);
QPHIX_destroy_V(shiftx[1]);
QPHIX_destroy_M(temp);

info->status = QPHIX_SUCCESS;
#ifdef FF_PROFILE
printf("shift inside get_mid = %f\n", sv_total);
#endif
}

/* Smearing */
  static void
staple(QPHIX_ColorMatrix *out, QPHIX_ColorMatrix *in0, QPHIX_ColorMatrix *link0, int mu, int nu)
{
#define link     ftmp0[0]
#define linkmu   ftmp[0][mu]
#define in       ftmp0[1]
#define innu     ftmp[1][nu]
#define linkin   mtmp[0]
#define back     btmp0[0]
#define backnu   btmp[0][nu]
#define linkinnu mtmp[1]
  struct timeval staple_end, staple_start;
  QPHIX_M_eq_sM(linkmu, link0, mu, FORWARDS);
  QPHIX_M_eq_sM(innu, in0, nu, FORWARDS);
  QPHIX_M_eq_Ma_times_M(linkin, link0, in0);
  QPHIX_M_eq_M_times_M(back, linkin, linkmu);
  QPHIX_M_eq_sM(backnu, back, nu, BACKWARDS); 
  QPHIX_M_eq_M_times_M(linkinnu, link0, innu);
  QPHIX_M_peq_M_times_Ma(out, linkinnu, linkmu);
  QPHIX_M_peq_M(out, backnu);
#define STAPLE_FLOPS (3*198+216+18)

#undef link
#undef linkmu
#undef in
#undef innu
#undef linkin
#undef back
#undef backnu
#undef linkinnu
}

  static void
side_force(QPHIX_ColorMatrix *force, QPHIX_ColorMatrix *bot0, QPHIX_ColorMatrix *side0,
    QPHIX_ColorMatrix *top0, int mu, int nu, QPHIX_ColorMatrix *stpl)
{
#define side     ftmp0[0]
#define sidemu   ftmp[0][mu]
#define top      ftmp0[1]
#define topnu    ftmp[1][nu]
#define bot      ftmp0[2]
#define botnu    ftmp[2][nu]
#define sidebot  mtmp[0]
#define sidetop  mtmp[1]
#define topnusidebot  btmp0[0]
#define fbmu          btmp[0][mu]
#define botnusidetop  btmp0[1]
#define fmbmu         btmp[1][mu]
#define sidebotsidemu btmp0[2]
#define stm           btmp[2][nu]
#define botnusidemu   mtmp[2]
#define botsidemu     mtmp[3]

  // force += bot * sidemu * topnu+
  // force -= bot-mu+ * side-mu * topnu-mu
  // -= top <-> bot
  // stpl += side * botnu * sidemu+
  // stpl += side-nu+ * bot-nu * sidemu-nu
  QPHIX_M_eq_sM(sidemu, side0, mu, FORWARDS); 
  QPHIX_M_eq_sM(topnu, top0, nu, FORWARDS); 
  QPHIX_M_eq_sM(botnu, bot0, nu, FORWARDS); 
  QPHIX_M_eq_Ma_times_M(sidebot, side0, bot0);
  QPHIX_M_eq_Ma_times_M(sidetop, side0, top0);
  QPHIX_M_eq_Ma_times_M(topnusidebot, topnu, sidebot);
  QPHIX_M_eq_sM(fbmu, topnusidebot, mu, BACKWARDS); 
  QPHIX_M_eq_Ma_times_M(botnusidetop, botnu, sidetop);
  QPHIX_M_eq_sM(fmbmu, botnusidetop, mu, BACKWARDS); 
  QPHIX_M_eq_M_times_M(sidebotsidemu, sidebot, sidemu);
  QPHIX_M_eq_sM(stm, sidebotsidemu, nu, BACKWARDS);
  QPHIX_M_eq_M_times_Ma(botnusidemu, botnu, sidemu);
  QPHIX_M_peq_M_times_M(stpl, side0, botnusidemu);
  QPHIX_M_peq_M_times_Ma(force, top0, botnusidemu);
  QPHIX_M_eq_M_times_M(botsidemu, bot0, sidemu);
  QPHIX_M_peq_M_times_Ma(force, botsidemu, topnu);
  QPHIX_M_peq_Ma(force, fbmu);
  QPHIX_M_peq_Ma(force, fmbmu);
  QPHIX_M_peq_M(stpl, stm);
#define SIDE_FORCE_FLOPS (7*198+3*216+3*18)

#undef side
#undef sidemu
#undef top
#undef topnu
#undef bot
#undef botnu
#undef sidebot
#undef sidetop
#undef topnusidebot
#undef fbmu
#undef botnusidetop
#undef fmbmu
#undef sidebotsidemu
#undef stm
#undef botnusidemu
#undef botsidemu
}

  static void
QPHIX_fat_deriv(QPHIX_info_t *info, QPHIX_ColorMatrix *gauge[],
    QPHIX_ColorMatrix *deriv[], QPHIX_asqtad_coeffs_t *coef,
    QPHIX_ColorMatrix *mid[])
{
  double dtime = QPHIX_time();
  double nflops = 0;
#ifdef FF_PROFILE 
  struct timeval staple_start, staple_end, sideforce_start, sideforce_end, fatshift_start, fatshift_end;
  static double staple_total  = 0.0, sideforce_total=0.0, fatshift_total=0.0;
#endif
  QPHIX_Real coef1 = coef->one_link;
  QPHIX_Real coef3 = coef->three_staple;
  QPHIX_Real coef5 = coef->five_staple;
  QPHIX_Real coef7 = coef->seven_staple;
  QPHIX_Real coefL = coef->lepage;
  coef1 -= 6*coefL;
  int have5 = coef5 || coef7 || coefL;
  int have3 = coef3 || have5;
  QPHIX_ColorMatrix *stpl3[4];
  QPHIX_ColorMatrix *mid5[4];
  QPHIX_ColorMatrix *stpl5=NULL, *mid3=NULL;
  if(have3) {
    if(have5) {
      for(int mu=0; mu<4; mu++) {
        QPHIX_create_M(&stpl3[mu], QPHIX_EVENODD);
        QPHIX_create_M(&mid5[mu], QPHIX_EVENODD);
      }
    }
    QPHIX_create_M(&stpl5, QPHIX_EVENODD);
    QPHIX_create_M(&mid3, QPHIX_EVENODD);
    set_temps();
  }

  for(int sig=0; sig<4; sig++) {
    if(have5) {
      for(int mu=0; mu<4; mu++) {
        if(mu==sig) continue;
        QPHIX_M_eq_zero(stpl3[mu]);
        staple(stpl3[mu], gauge[sig], gauge[mu], sig, mu);
        QPHIX_M_eq_zero(mid5[mu]);
        nflops += STAPLE_FLOPS;
      } // end of mu loop
      for(int rho=0; rho<4; rho++) {
        if(rho==sig) continue;
        QPHIX_M_eq_zero(stpl5);

        if(coef7) {
          for(int mu=0; mu<4; mu++) {
            if(mu==sig||mu==rho) continue;
            for(int nu=0; nu<4; nu++) {
              if(nu==sig||nu==rho||nu==mu) continue;
              staple(stpl5, stpl3[mu], gauge[nu], sig, nu);
              nflops += STAPLE_FLOPS;
            }
          }
          QPHIX_M_eq_r_times_M(stpl5, coef7, stpl5);
          nflops += 18;
        } // end of coef7 loop

        QPHIX_M_eq_zero(mid3);
        if(coefL) {
          QPHIX_M_peq_r_times_M(stpl5, coefL, stpl3[rho]);
          nflops += 36;
        }
        if(coefL || coef7) {
          side_force(deriv[rho], mid[sig], gauge[rho], stpl5, sig, rho, mid3);
          nflops += SIDE_FORCE_FLOPS;
        }
        if(coefL) {
          QPHIX_M_peq_r_times_M(mid5[rho], coefL, mid3);
          nflops += 36;
        }

        QPHIX_M_eq_r_times_M(mid3, coef7, mid3);
        QPHIX_M_peq_r_times_M(mid3, coef5, mid[sig]);
        nflops += 18+36;
        for(int mu=0; mu<4; mu++) {
          if(mu==sig||mu==rho) continue;
          for(int nu=0; nu<4; nu++) {
            if(nu==sig||nu==rho||nu==mu) continue;
            side_force(deriv[mu], mid3, gauge[mu], stpl3[nu], sig, mu, mid5[nu]);
            nflops += SIDE_FORCE_FLOPS;
          } // end of nu loop
        } // end of mu loop
      } // end of rho loop
    } // end of have5 condition

    if(have3) {
      QPHIX_M_eq_zero(mid3);
      for(int mu=0; mu<4; mu++) {
        if(mu==sig) continue;
        if(have5) {
          QPHIX_M_eq_r_times_M_plus_M(stpl5, coef3, mid[sig], mid5[mu]);
          nflops += 36;
        } else {
          QPHIX_M_eq_r_times_M(stpl5, coef3, mid[sig]);
          nflops += 18;
        }
        side_force(deriv[mu], stpl5, gauge[mu], gauge[sig], sig, mu, mid3);
        nflops += SIDE_FORCE_FLOPS;
      } // end of mu loop
      QPHIX_M_peq_M(deriv[sig], mid3);
      nflops += 18;
    } // end of have3 condition
    if(coef1) {
      QPHIX_M_peq_r_times_M(deriv[sig], coef1, mid[sig]);
      nflops += 36;
    }
  } // end of sig loop

#if 1
  if(coefL) {
    // fix up Lepage term
    QPHIX_Real fixL = -coefL;
    for(int mu=0; mu<4; mu++) {
      QPHIX_M_eq_zero(ftmp0[0]);
      for(int nu=0; nu<4; nu++) {
        if(nu==mu) continue;
        QPHIX_M_eq_Ma_times_M(btmp0[0], mid[nu], gauge[nu]);
        QPHIX_M_eq_sM(btmp[0][nu], btmp0[0], nu, BACKWARDS);
        QPHIX_M_eq_M_times_Ma(stpl5, mid[nu], gauge[nu]);
        QPHIX_M_meq_M(stpl5, btmp[0][nu]);
        QPHIX_M_peq_M(ftmp0[0], stpl5);
        QPHIX_M_peq_Ma(ftmp0[0], stpl5);
      }  // end of nu loop
      QPHIX_M_eq_sM(ftmp[0][mu], ftmp0[0], mu, FORWARDS); 
      QPHIX_M_eq_M_times_M(stpl5, ftmp0[0], gauge[mu]);
      QPHIX_M_meq_M_times_M(stpl5, gauge[mu], ftmp[0][mu]);
      QPHIX_M_peq_r_times_M(deriv[mu], fixL, stpl5);
    } // end of mu loop
    nflops += 4*(3*(2*198+3*18)+198+216+36);
  } // end of coefL condition
#endif

  if(have3) {
    if(have5) {
      for(int mu=0; mu<4; mu++) {
        QPHIX_destroy_M(stpl3[mu]);
        QPHIX_destroy_M(mid5[mu]);
      }
    } // end of mu loop
    QPHIX_destroy_M(stpl5);
    QPHIX_destroy_M(mid3);
    free_temps();
  }

  info->final_sec = QPHIX_time() - dtime;
  info->final_flop = nflops*qphix_sites_on_node;
  info->status = QPHIX_SUCCESS;
#ifdef FF_PROFILE
  printf("staple time = %f\n", staple_total);
  printf("side force time = %f\n", sideforce_total);
#endif
}

void QPHIX_asqtad_deriv(QPHIX_info_t *info, QPHIX_ColorMatrix *gauge[],
    QPHIX_ColorMatrix *deriv[], QPHIX_asqtad_coeffs_t *coef,
    QPHIX_ColorMatrix *mid_fat[],
    QPHIX_ColorMatrix *mid_naik[])
{
  SETNCF(deriv[0]);
  double dtime = QPHIX_time();
  double nflops = 0;
  QPHIX_info_t tinfo;
  set_temps();
#ifdef FF_PROFILE
  struct timeval qphix_start, qphix_end;
  double qphix_total;
#endif
  // fat
  if(mid_fat) {
    QPHIX_fat_deriv(&tinfo, gauge, deriv, coef, mid_fat);
    nflops += tinfo.final_flop;
  }

  // Naik
  QPHIX_Real coefN = coef->naik;
  if(coefN && mid_naik) {
#define U       gauge[mu]
#define mid     mid_naik[mu]
#define Uf      ftmp0[0]
#define Umu     ftmp[0][mu]
#define Umid    btmp0[0]
#define Umidbmu btmp[0][mu]
#define UmuU    ftmp0[1]
#define UmuUs   ftmp[1][mu]
#define f3b     btmp0[1]
#define f3      btmp[1][mu]
#define f       mtmp[0]
    for(int mu=0; mu<4; mu++) {
      QPHIX_M_eq_sM(Umu, U, mu, FORWARDS); 
      QPHIX_M_eq_Ma_times_M(Umid, U, mid);
      QPHIX_M_eq_sM(Umidbmu, Umid, mu, BACKWARDS); 
      QPHIX_M_eq_Ma_times_Ma(UmuU, Umu, U);
      QPHIX_M_eq_sM(UmuUs, UmuU, mu, FORWARDS); 
      QPHIX_M_eq_Ma_times_M(f3b, U, Umidbmu);
      QPHIX_M_eq_sM(f3, f3b, mu, BACKWARDS); 
      QPHIX_M_eq_M_times_M(f, mid, UmuUs);
      QPHIX_M_peq_M_times_Ma(f, Umidbmu, Umu);
      QPHIX_M_peq_M(f, f3);
      QPHIX_M_peq_r_times_M(deriv[mu], coefN, f);
    }
#undef U
#undef mid
#undef Uf
#undef Umu
#undef Umid
#undef Umidbmu
#undef UmuU
#undef UmuUs
#undef f3b
#undef f3
#undef f
    nflops += 4*(4*198+216+18+36)*qphix_sites_on_node;
  }
  free_temps();

  info->final_sec = QPHIX_time() - dtime;
  info->final_flop = nflops;
  info->status = QPHIX_SUCCESS;
}

  void
QPHIX_hisq_force_multi_reunit(QPHIX_info_t *info,
    QPHIX_ColorMatrix *V[4],
    QPHIX_ColorMatrix *Force[4],
    QPHIX_ColorMatrix *Force_old[4])
{

  QPHIX_ColorTensor4 *dwdv = (QPHIX_ColorTensor4*)malloc(sizeof(QPHIX_ColorTensor4)*qphix_sites_on_node);
  QPHIX_ColorTensor4 *dwdagdv = (QPHIX_ColorTensor4*)malloc(sizeof(QPHIX_ColorTensor4)*qphix_sites_on_node);
  for(int i=0; i<4; i++ ) {
    //expose QPhiX links
    qphix_su3_matrix *Vlinks = QPHIX_expose_M( V[i] );
    qphix_su3_matrix *ff_new = QPHIX_expose_M( Force[i] );
    qphix_su3_matrix *ff_old = QPHIX_expose_M( Force_old[i] );

    //AB NO NEED TO REPHASE V LINKS???
#pragma omp parallel
    {
#pragma omp for
      for(int j=0; j<qphix_sites_on_node; j++ ) {
  //      QPHIX_ColorTensor4 dwdv, dwdagdv;
        // derivative with respect to V and V^+
        u3_un_der_analytic(info, &(Vlinks[j]), &dwdv[j], &dwdagdv[j]);
        // adjoint piece of force from the previous level
        qphix_su3_matrix ff_old_adj;
        QSU3_M_eq_Ma( &ff_old_adj, &( ff_old[j] ) ); 
        // see LONG COMMENT in fermion_force_fn_multi_hisq
        QSU3_M_eq_zero( &( ff_new[j] ) ); 
        for(int m=0; m<3; m++) {
          for(int n=0; n<3; n++) {
            for(int k=0; k<3; k++) {
              for(int l=0; l<3; l++) {
                // direct part
                QPHIX_Complex ftmp;
                QPHIX_c_eq_c_times_c( ftmp, dwdv[j].t4[k][m][n][l], QPHIX_elem_M(ff_old[j],l,k) );
                QPHIX_c_peq_c( QPHIX_elem_M(ff_new[j],n,m), ftmp );
                // adjoint part
                QPHIX_c_eq_c_times_c( ftmp, dwdagdv[j].t4[k][m][n][l], QPHIX_elem_M(ff_old_adj,l,k) );
                QPHIX_c_peq_c( QPHIX_elem_M(ff_new[j],n,m), ftmp );
              }
            }
          }
        }
      }
    }
    // resume QDP operations on links
    QPHIX_reset_M(Force[i], ff_new );
  }
}

// Determinant of 3x3 complex matrix
QPHIX_Complex su3_mat_det(qphix_su3_matrix *U)
{
  QPHIX_Complex a, b, m0, m1, m2, cdet;
  QPHIX_c_eq_c_times_c( a, Uelem(1,1), Uelem(2,2) );
  QPHIX_c_eq_c_times_c( b, Uelem(1,2), Uelem(2,1) );
  QPHIX_c_eq_c_minus_c( m0, a, b );

  QPHIX_c_eq_c_times_c( a, Uelem(1,0), Uelem(2,2) );
  QPHIX_c_eq_c_times_c( b, Uelem(1,2), Uelem(2,0) );
  QPHIX_c_eq_c_minus_c( m1, a, b );

  QPHIX_c_eq_c_times_c( a, Uelem(1,0), Uelem(2,1) );
  QPHIX_c_eq_c_times_c( b, Uelem(1,1), Uelem(2,0) );
  QPHIX_c_eq_c_minus_c( m2, a, b );

  QPHIX_c_eq_c_times_c( a, Uelem(0,0), m0 );
  QPHIX_c_eq_c_times_c( b, Uelem(0,1), m1 );
  QPHIX_c_eq_c_minus_c( cdet, a, b );
  QPHIX_c_eq_c_times_c( a, Uelem(0,2), m2 );
  QPHIX_c_eq_c_plus_c( cdet, cdet, a );

  return cdet;
}

/* Analytic derivative of the unitarized matrix with respect to the original:
   dW/dV and d(W^+)/dV, where W=V(V^+V)^-1/2
   */
void u3_un_der_analytic(QPHIX_info_t *info, qphix_su3_matrix *V,
    QPHIX_ColorTensor4 *dwdv, QPHIX_ColorTensor4 *dwdagdv)
{
  int i, j, m, n, perform_svd, this_node=0;
  QPHIX_Complex det, der, ctmp, ctmp2;
  QPHIX_Real det_check=0, c0, c1, c2, S, g0, g1, g2, R, S3, RoS, theta, theta3, pi23;
  QPHIX_Real g0sq, g1sq, g2sq, us, vs, ws, f0, f1, f2, denom;
  QPHIX_Real u2, u3, u4, u5, u6, u7, u8, v2, v3, v4, v5, v6, w2, w3, w4, w5;
  QPHIX_Real b00, b01, b02, b11, b12, b22, denom3;
  qphix_su3_matrix Q12, Q, Q2, Q3, Uleft, Vright, W, S1, S2;
  qphix_su3_matrix VVd, VQ, QVd, QQVd, VQQ, VQVd, PVd, RVd, SVd, Vd;
  QPHIX_Real sigma[3];
  size_t nflops = 0;

  //this_node = layout.this_node;

  // No SVD initially
  perform_svd=0;

  if(QPHIX_hisq_links.reunit_allow_svd)
  {
    // get determinant for future comparison
    det = su3_mat_det(V);
    // det_check = |det|^2
    det_check = (det.real * det.real) + (det.imag * det.imag);
  }

  /* adjoint */
  QSU3_M_eq_Ma(&Vd, V);

  /* Hermitian matrix: Q=V^+ V */
  QSU3_M_eq_Ma_times_M(&Q, V, V);
  nflops += 198;

  /* Q^2 */
  QSU3_M_eq_M_times_M(&Q2, &Q, &Q);
  nflops += 198;

  /* Q^3 */
  QSU3_M_eq_M_times_M(&Q3, &Q2, &Q);
  nflops += 198;

  /* (real) traces */
  QSU3_R_eq_re_trace_M(&c0, &Q);
  QSU3_R_eq_re_trace_M(&c1, &Q2);
  c1 /= 2;
  QSU3_R_eq_re_trace_M(&c2, &Q3);
  c2 /= 3;

  S = c1/3 - c0 * (c0/18);
  nflops += 24;

  if (fabs(S) < QPHIX_U3_UNIT_ANALYTIC_EPS)
  {
    /* eigenvalues of Q */
    g0 = c0/3;
    g1 = c0/3;
    g2 = c0/3;
    nflops += 3;
  }
  else
  {
    R = c2/2 - c0 * (c1/3) + c0 * c0 * (c0/27);
    //R2 = R*R;
    //CQ3 = S*S*S;
    S = sqrt(S);
    S3 = S * S * S;
    /* treat possible underflow: R/S^3/2>1.0 leads to acos giving NaN */
    RoS = R/S3;

    nflops += 14;

    if (!(fabs(RoS) < 1.0))
    {
      if (R > 0)
      {
        theta = 0.0;
      }
      else
      {
        theta = QPHIX_PI;
      }
    } 
    else
    {
      theta = acos(RoS);
      if (isnan(theta))
      {
        printf("Hit NaN in QOPPC(u3_un_der_analytic)\n");
        printf("RoS=%24.18g\n",RoS);
        printf("Matrix V (row-wise):\n");
        for (i=0; i<3; i++)
        {
          for (j=0; j<3; j++)
          {
            printf("%24.18g %24.18g\n", QPHIX_real(QPHIX_elem_M(*V, i, j)), QPHIX_imag(QPHIX_elem_M(*V, i, j)));
          }
          nflops += 1;
        }
        abort();
      }
    }

    /* eigenvalues of Q */
    theta3 = theta / 3;
    pi23 = (2 * QPHIX_PI) / 3;
    g0 = c0 / 3 + 2 * S * cos(theta3);
    g1 = c0 / 3 + 2 * S * cos(theta3 + pi23);
    g2 = c0 / 3 + 2 * S * cos(theta3 + 2 * pi23);

    nflops += 20;
  }

  if (QPHIX_hisq_links.reunit_allow_svd)
  {
    if (!QPHIX_hisq_links.reunit_svd_only)
    {
      /* conditions to call SVD */
      if (det_check != 0)
      {
        if (fabs(det_check - g0 * g1 * g2) / fabs(det_check) > QPHIX_hisq_links.reunit_svd_rel_error)
        {
          perform_svd = 1;
          nflops += 4;
        }
      }
      if (det_check < QPHIX_hisq_links.reunit_svd_abs_error)
      {
        perform_svd = 1;
      }
    }
    else
    {
      /* exclusively use SVD for finding eigenvalues and reunitarization,
       * this is slow since Q, Q^2 and Q^3 are calculated anyway;
       * this option is for testing: under normal circumstances SVD is
       * rarely used, which makes it harder to test, therefore one can
       * SVD with this switch
       */
      perform_svd = 1;
    }
  } /* reunit_allow_svd */

  if (perform_svd != 0)
  {
    /* call SVD */
    QPHIX_svd3x3(info, V, sigma, &Uleft, &Vright);

    if (!QPHIX_hisq_links.reunit_svd_only && QPHIX_hisq_links.svd_values_info)
    {
      printf("*** QOPQDP reunitarization while calculating fat links ***\n");
      printf("*** Resort to using svd (force) ***\n");
      printf("*** printing from node %d ***\n",this_node);
      printf("Eigenvalues from cubic equation:\n");
      printf("  g0 = %28.18f\n", g0);
      printf("  g1 = %28.18f\n", g1);
      printf("  g2 = %28.18f\n", g2);
      printf("Eigenvalues from singular value decomposition:\n");
      printf("  g0 = %28.18f\n", sigma[0] * sigma[0]);
      printf("  g1 = %28.18f\n", sigma[1] * sigma[1]);
      printf("  g2 = %28.18f\n", sigma[2] * sigma[2]);
    }

    g0 = sigma[0] * sigma[0];
    g1 = sigma[1] * sigma[1];
    g2 = sigma[2] * sigma[2];

    nflops += 3;
  }

  if (QPHIX_hisq_ff.force_filter > 0)
  {
    QPHIX_Real gmin, g_epsilon;
    gmin = g0;
    if (g1 < gmin) gmin = g1;
    if (g2 < gmin) gmin = g2;
    g_epsilon = QPHIX_hisq_ff.force_filter;
    if (gmin < g_epsilon)
    {
      g0 += g_epsilon;
      g1 += g_epsilon;
      g2 += g_epsilon;
      nflops += 3;

      // modify also Q and Q2 matrices
      for (i=0; i<3; i++)
      {
        //Q.e[i][i].real+=g_epsilon;
        QPHIX_c_peq_r(QPHIX_elem_M(Q, i, i), g_epsilon);
        nflops += 1;
      }
      QSU3_M_eq_M_times_M(&Q2, &Q, &Q);
      nflops += 198;
      //QOP_info_hisq_force_filter_counter(info)++;
    }
  }

  /* roots of eigenvalues */
  g0sq = sqrt(g0);
  g1sq = sqrt(g1);
  g2sq = sqrt(g2);
  // +3 flops (sqrt counted as 1 flop)

  /* symmetric combinations */
  us = g1sq + g2sq;
  ws = g1sq * g2sq;
  vs = g0sq * us + ws;
  us += g0sq; 
  ws *= g0sq;
  // +6 flops

  if (ws < QPHIX_U3_UNIT_ANALYTIC_EPS)
  {
    printf( "WARNING: u3_un_der_analytic: ws is too small!\n" );
  }

  denom = ws * (us * vs - ws);
  // +3 flops

  /* constants in inverse root expression */
  f0 = (us * vs * vs - ws * (us * us + vs)) / denom;
  f1 = (2 * us * vs - ws - us * us * us) / denom;
  f2 = us / denom;
  // +15 flops

  nflops += 12+15;

  /* assemble inverse root: Q^-1/2 = f0 + f1*Q + f2*Q^2 */
  QSU3_M_eq_r_times_M(&S1, &f2, &Q2);
  QSU3_M_eq_r_times_M_plus_M(&Q12, &f1, &Q, &S1);
  QPHIX_c_peq_r(QPHIX_elem_M(Q12,0,0), f0);
  QPHIX_c_peq_r(QPHIX_elem_M(Q12,1,1), f0);
  QPHIX_c_peq_r(QPHIX_elem_M(Q12,2,2), f0);
  nflops += 18+36+3;

  /* W = V*S2 */
  QSU3_M_eq_M_times_M(&W, V, &Q12);

  denom3 = 2 * denom * denom * denom;
  // +3 flops
  nflops += 198 + 3;

  /* derivatives of coefficients: B_ij=df_i/dc_j */
  u2 = us * us;
  u3 = u2 * us;
  u4 = u3 * us;
  u5 = u4 * us;
  u6 = u5 * us;
  u7 = u6 * us;
  u8 = u7 * us;
  v2 = vs * vs;
  v3 = v2 * vs;
  v4 = v3 * vs;
  v5 = v4 * vs;
  v6 = v5 * vs;
  w2 = ws * ws;
  w3 = w2 * ws;
  w4 = w3 * ws;
  w5 = w4 * ws;
  // +16 flops
  b00 = - w3 * u6 + (3 * u4) * (vs * w3 + v4 * ws)
    - u3 * (v6 + 4 * w4 + 12 * v3 * w2) + u2 * (16 * v2 * w3 + 3 * v5 * ws)
    - us * ((8 * vs * w4) + (3 * v4 * w2)) + w5 + (v3 * w3);
  b00 /= denom3; // + 33 flops
  b01 = -w2*u7 -v2*ws*u6 +u5*( v4 +6*vs*w2 ) -u4*( 5*w3 +v3*ws )
    -u3*( 2*v5 +6*v2*w2 ) +u2*( 10*vs*w3 +6*v4*ws )
    -us*( 3*w4 +6*v3*w2 ) +2*v2*w3;
  b01 /= denom3; // +38 flops
  b02 = w2*u5 +v2*ws*u4 -u3*( v4 +4*vs*w2 )
    +u2*( 4*w3 +3*v3*ws ) -3*v2*w2*us +vs*w3;
  b02 /= denom3; // +22 flops
  b11 = -ws*u8 -v2*u7 +7*vs*ws*u6 +u5*( 4*v3 -5*w2 ) -16*v2*ws*u4
    -u3*( 4*v4 -16*vs*w2 ) -u2*( 3*w3 -12*v3*ws ) -12*v2*w2*us +3*vs*w3;
  b11 /= denom3; // +37 flops
  b12 = ws*u6 +v2*u5 -5*vs*ws*u4 -u3*( 2*v3 -4*w2 ) +6*v2*ws*u2 -6*vs*w2*us+w3;
  b12 /= denom3; // +22 flops
  b22 = -ws*u4 -v2*u3 +3*vs*ws*u2 -3*w2*us;
  b22 /= denom3; // +12 flops

  nflops += 16+33+38+22+37+22+12;

  /* ** create several building blocks for derivative ** */
  QSU3_M_eq_M_times_M(&VQ, V, &Q);
  QSU3_M_eq_M_times_Ma(&QVd, &Q, V);
  QSU3_M_eq_M_times_Ma(&VVd, V, V);
  QSU3_M_eq_M_times_M(&VQQ, V, &Q2);
  QSU3_M_eq_M_times_Ma(&QQVd, &Q2, V);
  QSU3_M_eq_M_times_Ma(&VQVd, &VQ, V);
  QSU3_M_eq_r_times_M(&S1, &b01, &QVd);
  QSU3_M_eq_r_times_M_plus_M(&S2, &b02, &QQVd, &S1);
  QSU3_M_eq_r_times_M_plus_M(&PVd, &b00, &Vd, &S2);
  QSU3_M_eq_r_times_M(&S1, &b11, &QVd);
  QSU3_M_eq_r_times_M_plus_M(&S2, &b12, &QQVd, &S1);
  QSU3_M_eq_r_times_M_plus_M(&RVd, &b01, &Vd, &S2);
  QSU3_M_eq_r_times_M(&S1, &b12, &QVd);
  QSU3_M_eq_r_times_M_plus_M(&S2, &b22, &QQVd, &S1);
  QSU3_M_eq_r_times_M_plus_M(&SVd, &b02, &Vd, &S2);

  nflops += 198*6+4*18+5*36;

  /* assemble the derivative rank 4 tensor */
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      for (m=0; m<3; m++)
      {
        for (n=0; n<3; n++)
        {
          QPHIX_c_eq_r(der, 0);
          /* dW/dV */
          if (i == m)
          {
            QPHIX_c_peq_c(der, QPHIX_elem_M(Q12, n, j));
            nflops += 2;
          }
          if (j == n)
          {
            QPHIX_c_peq_r_times_c(der, f1, QPHIX_elem_M(VVd, i, m));
            QPHIX_c_peq_r_times_c(der, f2, QPHIX_elem_M(VQVd, i, m));
            nflops += 4;
          }
          QPHIX_c_eq_c_times_c(ctmp, QPHIX_elem_M(*V, i, j), QPHIX_elem_M(PVd, n, m));
          QPHIX_c_peq_c(der, ctmp);
          QPHIX_c_eq_c_times_c(ctmp, QPHIX_elem_M(VQ, i, j), QPHIX_elem_M(RVd, n, m));
          QPHIX_c_peq_c(der, ctmp);
          QPHIX_c_eq_c_times_c(ctmp, QPHIX_elem_M(VQQ, i, j), QPHIX_elem_M(SVd, n, m));
          QPHIX_c_peq_c(der, ctmp);
          QPHIX_c_eq_c_times_c(ctmp, QPHIX_elem_M(VVd, i, m), QPHIX_elem_M(Q, n, j));
          QPHIX_c_peq_r_times_c(der, f2, ctmp);
          dwdv->t4[i][m][n][j].real = QPHIX_real(der);
          dwdv->t4[i][m][n][j].imag = QPHIX_imag(der);
          /* dW^+/dV */
          QPHIX_c_eq_c_times_c(der, QPHIX_elem_M(Vd, i, j), QPHIX_elem_M(PVd, n, m));
          QPHIX_c_eq_c_times_c(ctmp, QPHIX_elem_M(QVd, i, j), QPHIX_elem_M(RVd, n, m));
          QPHIX_c_peq_c(der, ctmp);
          QPHIX_c_eq_c_times_c(ctmp, QPHIX_elem_M(Vd, i, m), QPHIX_elem_M(Vd, n, j));
          QPHIX_c_peq_r_times_c(der, f1, ctmp);
          QPHIX_c_eq_c_times_c(ctmp, QPHIX_elem_M(QQVd, i, j), QPHIX_elem_M(SVd, n, m));
          QPHIX_c_peq_c(der, ctmp);
          QPHIX_c_eq_c_times_c(ctmp, QPHIX_elem_M(Vd, i, m), QPHIX_elem_M(QVd, n, j));
          QPHIX_c_eq_c_times_c(ctmp2, QPHIX_elem_M(Vd, n, j), QPHIX_elem_M(QVd, i, m));
          QPHIX_c_peq_r_times_c(der, f2, ctmp);
          QPHIX_c_peq_r_times_c(der, f2, ctmp2);
          dwdagdv->t4[i][m][n][j].real = QPHIX_real(der);
          dwdagdv->t4[i][m][n][j].imag = QPHIX_imag(der);
          nflops += 10*6+5*2+4*4;
        }
      }
    }
  }

  info->final_flop += nflops;
}


/* **************************************************
   SVD stuff
 ************************************************** */
/* Singular value decomposition (SVD) for
   3x3 complex matrix
   A.Bazavov, Feb 20 2009


   Algorithm sketch:

   SVD is performed in two steps:
   1) 3x3 complex matrix is reduced to real bidiagonal form
   with Householder transformations,
   3x3 matrix requires three left and two right such
   transformations
   2) bidiagonal matrix has the form
   [ b00 b01   0 ]
   [   0 b11 b12 ]
   [   0   0 b22 ]
   It is iteratively diagonalized with QR algorithm with shifts
   (it constructs Given rotations).
   There are many special cases (such as b00==0, b11==0, etc.)
   that are handled separately. If b01==0 then an auxiliary
   routine QPHIX_svd2x2bidiag is used to decompose the lower 2x2 block,
   if b12==0 the same is done for upper block.
   QPHIX_svd2x2bidiag is a separate routine because there are special
   cases for 2x2 superdiagonal matrix that need to be handled.

   This routine needs to be stable for singular matrices. Therefore, most
   of the operations are done to avoid underflow/overflow, for example,
   if norm=sqrt(a^2+b^2) then the calculation proceeds as:
   min=min(a,b), max=max(a,b)
   norm=max*sqrt(1+(min/max)^2)
   and so on.
   */

/* debugging define: prints a lot(!) */
/*#define QPHIX_SVD3x3_DEBUG*/

/* define precision for chopping small values */

/* defines that allow to remap input arrays easily */


/* Input: A -- 3x3 complex matrix,
Output: sigma[3] -- singular values,
U,V -- U(3) matrices such, that
A=U Sigma V^+ */
int QPHIX_svd3x3(QPHIX_info_t *info, qphix_su3_matrix *A, QPHIX_Real *sigma, 
    qphix_su3_matrix *U, qphix_su3_matrix *V) {
#define b00 P[0][0][0]
#define b01 P[0][1][0]
#define b02 P[0][2][0]
#define b10 P[1][0][0]
#define b11 P[1][1][0]
#define b12 P[1][2][0]
#define b20 P[2][0][0]
#define b21 P[2][1][0]
#define b22 P[2][2][0]

  QPHIX_Real Ad[3][3][2], P[3][3][2], Q[3][3][2];
  QPHIX_Real U1[3][3][2], U2[3][3][2], U3[3][3][2], V1[3][3][2], V2[3][3][2];
  QPHIX_Real UO3[3][3], VO3[3][3], v[3][2];
  QPHIX_Real UO2[2][2], VO2[2][2];
  register QPHIX_Real a, b, c, d, factor, norm, min, max, taure, tauim, beta;
  register QPHIX_Real m11, m12, m22, dm, lambdamax, cosphi, sinphi, tanphi, cotphi;
  register int i, j, iter;
  size_t nflops = 0;


  /* format of external matrices A, U and V can be arbitrary,
     therefore this routine uses defines (above) to access them
     and never reads A, U and V directly */

  /* original matrix can be in single precision,
     so copy it into double */
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      Ad[i][j][0] = QPHIX_real(QPHIX_elem_M(*A,i,j));
      Ad[i][j][1] = QPHIX_imag(QPHIX_elem_M(*A,i,j));
    }
  }


  i=0; j=0;

  /* *** Step 1: build first left reflector v,
     calculate first left rotation U1,
     apply to the original matrix A *** */
  /* calculate norm of ( A[10] )
     ( A[20] ) vector
     with minimal loss of accuracy (similar to BLAS) */
  c = 1.;
  factor = fabs( Ad[1][0][0] );
  a = fabs( Ad[1][0][1] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + (factor/a)*(factor/a);
      factor = a;
    }
    else {
      c = 1 + (a/factor)*(a/factor);
    }
  }
  a = fabs( Ad[2][0][0] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + c*(factor/a)*(factor/a);
      factor = a;
    }
    else {
      c += (a/factor)*(a/factor);
    }
  }
  a = fabs( Ad[2][0][1] );
  if( a!=0 ) {
    if( factor < a ) {
      c = 1 + c*(factor/a)*(factor/a);
      factor = a;
    }
    else {
      c += (a/factor)*(a/factor);
    }
  }
  norm = factor*sqrt(c);

  nflops += 15;

  if( norm==0 && Ad[0][0][1]==0 ) { /* no rotation needed */
#ifdef QPHIX_SVD3x3_DEBUG
    printf("Step 1: no rotation needed\n");
#endif /* QPHIX_SVD3x3_DEBUG */
    U1[0][0][0]=1.; U1[0][0][1]=0.;
    U1[0][1][0]=0.; U1[0][1][1]=0.;
    U1[0][2][0]=0.; U1[0][2][1]=0.;
    U1[1][0][0]=0.; U1[1][0][1]=0.;
    U1[1][1][0]=1.; U1[1][1][1]=0.;
    U1[1][2][0]=0.; U1[1][2][1]=0.;
    U1[2][0][0]=0.; U1[2][0][1]=0.;
    U1[2][1][0]=0.; U1[2][1][1]=0.;
    U1[2][2][0]=1.; U1[2][2][1]=0.;
    P[0][0][0]=Ad[0][0][0]; P[0][0][1]=Ad[0][0][1];
    P[1][0][0]=Ad[1][0][0]; P[1][0][1]=Ad[1][0][1];
    P[2][0][0]=Ad[2][0][0]; P[2][0][1]=Ad[2][0][1];
    P[0][1][0]=Ad[0][1][0]; P[0][1][1]=Ad[0][1][1];
    P[1][1][0]=Ad[1][1][0]; P[1][1][1]=Ad[1][1][1];
    P[2][1][0]=Ad[2][1][0]; P[2][1][1]=Ad[2][1][1];
    P[0][2][0]=Ad[0][2][0]; P[0][2][1]=Ad[0][2][1];
    P[1][2][0]=Ad[1][2][0]; P[1][2][1]=Ad[1][2][1];
    P[2][2][0]=Ad[2][2][0]; P[2][2][1]=Ad[2][2][1];
  }
  else {


    /* get the norm of full first column of A matrix */
    c=1.;
    factor = norm;
    a = fabs( Ad[0][0][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( Ad[0][0][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of first column */
    if( Ad[0][0][0]>0 ) {
      beta = -beta;
    }

#ifdef QPHIX_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QPHIX_SVD3x3_DEBUG */


    /* a=Re(A_00-beta), b=Im(A_00-beta) */
    a=Ad[0][0][0]-beta; b=Ad[0][0][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;


    /* construct reflector (vector "v" for Householder transformation)
       v_0=1 */
    v[0][0]=1.; v[0][1]=0.;
    /* v_1=A_10/(A_00-beta)=A_10/(a+ib)=(A_10*(a-ib))/norm^2=(A_10/norm)*((a-ib)/norm)
       =(A_10/norm)*(c-id)=|a=Re(A_10)/norm,b=Im(A_10)/norm|=(a+ib)*(c-id)
       =(a*c+b*d)+i(b*c-a*d) */
    a=Ad[1][0][0]/norm; b=Ad[1][0][1]/norm;
    v[1][0]=a*c+b*d;
    v[1][1]=b*c-a*d;
    /* v_2=A_20/(A_00-beta)=A_20/(a+ib)=(A_20*(a-ib))/norm^2=(A_20/norm)*((a-ib)/norm)
       =(A_20/norm)*(c-id)=|a=Re(A_20)/norm,b=Im(A_20)/norm|=(a+ib)*(c-id)
       =(a*c+b*d)+i(b*c-a*d) */
    a=Ad[2][0][0]/norm; b=Ad[2][0][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;
#ifdef QPHIX_SVD3x3_DEBUG
    for(i=0;i<3;i++) {
      printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
    }
#endif /* QPHIX_SVD3x3_DEBUG */

    /* calcualate tau (coefficient for reflector) */
    taure=(beta-Ad[0][0][0])/beta;
    tauim=(Ad[0][0][1])/beta;


    /* assemble left unitary matrix U1=I-tau^+*v*v^+ (store in U1[3][3][2])
       U1_00=A_00/beta */
    U1[0][0][0]=(Ad[0][0][0])/beta;
    U1[0][0][1]=(Ad[0][0][1])/beta;
    /* U1_10=A_10/beta */
    U1[1][0][0]=(Ad[1][0][0])/beta;
    U1[1][0][1]=(Ad[1][0][1])/beta;
    /* U1_20=A_20/beta */
    U1[2][0][0]=(Ad[2][0][0])/beta;
    U1[2][0][1]=(Ad[2][0][1])/beta;
    /* U1_01=-tau^+*v_1^+=-(tau*v_1)^+ */
    U1[0][1][0]=-(taure*v[1][0]-tauim*v[1][1]);
    U1[0][1][1]=taure*v[1][1]+tauim*v[1][0];
    /* U1_11=1-tau^+*v_1*v_1^+ */
    a=v[1][0]*v[1][0]+v[1][1]*v[1][1];
    U1[1][1][0]=1-taure*a;
    U1[1][1][1]=tauim*a;
    /* U1_21=-tau^+*v_2*v_1^+ */
    /* v_2*v_1^+ */
    a=v[2][0]*v[1][0]+v[2][1]*v[1][1];
    b=-v[2][0]*v[1][1]+v[2][1]*v[1][0];
    U1[2][1][0]=-(taure*a+tauim*b);
    U1[2][1][1]=-(taure*b-tauim*a);
    /* U1_02=-tau^+*v_2^+=-(tau*v_2)^+ */
    U1[0][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    U1[0][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* U1_12=-tau^+*v_1*v_2^+ */
    /* v_1*v_2^+ */
    a=v[1][0]*v[2][0]+v[1][1]*v[2][1];
    b=-v[1][0]*v[2][1]+v[1][1]*v[2][0];
    U1[1][2][0]=-(taure*a+tauim*b);
    U1[1][2][1]=-(taure*b-tauim*a);
    /* U1_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    U1[2][2][0]=1-taure*a;
    U1[2][2][1]=tauim*a;

    nflops += 91;


    /* apply the transformation to A matrix and store the result in P
       P=U^+A */
    P[0][0][0]=beta;
    P[0][0][1]=0;
    P[1][0][0]=0;
    P[1][0][1]=0;
    P[2][0][0]=0;
    P[2][0][1]=0;
    /* P_01=U1_00^+*A_01+U1_10^+*A_11+U1_20^+*A_21 */
    P[0][1][0]=U1[0][0][0]*Ad[0][1][0]+U1[0][0][1]*Ad[0][1][1]
      +U1[1][0][0]*Ad[1][1][0]+U1[1][0][1]*Ad[1][1][1]
      +U1[2][0][0]*Ad[2][1][0]+U1[2][0][1]*Ad[2][1][1];
    P[0][1][1]=U1[0][0][0]*Ad[0][1][1]-U1[0][0][1]*Ad[0][1][0]
      +U1[1][0][0]*Ad[1][1][1]-U1[1][0][1]*Ad[1][1][0]
      +U1[2][0][0]*Ad[2][1][1]-U1[2][0][1]*Ad[2][1][0];
    /* P_02=U1_00^+*A_02+U1_10^+*A_12+U1_20^+*A_22 */
    P[0][2][0]=U1[0][0][0]*Ad[0][2][0]+U1[0][0][1]*Ad[0][2][1]
      +U1[1][0][0]*Ad[1][2][0]+U1[1][0][1]*Ad[1][2][1]
      +U1[2][0][0]*Ad[2][2][0]+U1[2][0][1]*Ad[2][2][1];
    P[0][2][1]=U1[0][0][0]*Ad[0][2][1]-U1[0][0][1]*Ad[0][2][0]
      +U1[1][0][0]*Ad[1][2][1]-U1[1][0][1]*Ad[1][2][0]
      +U1[2][0][0]*Ad[2][2][1]-U1[2][0][1]*Ad[2][2][0];
    /* P_11=U1_01^+*A_01+U1_11^+*A_11+U1_21^+*A_21 */
    P[1][1][0]=U1[0][1][0]*Ad[0][1][0]+U1[0][1][1]*Ad[0][1][1]
      +U1[1][1][0]*Ad[1][1][0]+U1[1][1][1]*Ad[1][1][1]
      +U1[2][1][0]*Ad[2][1][0]+U1[2][1][1]*Ad[2][1][1];
    P[1][1][1]=U1[0][1][0]*Ad[0][1][1]-U1[0][1][1]*Ad[0][1][0]
      +U1[1][1][0]*Ad[1][1][1]-U1[1][1][1]*Ad[1][1][0]
      +U1[2][1][0]*Ad[2][1][1]-U1[2][1][1]*Ad[2][1][0];
    /* P_12=U1_01^+*A_02+U1_11^+*A_12+U1_21^+*A_22 */
    P[1][2][0]=U1[0][1][0]*Ad[0][2][0]+U1[0][1][1]*Ad[0][2][1]
      +U1[1][1][0]*Ad[1][2][0]+U1[1][1][1]*Ad[1][2][1]
      +U1[2][1][0]*Ad[2][2][0]+U1[2][1][1]*Ad[2][2][1];
    P[1][2][1]=U1[0][1][0]*Ad[0][2][1]-U1[0][1][1]*Ad[0][2][0]
      +U1[1][1][0]*Ad[1][2][1]-U1[1][1][1]*Ad[1][2][0]
      +U1[2][1][0]*Ad[2][2][1]-U1[2][1][1]*Ad[2][2][0];
    /* P_21=U1_02^+*A_01+U1_12^+*A_11+U1_22^+*A_21 */
    P[2][1][0]=U1[0][2][0]*Ad[0][1][0]+U1[0][2][1]*Ad[0][1][1]
      +U1[1][2][0]*Ad[1][1][0]+U1[1][2][1]*Ad[1][1][1]
      +U1[2][2][0]*Ad[2][1][0]+U1[2][2][1]*Ad[2][1][1];
    P[2][1][1]=U1[0][2][0]*Ad[0][1][1]-U1[0][2][1]*Ad[0][1][0]
      +U1[1][2][0]*Ad[1][1][1]-U1[1][2][1]*Ad[1][1][0]
      +U1[2][2][0]*Ad[2][1][1]-U1[2][2][1]*Ad[2][1][0];
    /* P_22=U1_02^+*A_02+U1_12^+*A_12+U1_22^+*A_22 */
    P[2][2][0]=U1[0][2][0]*Ad[0][2][0]+U1[0][2][1]*Ad[0][2][1]
      +U1[1][2][0]*Ad[1][2][0]+U1[1][2][1]*Ad[1][2][1]
      +U1[2][2][0]*Ad[2][2][0]+U1[2][2][1]*Ad[2][2][1];
    P[2][2][1]=U1[0][2][0]*Ad[0][2][1]-U1[0][2][1]*Ad[0][2][0]
      +U1[1][2][0]*Ad[1][2][1]-U1[1][2][1]*Ad[1][2][0]
      +U1[2][2][0]*Ad[2][2][1]-U1[2][2][1]*Ad[2][2][0];

    nflops += 9*12;

  }
#ifdef QPHIX_SVD3x3_DEBUG
  printf("Left unitary matrix U1:\n");
  for(i=0;i<3;i++)for(j=0;j<3;j++) {
    printf( "U1[%d][%d].re=%26.18e  U1[%d][%d].im=%26.18e\n",
        i, j, U1[i][j][0], i, j, U1[i][j][1] );
  }
#endif /* QPHIX_SVD3x3_DEBUG */



  /* *** Step 2: build first right reflector v,
     calculate first right rotation V1,
     apply to the matrix P from step 1 *** */
  /* calculate norm of ( P[02] )
     with minimal loss of accuracy */
  a=fabs( P[0][2][0] ); b=fabs( P[0][2][1] );
  /* norm=sqrt(a^2+b^2) */
  if( a>b ) {
    max=a; min=b;
  }
  else {
    max=b; min=a;
  }
  if( min==0 ) {
    norm = max;
  }
  else {
    c = min/max;
    norm = max*sqrt(1+c*c);
  }

  if( norm==0 && P[0][1][1]==0 ) { /* no rotation needed */
#ifdef QPHIX_SVD3x3_DEBUG
    printf("Step 2: no rotation needed\n");
#endif /* QPHIX_SVD3x3_DEBUG */
    V1[0][0][0]=1.; V1[0][0][1]=0.;
    V1[0][1][0]=0.; V1[0][1][1]=0.;
    V1[0][2][0]=0.; V1[0][2][1]=0.;
    V1[1][0][0]=0.; V1[1][0][1]=0.;
    V1[1][1][0]=1.; V1[1][1][1]=0.;
    V1[1][2][0]=0.; V1[1][2][1]=0.;
    V1[2][0][0]=0.; V1[2][0][1]=0.;
    V1[2][1][0]=0.; V1[2][1][1]=0.;
    V1[2][2][0]=1.; V1[2][2][1]=0.;
    Q[0][0][0]=P[0][0][0]; Q[0][0][1]=P[0][0][1];
    Q[1][0][0]=P[1][0][0]; Q[1][0][1]=P[1][0][1];
    Q[2][0][0]=P[2][0][0]; Q[2][0][1]=P[2][0][1];
    Q[0][1][0]=P[0][1][0]; Q[0][1][1]=P[0][1][1];
    Q[1][1][0]=P[1][1][0]; Q[1][1][1]=P[1][1][1];
    Q[2][1][0]=P[2][1][0]; Q[2][1][1]=P[2][1][1];
    Q[0][2][0]=P[0][2][0]; Q[0][2][1]=P[0][2][1];
    Q[1][2][0]=P[1][2][0]; Q[1][2][1]=P[1][2][1];
    Q[2][2][0]=P[2][2][0]; Q[2][2][1]=P[2][2][1];
  }
  else {
    /* get the norm of (P_01 P_02) row vector */
    c=1.;
    factor = norm;
    a = fabs( P[0][1][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( P[0][1][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of (P_01 P_02) row vector */
    if( P[0][1][0]>0 ) {
      beta = -beta;
    }

#ifdef QPHIX_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QPHIX_SVD3x3_DEBUG */


    /* a=Re(P_01^+-beta), b=Im(P_01^+-beta) */
    a=P[0][1][0]-beta; b=-P[0][1][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;


    /* construct reflector (vector "v" for Householder transformation) */
    /* v_0=0 */
    v[0][0]=0.; v[0][1]=0.;
    /* v_1=1 */
    v[1][0]=1.; v[1][1]=0.;
    /* v_2=P_02^+/(P_01^+-beta)=P_02^+/(a+ib)=(P_02^+*(a-ib))/norm^2=(P_02^+/norm)*((a-ib)/norm)
       =(P_02^+/norm)*(c-id)=|a=Re(P_02^+)/norm,b=Im(P_02^+)/norm|=(a+ib)*(c-id)
       =(a*c+b*d)+i(b*c-a*d) */
    a=P[0][2][0]/norm; b=-P[0][2][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;

    nflops += 27;

#ifdef QPHIX_SVD3x3_DEBUG
    for(i=0;i<3;i++) {
      printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
    }
#endif /* QPHIX_SVD3x3_DEBUG */

    /* calcualate tau (coefficient for reflector) */
    taure=(beta-P[0][1][0])/beta;
    tauim=-P[0][1][1]/beta;

    /* assemble right unitary matrix V1=I-tau^+*v*v^+ (store in V1[3][3][2]) */
    V1[0][0][0]=1.;
    V1[0][0][1]=0.;
    V1[1][0][0]=0.;
    V1[1][0][1]=0.;
    V1[2][0][0]=0.;
    V1[2][0][1]=0.;
    V1[0][1][0]=0.;
    V1[0][1][1]=0.;
    V1[0][2][0]=0.;
    V1[0][2][1]=0.;
    /* V1_11=P_01^+/beta */
    V1[1][1][0]=P[0][1][0]/beta;
    V1[1][1][1]=-P[0][1][1]/beta;
    /* V1_21=P_02^+/beta */
    V1[2][1][0]=P[0][2][0]/beta;
    V1[2][1][1]=-P[0][2][1]/beta;
    /* V1_12=-tau^+*v_2^+=-(tau*v_2)^+ */
    V1[1][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    V1[1][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* V1_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    V1[2][2][0]=1-taure*a;
    V1[2][2][1]=tauim*a;


    /* apply the transformation to P matrix and store the result in Q
       Q=PV */
    Q[0][0][0]=P[0][0][0];
    Q[0][0][1]=0.;
    Q[1][0][0]=0.;
    Q[1][0][1]=0.;
    Q[2][0][0]=0.;
    Q[2][0][1]=0.;
    Q[0][1][0]=beta;
    Q[0][1][1]=0.;
    Q[0][2][0]=0.;
    Q[0][2][1]=0.;
    /* Q_11=P_11*V1_11+P_12*V_21 */
    Q[1][1][0]=P[1][1][0]*V1[1][1][0]-P[1][1][1]*V1[1][1][1]
      +P[1][2][0]*V1[2][1][0]-P[1][2][1]*V1[2][1][1];
    Q[1][1][1]=P[1][1][0]*V1[1][1][1]+P[1][1][1]*V1[1][1][0]
      +P[1][2][0]*V1[2][1][1]+P[1][2][1]*V1[2][1][0];
    /* Q_12=P_11*V1_12+P_12*V_22 */
    Q[1][2][0]=P[1][1][0]*V1[1][2][0]-P[1][1][1]*V1[1][2][1]
      +P[1][2][0]*V1[2][2][0]-P[1][2][1]*V1[2][2][1];
    Q[1][2][1]=P[1][1][0]*V1[1][2][1]+P[1][1][1]*V1[1][2][0]
      +P[1][2][0]*V1[2][2][1]+P[1][2][1]*V1[2][2][0];
    /* Q_21=P_21*V1_11+P_22*V_21 */
    Q[2][1][0]=P[2][1][0]*V1[1][1][0]-P[2][1][1]*V1[1][1][1]
      +P[2][2][0]*V1[2][1][0]-P[2][2][1]*V1[2][1][1];
    Q[2][1][1]=P[2][1][0]*V1[1][1][1]+P[2][1][1]*V1[1][1][0]
      +P[2][2][0]*V1[2][1][1]+P[2][2][1]*V1[2][1][0];
    /* Q_22=P_21*V1_12+P_22*V_22 */
    Q[2][2][0]=P[2][1][0]*V1[1][2][0]-P[2][1][1]*V1[1][2][1]
      +P[2][2][0]*V1[2][2][0]-P[2][2][1]*V1[2][2][1];
    Q[2][2][1]=P[2][1][0]*V1[1][2][1]+P[2][1][1]*V1[1][2][0]
      +P[2][2][0]*V1[2][2][1]+P[2][2][1]*V1[2][2][0];

    nflops += 15 + 7*8;
  }
#ifdef QPHIX_SVD3x3_DEBUG
  printf("Right unitary matrix V1:\n");
  for(i=0;i<3;i++)for(j=0;j<3;j++) {
    printf( "V1[%d][%d].re=%26.18e  V1[%d][%d].im=%26.18e\n",
        i, j, V1[i][j][0], i, j, V1[i][j][1] );
  }
#endif /* QPHIX_SVD3x3_DEBUG */



  /* *** Step 3: build second left reflector v,
     calculate second left rotation U2,
     apply to the matrix Q *** */
  /* calculate norm of ( Q[21] )
     with minimal loss of accuracy (similar to BLAS) */
  c=fabs(Q[2][1][0]); d=fabs(Q[2][1][1]);
  if( c>d ) {
    max=c; min=d;
  }
  else {
    max=d; min=c;
  }
  if( min==0 ) {
    norm = max;
  }
  else {
    c = min/max;
    norm = max*sqrt(1+c*c);
  }

  if( norm==0 && Q[1][1][1]==0 ) { /* no rotation needed */
#ifdef QPHIX_SVD3x3_DEBUG
    printf("Step 3: no rotation needed\n");
#endif /* QPHIX_SVD3x3_DEBUG */
    U2[0][0][0]=1.; U2[0][0][1]=0.;
    U2[0][1][0]=0.; U2[0][1][1]=0.;
    U2[0][2][0]=0.; U2[0][2][1]=0.;
    U2[1][0][0]=0.; U2[1][0][1]=0.;
    U2[1][1][0]=1.; U2[1][1][1]=0.;
    U2[1][2][0]=0.; U2[1][2][1]=0.;
    U2[2][0][0]=0.; U2[2][0][1]=0.;
    U2[2][1][0]=0.; U2[2][1][1]=0.;
    U2[2][2][0]=1.; U2[2][2][1]=0.;
    P[0][0][0]=Q[0][0][0]; P[0][0][1]=Q[0][0][1];
    P[1][0][0]=Q[1][0][0]; P[1][0][1]=Q[1][0][1];
    P[2][0][0]=Q[2][0][0]; P[2][0][1]=Q[2][0][1];
    P[0][1][0]=Q[0][1][0]; P[0][1][1]=Q[0][1][1];
    P[1][1][0]=Q[1][1][0]; P[1][1][1]=Q[1][1][1];
    P[2][1][0]=Q[2][1][0]; P[2][1][1]=Q[2][1][1];
    P[0][2][0]=Q[0][2][0]; P[0][2][1]=Q[0][2][1];
    P[1][2][0]=Q[1][2][0]; P[1][2][1]=Q[1][2][1];
    P[2][2][0]=Q[2][2][0]; P[2][2][1]=Q[2][2][1];
  }
  else {
    /* get the norm of (Q_11 Q_21) column vector */
    c=1.;
    factor = norm;
    a = fabs( Q[1][1][0] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + (factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    a = fabs( Q[1][1][1] );
    if( a!=0 ) {
      if( factor < a ) {
        c = 1 + c*(factor/a)*(factor/a);
        factor = a;
      }
      else {
        c += (a/factor)*(a/factor);
      }
    }
    beta = factor*sqrt(c); /* norm of (Q_11 Q_21) column vector */
    if( Q[1][1][0]>0 ) {
      beta = -beta;
    }

#ifdef QPHIX_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QPHIX_SVD3x3_DEBUG */


    /* a=Re(Q_11-beta), b=Im(Q_11-beta) */
    a=Q[1][1][0]-beta; b=Q[1][1][1];
    /* norm=sqrt(a^2+b^2) */
    c=fabs(a); d=fabs(b);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      norm = max;
    }
    else {
      c = min/max;
      norm = max*sqrt(1+c*c);
    }
    /* c=a/norm, d=b/norm */
    c=a/norm; d=b/norm;

    /* construct reflector (vector "v" for Householder transformation) */
    /* v_0=0 */
    v[0][0]=0.; v[0][1]=0.;
    /* v_1=1 */
    v[1][0]=1.; v[1][1]=0.;
    /* v_2=Q_21/(Q_11-beta)=Q_21/(a+ib)=(Q_21*(a-ib))/norm^2=(Q_21/norm)*((a-ib)/norm)
       =(Q_21/norm)*(c-id)=|a=Re(Q_21)/norm,b=Im(Q_21)/norm|=(a+ib)*(c-id)
       =(a*c+b*d)+i(b*c-a*d) */
    a=Q[2][1][0]/norm; b=Q[2][1][1]/norm;
    v[2][0]=a*c+b*d;
    v[2][1]=b*c-a*d;

    nflops += 27;
#ifdef QPHIX_SVD3x3_DEBUG
    for(i=0;i<3;i++) {
      printf("v[%d].re=%28.18e  v[%d].im=%28.18e\n",i,v[i][0],i,v[i][1]);
    }
#endif /* QPHIX_SVD3x3_DEBUG */


    /* calcualate tau (coefficient for reflector) */
    taure=(beta-Q[1][1][0])/beta;
    tauim=Q[1][1][1]/beta;


    /* assemble right unitary matrix U2=I-tau^+*v*v^+ (store in U2[3][3][2]) */
    U2[0][0][0]=1.;
    U2[0][0][1]=0.;
    U2[1][0][0]=0.;
    U2[1][0][1]=0.;
    U2[2][0][0]=0.;
    U2[2][0][1]=0.;
    U2[0][1][0]=0.;
    U2[0][1][1]=0.;
    U2[0][2][0]=0.;
    U2[0][2][1]=0.;
    /* U2_11=Q_11/beta */
    U2[1][1][0]=Q[1][1][0]/beta;
    U2[1][1][1]=Q[1][1][1]/beta;
    /* U2_21=Q_21/beta */
    U2[2][1][0]=Q[2][1][0]/beta;
    U2[2][1][1]=Q[2][1][1]/beta;
    /* U2_12=-tau^+*v_2^+=-(tau*v_2)^+ */
    U2[1][2][0]=-(taure*v[2][0]-tauim*v[2][1]);
    U2[1][2][1]=taure*v[2][1]+tauim*v[2][0];
    /* U2_22=1-tau^+*v_2*v_2^+ */
    a=v[2][0]*v[2][0]+v[2][1]*v[2][1];
    U2[2][2][0]=1-taure*a;
    U2[2][2][1]=tauim*a;
#ifdef QPHIX_SVD3x3_DEBUG
    printf("Left unitary matrix U2:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "U2[%d][%d].re=%26.18e  U2[%d][%d].im=%26.18e\n",
          i, j, U2[i][j][0], i, j, U2[i][j][1] );
    }
#endif /* QPHIX_SVD3x3_DEBUG */


    /* apply the transformation to Q matrix and store the result in P
       P=U^+Q */
    P[0][0][0]=Q[0][0][0];
    P[0][0][1]=0.;
    P[1][0][0]=0.;
    P[1][0][1]=0.;
    P[2][0][0]=0.;
    P[2][0][1]=0.;
    P[0][1][0]=Q[0][1][0];
    P[0][1][1]=0.;
    P[0][2][0]=0.;
    P[0][2][1]=0.;
    P[1][1][0]=beta;
    P[1][1][1]=0.;
    P[2][1][0]=0.;
    P[2][1][1]=0.;
    /* P_12=U2_11^+*Q_12+U2_21^+*Q_22 */
    P[1][2][0]=U2[1][1][0]*Q[1][2][0]+U2[1][1][1]*Q[1][2][1]
      +U2[2][1][0]*Q[2][2][0]+U2[2][1][1]*Q[2][2][1];
    P[1][2][1]=U2[1][1][0]*Q[1][2][1]-U2[1][1][1]*Q[1][2][0]
      +U2[2][1][0]*Q[2][2][1]-U2[2][1][1]*Q[2][2][0];
    /* P_22=U2_12^+*Q_12+U2_22^+*Q_22 */
    P[2][2][0]=U2[1][2][0]*Q[1][2][0]+U2[1][2][1]*Q[1][2][1]
      +U2[2][2][0]*Q[2][2][0]+U2[2][2][1]*Q[2][2][1];
    P[2][2][1]=U2[1][2][0]*Q[1][2][1]-U2[1][2][1]*Q[1][2][0]
      +U2[2][2][0]*Q[2][2][1]-U2[2][2][1]*Q[2][2][0];

    nflops += 15 + 7*8;

  }



  /* *** Step 4: build second right reflector v,
     calculate second right rotation V2,
     apply to the matrix P *** */
  if( P[1][2][1]==0 ) { /* no rotation needed */
#ifdef QPHIX_SVD3x3_DEBUG
    printf("Step 4: no rotation needed\n");
#endif /* QPHIX_SVD3x3_DEBUG */
    V2[0][0][0]=1.; V2[0][0][1]=0.;
    V2[0][1][0]=0.; V2[0][1][1]=0.;
    V2[0][2][0]=0.; V2[0][2][1]=0.;
    V2[1][0][0]=0.; V2[1][0][1]=0.;
    V2[1][1][0]=1.; V2[1][1][1]=0.;
    V2[1][2][0]=0.; V2[1][2][1]=0.;
    V2[2][0][0]=0.; V2[2][0][1]=0.;
    V2[2][1][0]=0.; V2[2][1][1]=0.;
    V2[2][2][0]=1.; V2[2][2][1]=0.;
    Q[0][0][0]=P[0][0][0]; Q[0][0][1]=P[0][0][1];
    Q[1][0][0]=P[1][0][0]; Q[1][0][1]=P[1][0][1];
    Q[2][0][0]=P[2][0][0]; Q[2][0][1]=P[2][0][1];
    Q[0][1][0]=P[0][1][0]; Q[0][1][1]=P[0][1][1];
    Q[1][1][0]=P[1][1][0]; Q[1][1][1]=P[1][1][1];
    Q[2][1][0]=P[2][1][0]; Q[2][1][1]=P[2][1][1];
    Q[0][2][0]=P[0][2][0]; Q[0][2][1]=P[0][2][1];
    Q[1][2][0]=P[1][2][0]; Q[1][2][1]=P[1][2][1];
    Q[2][2][0]=P[2][2][0]; Q[2][2][1]=P[2][2][1];
  }
  else {
    /* calculate norm of ( P[12] ) */
    c=fabs(P[1][2][0]); d=fabs(P[1][2][1]);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      beta = max;
    }
    else {
      c = min/max;
      beta = max*sqrt(1+c*c);
    }

    if( P[1][2][0]>0 ) {
      beta = -beta;
    }

#ifdef QPHIX_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QPHIX_SVD3x3_DEBUG */

    /* assemble right unitary matrix V1=I-tau^+*v*v^+ (store in V1[3][3][2]) */
    V2[0][0][0]=1.;
    V2[0][0][1]=0.;
    V2[1][0][0]=0.;
    V2[1][0][1]=0.;
    V2[2][0][0]=0.;
    V2[2][0][1]=0.;
    V2[0][1][0]=0.;
    V2[0][1][1]=0.;
    V2[0][2][0]=0.;
    V2[0][2][1]=0.;
    V2[1][1][0]=1.;
    V2[1][1][1]=0.;
    V2[2][1][0]=0.;
    V2[2][1][1]=0.;
    V2[1][2][0]=0.;
    V2[1][2][1]=0.;
    /* V2_22=1-tau^+*v_2*v_2^+=1-tau^+ */
    V2[2][2][0]=P[1][2][0]/beta;
    V2[2][2][1]=-P[1][2][1]/beta;
#ifdef QPHIX_SVD3x3_DEBUG
    printf("Right unitary matrix V2:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "V2[%d][%d].re=%26.18e  V2[%d][%d].im=%26.18e\n",
          i, j, V2[i][j][0], i, j, V2[i][j][1] );
    }
#endif /* QPHIX_SVD3x3_DEBUG */


    /* apply the transformation to P matrix and store the result in Q
       Q=PV */
    Q[0][0][0]=P[0][0][0];
    Q[0][0][1]=0.;
    Q[1][0][0]=0.;
    Q[1][0][1]=0.;
    Q[2][0][0]=0.;
    Q[2][0][1]=0.;
    Q[0][1][0]=P[0][1][0];
    Q[0][1][1]=0.;
    Q[0][2][0]=0.;
    Q[0][2][1]=0.;
    Q[1][1][0]=P[1][1][0];
    Q[1][1][1]=0.;
    Q[1][2][0]=beta;
    Q[1][2][1]=0.;
    Q[2][1][0]=0.;
    Q[2][1][1]=0.;
    /* Q_22=P_22*V2_22 */
    Q[2][2][0]=P[2][2][0]*V2[2][2][0]-P[2][2][1]*V2[2][2][1];
    Q[2][2][1]=P[2][2][0]*V2[2][2][1]+P[2][2][1]*V2[2][2][0];

    nflops += 12;
  }



  /* *** Step 5: build third left reflector v,
     calculate third left rotation U3,
     apply to the matrix P *** */
  if( Q[2][2][1]==0 ) { /* no rotation needed */
#ifdef QPHIX_SVD3x3_DEBUG
    printf("Step 5: no rotation needed\n");
#endif /* QPHIX_SVD3x3_DEBUG */
    U3[0][0][0]=1.; U3[0][0][1]=0.;
    U3[0][1][0]=0.; U3[0][1][1]=0.;
    U3[0][2][0]=0.; U3[0][2][1]=0.;
    U3[1][0][0]=0.; U3[1][0][1]=0.;
    U3[1][1][0]=1.; U3[1][1][1]=0.;
    U3[1][2][0]=0.; U3[1][2][1]=0.;
    U3[2][0][0]=0.; U3[2][0][1]=0.;
    U3[2][1][0]=0.; U3[2][1][1]=0.;
    U3[2][2][0]=1.; U3[2][2][1]=0.;
    P[0][0][0]=Q[0][0][0]; P[0][0][1]=Q[0][0][1];
    P[1][0][0]=Q[1][0][0]; P[1][0][1]=Q[1][0][1];
    P[2][0][0]=Q[2][0][0]; P[2][0][1]=Q[2][0][1];
    P[0][1][0]=Q[0][1][0]; P[0][1][1]=Q[0][1][1];
    P[1][1][0]=Q[1][1][0]; P[1][1][1]=Q[1][1][1];
    P[2][1][0]=Q[2][1][0]; P[2][1][1]=Q[2][1][1];
    P[0][2][0]=Q[0][2][0]; P[0][2][1]=Q[0][2][1];
    P[1][2][0]=Q[1][2][0]; P[1][2][1]=Q[1][2][1];
    P[2][2][0]=Q[2][2][0]; P[2][2][1]=Q[2][2][1];
  }
  else {
    /* calculate norm of ( Q[22] ) */
    c=fabs(Q[2][2][0]); d=fabs(Q[2][2][1]);
    if( c>d ) {
      max=c; min=d;
    }
    else {
      max=d; min=c;
    }
    if( min==0 ) {
      beta = max;
    }
    else {
      c = min/max;
      beta = max*sqrt(1+c*c);
    }

    if( Q[2][2][0]>0 ) {
      beta = -beta;
    }

#ifdef QPHIX_SVD3x3_DEBUG
    printf("beta=%28.18e\n",beta);
#endif /* QPHIX_SVD3x3_DEBUG */

    /* assemble left unitary matrix U3=I-tau^+*v*v^+ (store in U3[3][3][2]) */
    U3[0][0][0]=1.;
    U3[0][0][1]=0.;
    U3[1][0][0]=0.;
    U3[1][0][1]=0.;
    U3[2][0][0]=0.;
    U3[2][0][1]=0.;
    U3[0][1][0]=0.;
    U3[0][1][1]=0.;
    U3[0][2][0]=0.;
    U3[0][2][1]=0.;
    U3[1][1][0]=1.;
    U3[1][1][1]=0.;
    U3[2][1][0]=0.;
    U3[2][1][1]=0.;
    U3[1][2][0]=0.;
    U3[1][2][1]=0.;
    /* U3_22=1-tau^+*v_2*v_2^+=1-tau^+ */
    U3[2][2][0]=Q[2][2][0]/beta;
    U3[2][2][1]=Q[2][2][1]/beta;
#ifdef QPHIX_SVD3x3_DEBUG
    printf("Left unitary matrix U3:\n");
    for(i=0;i<3;i++)for(j=0;j<3;j++) {
      printf( "U3[%d][%d].re=%26.18e  U3[%d][%d].im=%26.18e\n",
          i, j, U3[i][j][0], i, j, U3[i][j][1] );
    }
#endif /* QPHIX_SVD3x3_DEBUG */


    /* apply the transformation to Q matrix and store the result in P
       P=U^+Q */
    P[0][0][0]=Q[0][0][0];
    P[0][0][1]=0.;
    P[1][0][0]=0.;
    P[1][0][1]=0.;
    P[2][0][0]=0.;
    P[2][0][1]=0.;
    P[0][1][0]=Q[0][1][0];
    P[0][1][1]=0.;
    P[0][2][0]=0.;
    P[0][2][1]=0.;
    P[1][1][0]=Q[1][1][0];
    P[1][1][1]=0.;
    P[1][2][0]=Q[1][2][0];
    P[1][2][1]=0.;
    P[2][1][0]=0.;
    P[2][1][1]=0.;
    P[2][2][0]=beta;
    P[2][2][1]=0.;

    nflops += 6;

  }




  /* *** This part starts with a bidiagonal matrix and uses
     QR algorithm with shifts to eliminate the superdiagonal *** */
  /* prepare left and right real orthogonal matrices that
     accumulate Givens rotations from QR algorithm */
  UO3[0][0]=1.; UO3[0][1]=0.; UO3[0][2]=0.;
  UO3[1][0]=0.; UO3[1][1]=1.; UO3[1][2]=0.;
  UO3[2][0]=0.; UO3[2][1]=0.; UO3[2][2]=1.;
  VO3[0][0]=1.; VO3[0][1]=0.; VO3[0][2]=0.;
  VO3[1][0]=0.; VO3[1][1]=1.; VO3[1][2]=0.;
  VO3[2][0]=0.; VO3[2][1]=0.; VO3[2][2]=1.;

  iter=0;

#ifdef QPHIX_SVD3x3_DEBUG
  printf( "QR iteration: %d\n", iter );
  printf( "%+20.16e %+20.16e %+20.16e\n", b00, b01, b02 );
  printf( "%+20.16e %+20.16e %+20.16e\n", b10, b11, b12 );
  printf( "%+20.16e %+20.16e %+20.16e\n", b20, b21, b22 );
#endif /* QPHIX_SVD3x3_DEBUG */

  do {

    iter++;
    if(iter>300) return 1;

    /* chop small superdiagonal elements */
    if( fabs(b01) < QPHIX_SVD3x3_PREC*(fabs(b00)+fabs(b11)) ) {
      b01=0;
    }
    if( fabs(b12) < QPHIX_SVD3x3_PREC*(fabs(b00)+fabs(b22)) ) {
      b12=0;
    }

    nflops += 4;

    /* Cases:
       b01=b12=0 -- matrix is already diagonalized,
       b01=0 -- need to work with 2x2 lower block,
       b12=0 -- need to work with 2x2 upper block,
       else -- normal iteration */
    if( !(b01==0 && b12==0) ) {
      if( b01==0 ) {
#ifdef QPHIX_SVD3x3_DEBUG
        printf( "Entering case b01==0\n" );
#endif /* QPHIX_SVD3x3_DEBUG */
        /* need to diagonalize 2x2 lower block */
        QPHIX_svd2x2bidiag(info, &b11, &b12, &b22, UO2, VO2 );

        /* multiply left UO3 matrix */
        for(i=0;i<3;i++) {
          a=UO3[i][1]; b=UO3[i][2];
          UO3[i][1]=a*UO2[0][0]+b*UO2[1][0];
          UO3[i][2]=a*UO2[0][1]+b*UO2[1][1];
        }
        /* multiply right VO3 matrix */
        for(i=0;i<3;i++) {
          a=VO3[i][1]; b=VO3[i][2];
          VO3[i][1]=a*VO2[0][0]+b*VO2[1][0];
          VO3[i][2]=a*VO2[0][1]+b*VO2[1][1];
        }

        nflops += 36;

      }
      else {
        if( b12==0 ) {
#ifdef QPHIX_SVD3x3_DEBUG
          printf( "Entering case b12==0\n" );
#endif /* QPHIX_SVD3x3_DEBUG */
          /* need to diagonalize 2x2 upper block */
          QPHIX_svd2x2bidiag(info, &b00, &b01, &b11, UO2, VO2 );

          /* multiply left UO3 matrix */
          for(i=0;i<3;i++) {
            a=UO3[i][0]; b=UO3[i][1];
            UO3[i][0]=a*UO2[0][0]+b*UO2[1][0];
            UO3[i][1]=a*UO2[0][1]+b*UO2[1][1];
          }
          /* multiply right VO3 matrix */
          for(i=0;i<3;i++) {
            a=VO3[i][0]; b=VO3[i][1];
            VO3[i][0]=a*VO2[0][0]+b*VO2[1][0];
            VO3[i][1]=a*VO2[0][1]+b*VO2[1][1];
          }

          nflops += 36;
        }
        else {
          /* normal 3x3 iteration */

          /* QR shift does not work if there are zeros
             on the diagonal, therefore first check
             for special cases: b00==0 or b11==0 or b22==0 */

          if( b00==0 ) {
#ifdef QPHIX_SVD3x3_DEBUG
            printf( "Entering case b00==0\n" );
#endif /* QPHIX_SVD3x3_DEBUG */
            /* b01 can be rotated away to create b02,
               and then b02 can be rotated away
               (both are left rotations) */
            if( fabs(b01)>fabs(b11) ) {
              cotphi=b11/b01;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b01/b11;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][1];
              UO3[i][0]=a*cosphi-b*sinphi;
              UO3[i][1]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b02 */
            b11=b01*sinphi+b11*cosphi;
            b02=-b12*sinphi;
            b12=b12*cosphi;
            b01=0.;
            if( fabs(b02)>fabs(b22) ) {
              cotphi=b22/b02;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b02/b22;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][2];
              UO3[i][0]=a*cosphi-b*sinphi;
              UO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b22=b02*sinphi+b22*cosphi;
            b02=0.;

            nflops += 56;
          }
          else if( b11==0 ) {
#ifdef QPHIX_SVD3x3_DEBUG
            printf( "Entering case b11==0\n" );
#endif /* QPHIX_SVD3x3_DEBUG */
            /* b12 is rotated away with left rotation */
            if( fabs(b12)>fabs(b22) ) {
              cotphi=b22/b12;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b12/b22;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][1]; b=UO3[i][2];
              UO3[i][1]=a*cosphi-b*sinphi;
              UO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b22=b12*sinphi+b22*cosphi;
            b12=0.;

            nflops += 27;
          }
          else if( b22==0 ) {
#ifdef QPHIX_SVD3x3_DEBUG
            printf( "Entering case b22==0\n" );
#endif /* QPHIX_SVD3x3_DEBUG */
            /* b12 is rotated away and b02 appears,
               then b02 is rotated away, both are
               right rotations */
            if( fabs(b12)>fabs(b11) ) {
              cotphi=b11/b12;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b12/b11;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][1]; b=VO3[i][2];
              VO3[i][1]= a*cosphi+b*sinphi;
              VO3[i][2]=-a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b02=-b01*sinphi;
            b01=b01*cosphi;
            b11=b11*cosphi+b12*sinphi;
            b12=0.;
            /* second rotation removes b02 */
            if( fabs(b02)>fabs(b00) ) {
              cotphi=b00/b02;
              sinphi=1/sqrt(1+cotphi*cotphi);
              cosphi=cotphi*sinphi;
            }
            else {
              tanphi=b02/b00;
              cosphi=1/sqrt(1+tanphi*tanphi);
              sinphi=tanphi*cosphi;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][0]; b=VO3[i][2];
              VO3[i][0]= a*cosphi+b*sinphi;
              VO3[i][2]=-a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix */
            b00=b00*cosphi+b02*sinphi;
            b02=0.;

            nflops += 64;
          }
          else {
            /* full iteration with QR shift */
#ifdef QPHIX_SVD3x3_DEBUG
            printf( "Entering case of normal QR iteration\n" );
#endif /* QPHIX_SVD3x3_DEBUG */

            /* find max eigenvalue of bottom 2x2 minor */
            m11=b11*b11+b01*b01;
            m22=b22*b22+b12*b12;
            m12=b11*b12;
            dm=(m11-m22)/2;

            /* safely calculate sqrt */
            c=fabs(dm); d=fabs(m12);
            if( c>d ) {
              max=c; min=d;
            }
            else {
              max=d; min=c;
            }
            if( min==0 ) {
              norm = max;
            }
            else {
              c = min/max;
              norm = max*sqrt(1+c*c);
            }

            if( dm>=0 ) {
              lambdamax=m22-(m12*m12)/(dm+norm);
            }
            else {
              lambdamax=m22+(m12*m12)/(norm-dm);
            }

            /* calculate first Givens rotation (on the right) */
            a=b00*b00-lambdamax;
            b=b00*b01;
            if( 0==b ) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b)>fabs(a) ) {
                cotphi=-a/b;
                sinphi=1./sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b/a;
                cosphi=1./sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }
              nflops += 7;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][0]; b=VO3[i][1];
              VO3[i][0]=a*cosphi-b*sinphi;
              VO3[i][1]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generate b10 */
            a=b00; b=b01;
            b00=a*cosphi-b*sinphi;
            b01=a*sinphi+b*cosphi;
            b10=-b11*sinphi;
            b11=b11*cosphi; 

            /* calculate second Givens rotation (on the left) */
            if(0==b10) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b10)>fabs(b00) ) {
                cotphi=-b00/b10;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b10/b00;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

              nflops += 7;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][0]; b=UO3[i][1];
              UO3[i][0]= a*cosphi-b*sinphi;
              UO3[i][1]= a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b02 */
            b00=b00*cosphi-b10*sinphi;
            a=b01; b=b11;
            b01=a*cosphi-b*sinphi;
            b11=a*sinphi+b*cosphi;
            b02=-b12*sinphi;
            b12=b12*cosphi;
            b10=0.;

            /* calculate third Givens rotation (on the right) */
            if(0==b02) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b02)>fabs(b01) ) {
                cotphi=-b01/b02;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b02/b01;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

              nflops += 7;
            }
            /* multiply right VO3 matrix */
            for(i=0;i<3;i++) {
              a=VO3[i][1]; b=VO3[i][2];
              VO3[i][1]=a*cosphi-b*sinphi;
              VO3[i][2]=a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this generates b21 */
            b01=b01*cosphi-b02*sinphi;
            a=b11; b=b12;
            b11=a*cosphi-b*sinphi;
            b12=a*sinphi+b*cosphi;
            b21=-b22*sinphi;
            b22=b22*cosphi;
            b02=0.;

            /* calculate fourth Givens rotation (on the left) */
            if(0==b21) {
              cosphi=1.;
              sinphi=0.;
            }
            else {
              if( fabs(b21)>fabs(b11) ) {
                cotphi=-b11/b21;
                sinphi=1/sqrt(1+cotphi*cotphi);
                cosphi=cotphi*sinphi;
              }
              else {
                tanphi=-b21/b11;
                cosphi=1/sqrt(1+tanphi*tanphi);
                sinphi=tanphi*cosphi;
              }

              nflops += 7;
            }
            /* multiply left UO3 matrix */
            for(i=0;i<3;i++) {
              a=UO3[i][1]; b=UO3[i][2];
              UO3[i][1]= a*cosphi-b*sinphi;
              UO3[i][2]= a*sinphi+b*cosphi;
            }
            /* apply to bidiagonal matrix, this eliminates b21 */
            b11=b11*cosphi-b21*sinphi;
            a=b12; b=b22;
            b12=a*cosphi-b*sinphi;
            b22=a*sinphi+b*cosphi;
            b21=0.;

            nflops += 127;
          }
        } /* end of normal 3x3 iteration */
      }
    }
#ifdef QPHIX_SVD3x3_DEBUG
    printf( "QR iteration: %d\n", iter );
    printf( "%+20.16e %+20.16e %+20.16e\n", b00, b01, b02 );
    printf( "%+20.16e %+20.16e %+20.16e\n", b10, b11, b12 );
    printf( "%+20.16e %+20.16e %+20.16e\n", b20, b21, b22 );
#endif /* QPHIX_SVD3x3_DEBUG */
  }
  while( b01!=0 || b12!=0 );


  /* make singular values positive */
  if(b00<0) {
    b00=-b00;
    VO3[0][0]=-VO3[0][0];
    VO3[1][0]=-VO3[1][0];
    VO3[2][0]=-VO3[2][0];
  }
  if(b11<0) {
    b11=-b11;
    VO3[0][1]=-VO3[0][1];
    VO3[1][1]=-VO3[1][1];
    VO3[2][1]=-VO3[2][1];
  }
  if(b22<0) {
    b22=-b22;
    VO3[0][2]=-VO3[0][2];
    VO3[1][2]=-VO3[1][2];
    VO3[2][2]=-VO3[2][2];
  }



  /* Q=U1*U2 (U2 is block diagonal with U2_00=1) */
  Q[0][0][0]=U1[0][0][0]; Q[0][0][1]=U1[0][0][1];
  Q[1][0][0]=U1[1][0][0]; Q[1][0][1]=U1[1][0][1];
  Q[2][0][0]=U1[2][0][0]; Q[2][0][1]=U1[2][0][1];
  /* Q_01=U1_01*U2_11+U1_02*U2_21 */
  Q[0][1][0]=U1[0][1][0]*U2[1][1][0]-U1[0][1][1]*U2[1][1][1]
    +U1[0][2][0]*U2[2][1][0]-U1[0][2][1]*U2[2][1][1];
  Q[0][1][1]=U1[0][1][0]*U2[1][1][1]+U1[0][1][1]*U2[1][1][0]
    +U1[0][2][0]*U2[2][1][1]+U1[0][2][1]*U2[2][1][0];
  /* Q_02=U1_01*U2_12+U1_02*U2_22 */
  Q[0][2][0]=U1[0][1][0]*U2[1][2][0]-U1[0][1][1]*U2[1][2][1]
    +U1[0][2][0]*U2[2][2][0]-U1[0][2][1]*U2[2][2][1];
  Q[0][2][1]=U1[0][1][0]*U2[1][2][1]+U1[0][1][1]*U2[1][2][0]
    +U1[0][2][0]*U2[2][2][1]+U1[0][2][1]*U2[2][2][0];
  /* Q_11=U1_11*U2_11+U1_12*U2_21 */
  Q[1][1][0]=U1[1][1][0]*U2[1][1][0]-U1[1][1][1]*U2[1][1][1]
    +U1[1][2][0]*U2[2][1][0]-U1[1][2][1]*U2[2][1][1];
  Q[1][1][1]=U1[1][1][0]*U2[1][1][1]+U1[1][1][1]*U2[1][1][0]
    +U1[1][2][0]*U2[2][1][1]+U1[1][2][1]*U2[2][1][0];
  /* Q_12=U1_11*U2_12+U1_12*U2_22 */
  Q[1][2][0]=U1[1][1][0]*U2[1][2][0]-U1[1][1][1]*U2[1][2][1]
    +U1[1][2][0]*U2[2][2][0]-U1[1][2][1]*U2[2][2][1];
  Q[1][2][1]=U1[1][1][0]*U2[1][2][1]+U1[1][1][1]*U2[1][2][0]
    +U1[1][2][0]*U2[2][2][1]+U1[1][2][1]*U2[2][2][0];
  /* Q_21=U1_21*U2_11+U1_22*U2_21 */
  Q[2][1][0]=U1[2][1][0]*U2[1][1][0]-U1[2][1][1]*U2[1][1][1]
    +U1[2][2][0]*U2[2][1][0]-U1[2][2][1]*U2[2][1][1];
  Q[2][1][1]=U1[2][1][0]*U2[1][1][1]+U1[2][1][1]*U2[1][1][0]
    +U1[2][2][0]*U2[2][1][1]+U1[2][2][1]*U2[2][1][0];
  /* Q_22=U1_21*U2_12+U1_22*U2_22 */
  Q[2][2][0]=U1[2][1][0]*U2[1][2][0]-U1[2][1][1]*U2[1][2][1]
    +U1[2][2][0]*U2[2][2][0]-U1[2][2][1]*U2[2][2][1];
  Q[2][2][1]=U1[2][1][0]*U2[1][2][1]+U1[2][1][1]*U2[1][2][0]
    +U1[2][2][0]*U2[2][2][1]+U1[2][2][1]*U2[2][2][0];

  /* Q=Q*U3 (U3 is block diagonal with U3_00=1, U3_11=1)
     (this changes only third column of Q */
  a=Q[0][2][0]*U3[2][2][0]-Q[0][2][1]*U3[2][2][1];
  b=Q[0][2][0]*U3[2][2][1]+Q[0][2][1]*U3[2][2][0];
  Q[0][2][0]=a; Q[0][2][1]=b;
  a=Q[1][2][0]*U3[2][2][0]-Q[1][2][1]*U3[2][2][1];
  b=Q[1][2][0]*U3[2][2][1]+Q[1][2][1]*U3[2][2][0];
  Q[1][2][0]=a; Q[1][2][1]=b;
  a=Q[2][2][0]*U3[2][2][0]-Q[2][2][1]*U3[2][2][1];
  b=Q[2][2][0]*U3[2][2][1]+Q[2][2][1]*U3[2][2][0];
  Q[2][2][0]=a; Q[2][2][1]=b;

  nflops += 102;

  /* final U=Q*UO3
     (unitary times orthogonal that accumulated Givens rotations) */
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      QPHIX_Real tr,ti;
      tr = Q[i][0][0]*UO3[0][j]+Q[i][1][0]*UO3[1][j]+Q[i][2][0]*UO3[2][j];
      ti = Q[i][0][1]*UO3[0][j]+Q[i][1][1]*UO3[1][j]+Q[i][2][1]*UO3[2][j];
      QPHIX_c_eq_r_plus_ir(QPHIX_elem_M(*U,i,j), tr, ti);
    }
  }

  nflops += 90;

  /* Q=V1*V2 (V1 is block diagonal with V2_11=1,
     V2 is block diagonal with V2_11=1, V2_22=1) */
  Q[0][0][0]=V1[0][0][0]; Q[0][0][1]=V1[0][0][1];
  Q[1][0][0]=V1[1][0][0]; Q[1][0][1]=V1[1][0][1];
  Q[2][0][0]=V1[2][0][0]; Q[2][0][1]=V1[2][0][1];
  Q[0][1][0]=V1[0][1][0]; Q[0][1][1]=V1[0][1][1];
  Q[0][2][0]=V1[0][2][0]; Q[0][2][1]=V1[0][2][1];
  Q[1][1][0]=V1[1][1][0]; Q[1][1][1]=V1[1][1][1];
  Q[2][1][0]=V1[2][1][0]; Q[2][1][1]=V1[2][1][1];
  Q[1][2][0]=V1[1][2][0]*V2[2][2][0]-V1[1][2][1]*V2[2][2][1];
  Q[1][2][1]=V1[1][2][0]*V2[2][2][1]+V1[1][2][1]*V2[2][2][0];
  Q[2][2][0]=V1[2][2][0]*V2[2][2][0]-V1[2][2][1]*V2[2][2][1];
  Q[2][2][1]=V1[2][2][0]*V2[2][2][1]+V1[2][2][1]*V2[2][2][0];

  /* final V=Q*VO3
     (unitary times orthogonal that accumulated Givens rotations) */
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++) {
      QPHIX_Real tr,ti;
      tr = Q[i][0][0]*VO3[0][j]+Q[i][1][0]*VO3[1][j]+Q[i][2][0]*VO3[2][j];
      ti = Q[i][0][1]*VO3[0][j]+Q[i][1][1]*VO3[1][j]+Q[i][2][1]*VO3[2][j];
      QPHIX_c_eq_r_plus_ir(QPHIX_elem_M(*V,i,j), tr, ti);
    }
  }

  nflops += 102;

  /* singular values */
  sigma[0]=b00; sigma[1]=b11; sigma[2]=b22;

  info->final_flop += nflops;

  return 0;

#undef b00
#undef b01
#undef b02
#undef b10
#undef b11
#undef b12
#undef b20
#undef b21
#undef b22
}

/* SVD of 2x2 real matrix brought to the form:
   [ a00 a01]
   [   0 a11]
   This routine eliminates off-diagonal element, handling special cases */
static int QPHIX_svd2x2bidiag(QPHIX_info_t *info, QPHIX_Real *a00, QPHIX_Real *a01, 
    QPHIX_Real *a11, QPHIX_Real U2[2][2], QPHIX_Real V2[2][2]) {
#define b00 P[0][0][0]
#define b01 P[0][1][0]
#define b02 P[0][2][0]
#define b10 P[1][0][0]
#define b11 P[1][1][0]
#define b12 P[1][2][0]
#define b20 P[2][0][0]
#define b21 P[2][1][0]
#define b22 P[2][2][0]
  register double sinphi, cosphi, tanphi, cotphi;
  register double a, b, min, max, abs00, abs01, abs11;
  register double lna01a11, lna00, ln_num, tau, t;
  register double P00, P01, P10, P11;
  register int isign;
  size_t nflops = 0;

  U2[0][0]=1.; U2[0][1]=0.;
  U2[1][0]=0.; U2[1][1]=1.;
  V2[0][0]=1.; V2[0][1]=0.;
  V2[1][0]=0.; V2[1][1]=1.;

  if( *a00==0 ) {
    if( *a11==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(*a11)>fabs(*a01) ) {
        cotphi=-(*a01)/(*a11);
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-(*a11)/(*a01);
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 6;
    }
    /* multiply matrix A */
    (*a00)=cosphi*(*a01)-sinphi*(*a11);
    (*a01)=0.; (*a11)=0.;
    /* exchange columns in matrix V */
    V2[0][0]=0.; V2[0][1]=1.;
    V2[1][0]=1.; V2[1][1]=0.;
    /* U is just Givens rotation */
    U2[0][0]= cosphi; U2[0][1]= sinphi;
    U2[1][0]=-sinphi; U2[1][1]= cosphi;

    nflops += 3;
  }
  else if( *a11==0 ) {
    if( *a01==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(*a01)>fabs(*a00) ) {
        cotphi=-(*a00)/(*a01);
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-(*a01)/(*a00);
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 7;
    }
    /* multiply matrix A */
    (*a00)=cosphi*(*a00)-sinphi*(*a01);
    (*a01)=0.; (*a11)=0.;
    /* V is just Givens rotation */
    V2[0][0]= cosphi; V2[0][1]= sinphi;
    V2[1][0]=-sinphi; V2[1][1]= cosphi;
    nflops += 3;
  }
  else if( *a01==0 ){ /* nothing to be done */
    ;
  }
  else {
    /* need to calculate ( a11^2+a01^2-a00^2 )/( 2*a00*a01 )
       avoiding overflow/underflow,
       use logarithmic coding */
    abs01=fabs(*a01); abs11=fabs(*a11);
    if(abs01>abs11) {
      min=abs11; max=abs01;
    }
    else {
      min=abs01; max=abs11;
    }
    a=min/max;
    lna01a11=2*log(max)+log(1+a*a);

    abs00=fabs(*a00);
    lna00=2*log(abs00);
    if( lna01a11>lna00 ) {
      /* subtract smaller from larger, overall "+" */
      isign=1;
      ln_num=lna01a11+log(1.-exp(lna00-lna01a11));
    }
    else {
      /* subtract larger from smaller, need to change order, overall "-" */
      isign=-1;
      ln_num=lna00+log(1.-exp(lna01a11-lna00));
    }
    a=ln_num-log(2)-log(abs00)-log(abs01);
    tau=exp(a);
    tau*=isign;
    if(*a00<0.)
    {
      tau*=-1;
    }
    if(*a01<0.)
    {
      tau*=-1;
    }

    /* calculate b=sqrt(1+tau^2) */
    a=fabs(tau);
    if( a>1. ) {
      max=a; min=1.;
    }
    else {
      max=1.; min=a;
    }
    if( min==0 ) {
      b = max;
    }
    else {
      a = min/max;
      b = max*sqrt(1+a*a);
    }
    if(tau>=0.) {
      t = 1.0/(tau + b);
    }
    else {
      t = 1.0/(tau - b);
    }

    /* calculate b=sqrt(1+t^2) */
    a=fabs(t);
    if( a>1. ) {
      max=a; min=1.;
    }
    else {
      max=1.; min=a;
    }
    if( min==0 ) {
      b = max;
    }
    else {
      a = min/max;
      b = max*sqrt(1+a*a);
    }
    cosphi=1./b;
    sinphi=t*cosphi;

    /* transform matrix A so it has othogonal columns */
    P00= cosphi*(*a00)-sinphi*(*a01);
    P10=-sinphi*(*a11);
    P01= sinphi*(*a00)+cosphi*(*a01);
    P11= cosphi*(*a11);

    /* prepare V  */
    V2[0][0]= cosphi; V2[0][1]= sinphi;
    V2[1][0]=-sinphi; V2[1][1]= cosphi;

    /* make column with the largest norm first column */
    if( sqrt(P00*P00+P10*P10)<sqrt(P01*P01+P11*P11) ) {
      a=P00; P00=P01; P01=a;
      a=P10; P10=P11; P11=a;
      /* exchange columns in matrix V2 */
      a=V2[0][0]; V2[0][0]=V2[0][1]; V2[0][1]=a;
      a=V2[1][0]; V2[1][0]=V2[1][1]; V2[1][1]=a;
    }

    /* calculate left Givens rotation and diagonalize */
    if( P10==0 ) {
      cosphi=1.;
      sinphi=0.;
    }
    else {
      if( fabs(P10)>fabs(P00) ) {
        cotphi=-P00/P10;
        sinphi=1/sqrt(1+cotphi*cotphi);
        cosphi=cotphi*sinphi;
      }
      else {
        tanphi=-P10/P00;
        cosphi=1/sqrt(1+tanphi*tanphi);
        sinphi=tanphi*cosphi;
      }
      nflops += 7;
    }
    *a00=P00*cosphi-P10*sinphi;
    *a01=0.;
    *a11=P01*sinphi+P11*cosphi;

    /* U is just Givens rotation */
    U2[0][0]= cosphi; U2[0][1]= sinphi;
    U2[1][0]=-sinphi; U2[1][1]= cosphi;

    nflops += 56;
  }

  info->final_flop += nflops;

  return 0;
#undef b00
#undef b01
#undef b02
#undef b10
#undef b11
#undef b12
#undef b20
#undef b21
#undef b22
}


