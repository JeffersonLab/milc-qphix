#ifndef _GAUGE_FORCE_IMP_H
#define _GAUGE_FORCE_IMP_H

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
#include "gf_boundary.h"
#include "gauge_force_imp_complete_specialization.h"

/*
void 
gforce_imp_F( void *h_out, void *h_in, void *g_in, void *tmp, int cb, double epsilonH);

void
gforce_imp_D( void *h_out, void *h_in, void *g_in, void *tmp, int cb, double epsilonH);
*/
void *allocHermit_F();
void *allocHermit_D();
void *allocHermitHelper_F();
void *allocHermitHelper_D();
void *allocHermitHelperYZT_F();
void *allocHermitHelperYZT_D();
void update_h_u_gforce_imp_F( void *h_io, void *h_tmp, void *g_io, int cb, double epsilonH, double epsilonU, double kappaS, double kappaR, double kappaB );
void update_h_u_gforce_imp_D( void *h_io, void *h_tmp, void *g_io, int cb, double epsilonH, double epsilonU, double kappaS, double kappaR, double kappaB );

#if QPHIX_PrecisionInt==1
#define update_h_u_gforce_imp update_h_u_gforce_imp_F
#elif QPHIX_PrecisionInt==2
#define update_h_u_gforce_imp update_h_u_gforce_imp_D
#else
#error "QPHIX_PrecisionInt not defined/supported!"
#endif

#endif
