#ifdef ENABLE_MPI
#include <mpi.h>
#endif
#include "ks_blas_utils_c.h"
#include "ks_config.h"  /* macros */
#include "ks_globals.h" /* All global vars and macros (Nt, Nz, Nx, Ny etc.) */
#include <math.h>

/* Squaring the magnitude of a VECLEN array of spinors. */
double
magsq_su3vec (const KS *spinor)
{
    double source_norm[VECLEN] = {0};
    double norm = 0.0;
    for(int c = 0; c < 3; c++) {
        for(int x = 0; x < VECLEN; x++) {
            source_norm[x] += double((*spinor)[c][0][x] * (*spinor)[c][0][x] 
                                  + (*spinor)[c][1][x] * (*spinor)[c][1][x]
                                   );
        }
    }
    
    /* Reduction for soalen sized vector */
    for(int i = 0; i < VECLEN; i++) {
        norm += source_norm[i];
    }
    
    return norm;
}

/* axpy */
void 
scalar_mult_add_su3_vector (const KS *sp1, const KS *sp2, fptype s, KS *sp3)
{
    for(int c = 0; c < 3; c++) {
        for(int x = 0; x < VECLEN; x++) {
            (*sp3)[c][0][x] = (*sp1)[c][0][x] + s*(*sp2)[c][0][x];
            (*sp3)[c][1][x] = (*sp1)[c][1][x] + s*(*sp2)[c][1][x];
        }
    }  
}

void
scalar_mult_su3_vector (const KS *in, double scalar, KS *out)
{
    for(int c = 0; c < 3; c++) {
        for(int x = 0; x < VECLEN; x++) {
            (*out)[c][0][x] = (*in)[c][0][x] * scalar;
            (*out)[c][1][x] = (*in)[c][1][x] * scalar;
        }
    }  
}

/* Dot product of two KS* spinors */
double
su3_rdot (const KS *sp1, const KS *sp2)
{
    double ret[VECLEN] = {0};
    double rdot = 0;

    for(int c = 0; c < 3; c++) {
        for(int x = 0; x < VECLEN; x++) {
            ret[x] += double((*sp1)[c][0][x] * (*sp2)[c][0][x]);
            ret[x] += double((*sp1)[c][1][x] * (*sp2)[c][1][x]);
        }
    }      

    for(int i = 0; i < VECLEN; ++i)
        rdot+= ret[i];
    
    return rdot;
}

/* Add two KS* spinors */
void 
add_su3_vector (const KS *sp1, const KS *sp2, KS *sp3) {
    for(int c = 0; c < 3; c++) {
        for(int x = 0; x < VECLEN; x++) {
            (*sp3)[c][0][x] = (*sp1)[c][0][x] + (*sp2)[c][0][x];
            (*sp3)[c][1][x] = (*sp1)[c][1][x] + (*sp2)[c][1][x];
        }
    }    
}

/* Subtract two KS* spinors */
void 
sub_su3_vector (const KS *sp1, const KS *sp2, KS *sp3) {
    for(int c = 0; c < 3; c++) {
        for(int x = 0; x < VECLEN; x++) {
            (*sp3)[c][0][x] = (*sp1)[c][0][x] - (*sp2)[c][0][x];
            (*sp3)[c][1][x] = (*sp1)[c][1][x] - (*sp2)[c][1][x];
        }
    }    
}

/*! \fixme global allreduce is not needed if we are doing single threaded 
 * execution. 
 */
#if QPHIX_PrecisionInt==1 /* Avoid double definition */
void
g_doublesum (double *d)
{
#ifdef ENABLE_MPI
    double work = *d;
    MPI_Allreduce( &work, d, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_THISJOB );
#endif
}
#endif

static double
calc_soalen_rel_residue (KS *sp1, KS *sp2)
{
    double residue = 0.0;
    double num[VECLEN] = {0}, den[VECLEN] = {0};

    for(int c = 0; c < 3; c++) {
        for(int x = 0; x < VECLEN; x++) {
            num[x] += double((*sp1)[c][0][x] * (*sp1)[c][0][x] 
                          + (*sp1)[c][1][x] * (*sp1)[c][1][x]);
            den[x] += double((*sp2)[c][0][x] * (*sp2)[c][0][x] 
                          + (*sp2)[c][1][x] * (*sp2)[c][1][x]);
        }
    }
    
    /* get the residue */
    for(int i = 0; i < VECLEN; i++) {
        //printf("DEBUG: MBENCH num : %g den : %g\n", num[i], den[i]);
	if(num[i] == 0)continue;
        residue += (den[i] == 0) ? 1.0 : (num[i]/den[i]);
        //printf("DEBUG: MBENCH residue : %g\n", residue);
    }    
    
    return residue;
}

fptype 
relative_residue (KS *p, KS *q, int parity)
{
    double residue = 0.0;
#pragma omp parallel for num_threads (nThreads) reduction(+:residue)
    // for(int i = 0; i < Vxh*Vy*Vz*Vt; i++) {
    for(int i = 0; i < Pxyz*Vt; i++) {
        residue += calc_soalen_rel_residue(&p[i], &q[i]);
    }
    //printf("MBENCH Residue : %g\n", residue);
    g_doublesum(&residue);
    //printf("Ret Residue : %g\n", sqrt(2*residue/(Nx*Ny*Nz*Nt)));
    return sqrt(2*residue/(Gx*Gy*Gz*Gt));
}

/***************************** Unit test functions ***********************/
#include "ks_unit_tests.h"

#if 0
double
t_magsq_su3vec (void *ks_spinor)
{
    KS * spinor = (KS*) ks_spinor;
    return magsq_su3vec(spinor);
}

void 
t_scalar_mult_add_su3_vector (const void *sp1, const void *sp2 
                            , double s, void *sp3)
{
    const KS *spinor1 = (const KS*)sp1;
    const KS *spinor2 = (const KS*)sp2;
    KS *spinor3       = (KS*)sp3;
    scalar_mult_add_su3_vector(spinor1, spinor2, s, spinor3);
}

void
t_scalar_mult_su3_vector (const void *sp1, double s, void *sp2)
{
    const KS *spinor1 = (const KS*)sp1;
    KS *spinor2       = (KS*)sp2;
    scalar_mult_su3_vector (spinor1, s, spinor2);
}

double
t_su3_rdot (void *sp1, void *sp2)
{
    KS *spinor1 = (KS*)sp1;
    KS *spinor2 = (KS*)sp2;
    
    return su3_rdot(spinor1, spinor2);
}

void 
t_add_su3_vector (const void *sp1, const void *sp2, void *sp3)
{
    const KS* spinor1 = (const KS*)sp1;
    const KS* spinor2 = (const KS*)sp2;
    KS* spinor3 = (KS*)sp3;
    add_su3_vector(spinor1, spinor2, spinor3);
}
 
//void
//t_g_doublesum (double *d)
//{
//}

double 
t_relative_residue (void *p, void *q, int parity)
{
    KS* sp1 = (KS*)p;
    KS* sp2 = (KS*)q;
    
    return relative_residue(sp1, sp2, parity);
}

void 
t_sub_su3_vector (const void *sp1, const void *sp2, void *sp3)
{
    const KS* spinor1 = (const KS*)sp1;
    const KS* spinor2 = (const KS*)sp2;
    KS* spinor3 = (KS*)sp3;
    sub_su3_vector(spinor1, spinor2, spinor3);
}

#endif
