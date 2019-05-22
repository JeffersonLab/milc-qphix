#include <omp.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "qphix_internal.h"
#include "fermion_force.h"
#include "ks_globals.h"
#include "layout.h"
#include "qphix_su3_algebra.h"

extern int qphix_sites_on_node;
extern int qphix_even_sites_on_node;
extern int qphix_fused_sites_on_node;

extern "C"
{
/* M1 =- M2 */
/* QPHIX FF algorithm needs this operation on only odd sites
 * hence, writing it for odd sites alone right now. fix it later
 * */
void QPHIX_M_eqm_M(QPHIX_ColorMatrix *m1, QPHIX_ColorMatrix *m2, QPHIX_evenodd_t parity)
{
  int s, i, j, v;
  QSU3M *m1e, *m1o, *m2e, *m2o;

  m1o = (QSU3M *)m1->odd;
  m2o = (QSU3M *)m2->odd;

  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<qphix_fused_sites_on_node; s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1o[s][i][j][0][v] =- m2o[s][i][j][0][v];
          m1o[s][i][j][1][v] =- m2o[s][i][j][1][v];
        }
      }
    }
  }
}

/* M1 = V1 * adjoint(V2) */
  void QPHIX_M_eq_V_times_Va(QPHIX_ColorMatrix *m1, QPHIX_ColorVector *v1, QPHIX_ColorVector *v2)
  {
    register int s, i, j, v;
    QSU3V *v1e, *v1o, *v2e, *v2o;
    QSU3M *m1e, *m1o;

    v1e = (QSU3V *)v1->even;
    v1o = (QSU3V *)v1->odd;
    v2e = (QSU3V *)v2->even;
    v2o = (QSU3V *)v2->odd;
    m1e = (QSU3M *)m1->even;
    m1o = (QSU3M *)m1->odd;

    __assume_aligned(v1e, VECLENBYTES);
    __assume_aligned(v1o, VECLENBYTES);
    __assume_aligned(v2e, VECLENBYTES);
    __assume_aligned(v2o, VECLENBYTES);
    __assume_aligned(m1e, VECLENBYTES);
    __assume_aligned(m1o, VECLENBYTES);

#pragma omp parallel for private(i, j, v)
    for(s=0; s<qphix_fused_sites_on_node; s++)
    {
      for(i=0; i<3; i++)
      {
        for(j=0; j<3; j++)
        {
#pragma simd
          for(v=0; v<VECLEN; v++)
          {
            m1e[s][i][j][0][v] = v1e[s][i][0][v] * v2e[s][j][0][v] + v1e[s][i][1][v] * v2e[s][j][1][v];
            m1e[s][i][j][1][v] = v1e[s][i][1][v] * v2e[s][j][0][v] - v1e[s][i][0][v] * v2e[s][j][1][v];
            m1o[s][i][j][0][v] = v1o[s][i][0][v] * v2o[s][j][0][v] + v1o[s][i][1][v] * v2o[s][j][1][v];
            m1o[s][i][j][1][v] = v1o[s][i][1][v] * v2o[s][j][0][v] - v1o[s][i][0][v] * v2o[s][j][1][v];
          }
        }
      }
    }
  }

  /* make_anti_hermitian */
  void QPHIX_M_eq_antiherm_M( QPHIX_ColorMatrix *m2, QPHIX_ColorMatrix *m1 ) {
    QPHIX_Real temp;
    register int s, v;
    QSU3M *m1e, *m1o, *m2e, *m2o;

    m1e = (QSU3M *)m1->even;
    m1o = (QSU3M *)m1->odd;
    m2e = (QSU3M *)m2->even;
    m2o = (QSU3M *)m2->odd;

    __assume_aligned(m1e, VECLENBYTES);
    __assume_aligned(m1o, VECLENBYTES);
    __assume_aligned(m2e, VECLENBYTES);
    __assume_aligned(m2o, VECLENBYTES);

/*------------------
	mat_su3->e[0][0].imag=mat_antihermit->m00im;
	mat_su3->e[0][0].real=0.;
	mat_su3->e[1][1].imag=mat_antihermit->m11im;
	mat_su3->e[1][1].real=0.;
	mat_su3->e[2][2].imag=mat_antihermit->m22im;
	mat_su3->e[2][2].real=0.;
	mat_su3->e[0][1].imag=mat_antihermit->m01.imag;
	temp1=mat_antihermit->m01.real;
	mat_su3->e[0][1].real=temp1;
	mat_su3->e[1][0].real= -temp1;
	mat_su3->e[1][0].imag=mat_antihermit->m01.imag;
	mat_su3->e[0][2].imag=mat_antihermit->m02.imag;
	temp1=mat_antihermit->m02.real;
	mat_su3->e[0][2].real=temp1;
	mat_su3->e[2][0].real= -temp1;
	mat_su3->e[2][0].imag=mat_antihermit->m02.imag;
	mat_su3->e[1][2].imag=mat_antihermit->m12.imag;
	temp1=mat_antihermit->m12.real;
	mat_su3->e[1][2].real=temp1;
	mat_su3->e[2][1].real= -temp1;
	mat_su3->e[2][1].imag=mat_antihermit->m12.imag;
-------------------------*/
#pragma omp parallel for private(v)
    for(s=0; s<qphix_fused_sites_on_node; s++)
    {
//#pragma omp simd
      for(v=0; v<VECLEN; v++)
      {
        temp = \
          (m1e[s][0][0][1][v] + m1e[s][1][1][1][v] + m1e[s][2][2][1][v])*0.33333333333333333;
        m2e[s][0][0][0][v] = 0.0;
        m2e[s][0][0][1][v] = m1e[s][0][0][1][v] - temp;
        m2e[s][1][1][0][v] = 0.0;
        m2e[s][1][1][1][v] = m1e[s][1][1][1][v] - temp;
        m2e[s][2][2][0][v] = 0.0;
        m2e[s][2][2][1][v] = m1e[s][2][2][1][v] - temp;

/*        m2e[s][0][0][0][v] = m1e[s][0][0][0][v];
        m2e[s][0][0][1][v] = m1e[s][0][0][1][v];
        m2e[s][1][1][0][v] = m1e[s][1][1][0][v];
        m2e[s][1][1][1][v] = m1e[s][1][1][1][v];
        m2e[s][2][2][0][v] = m1e[s][2][2][0][v];
        m2e[s][2][2][1][v] = m1e[s][2][2][1][v];*/

        m2e[s][0][1][0][v] = (m1e[s][0][1][0][v] - m1e[s][1][0][0][v])*0.5;
        m2e[s][0][1][1][v] = (m1e[s][0][1][1][v] + m1e[s][1][0][1][v])*0.5;
        m2e[s][1][0][0][v] = -1.0*m2e[s][0][1][0][v];
        m2e[s][1][0][1][v] =      m2e[s][0][1][1][v];

        m2e[s][0][2][0][v] = (m1e[s][0][2][0][v] - m1e[s][2][0][0][v])*0.5;
        m2e[s][0][2][1][v] = (m1e[s][0][2][1][v] + m1e[s][2][0][1][v])*0.5;
        m2e[s][2][0][0][v] = -1.0*m2e[s][0][2][0][v];
	m2e[s][2][0][1][v] =      m2e[s][0][2][1][v];

        m2e[s][1][2][0][v] = (m1e[s][1][2][0][v] - m1e[s][2][1][0][v])*0.5;
        m2e[s][2][1][0][v] = -1.0*m2e[s][1][2][0][v];
        m2e[s][1][2][1][v] = (m1e[s][1][2][1][v] + m1e[s][2][1][1][v])*0.5;
        m2e[s][2][1][1][v] =      m2e[s][1][2][1][v];

        temp = \
          (m1o[s][0][0][1][v] + m1o[s][1][1][1][v] + m1o[s][2][2][1][v])*0.33333333333333333;
        m2o[s][0][0][0][v] = 0.0;
        m2o[s][0][0][1][v] = m1o[s][0][0][1][v] - temp;
        m2o[s][1][1][0][v] = 0.0;
        m2o[s][1][1][1][v] = m1o[s][1][1][1][v] - temp;
        m2o[s][2][2][0][v] = 0.0;
        m2o[s][2][2][1][v] = m1o[s][2][2][1][v] - temp;

/*        m2o[s][0][0][0][v] = m1o[s][0][0][0][v];
        m2o[s][0][0][1][v] = m1o[s][0][0][1][v];
        m2o[s][1][1][0][v] = m1o[s][1][1][0][v];
        m2o[s][1][1][1][v] = m1o[s][1][1][1][v];
        m2o[s][2][2][0][v] = m1o[s][2][2][0][v];
        m2o[s][2][2][1][v] = m1o[s][2][2][1][v];*/

        m2o[s][0][1][0][v] = (m1o[s][0][1][0][v] - m1o[s][1][0][0][v])*0.5;
        m2o[s][0][1][1][v] = (m1o[s][0][1][1][v] + m1o[s][1][0][1][v])*0.5;
        m2o[s][1][0][0][v] = -1.0*m2o[s][0][1][0][v];
        m2o[s][1][0][1][v] =      m2o[s][0][1][1][v];

        m2o[s][0][2][0][v] = (m1o[s][0][2][0][v] - m1o[s][2][0][0][v])*0.5;
        m2o[s][0][2][1][v] = (m1o[s][0][2][1][v] + m1o[s][2][0][1][v])*0.5;
        m2o[s][2][0][0][v] = -1.0*m2o[s][0][2][0][v];
	m2o[s][2][0][1][v] =      m2o[s][0][2][1][v];

        m2o[s][1][2][0][v] = (m1o[s][1][2][0][v] - m1o[s][2][1][0][v])*0.5;
        m2o[s][2][1][0][v] = -1.0*m2o[s][1][2][0][v];
        m2o[s][1][2][1][v] = (m1o[s][1][2][1][v] + m1o[s][2][1][1][v])*0.5;
        m2o[s][2][1][1][v] =      m2o[s][1][2][1][v];

      }
    }
  }

}

/* M1 = 0.0 */
void QPHIX_M_eq_zero(QPHIX_ColorMatrix *m1)
{
  int s, i, j, v;
  QSU3M *m1e, *m1o;

  m1e = (QSU3M *)m1->even;
  m1o = (QSU3M *)m1->odd;

  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<qphix_fused_sites_on_node; s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1e[s][i][j][0][v] = 0.0;
          m1e[s][i][j][1][v] = 0.0;
          m1o[s][i][j][0][v] = 0.0;
          m1o[s][i][j][1][v] = 0.0;
        }
      }
    }
  }
}

/* M1 = M2 */
void QPHIX_M_eq_M(QPHIX_ColorMatrix *m1, QPHIX_ColorMatrix *m2)
{
  int s, i, j, v;
  QSU3M *m1e, *m1o, *m2e, *m2o;

  m1e = (QSU3M *)m1->even;
  m1o = (QSU3M *)m1->odd;
  m2e = (QSU3M *)m2->even;
  m2o = (QSU3M *)m2->odd;

  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<qphix_fused_sites_on_node; s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1e[s][i][j][0][v] = m2e[s][i][j][0][v];
          m1e[s][i][j][1][v] = m2e[s][i][j][1][v];
          m1o[s][i][j][0][v] = m2o[s][i][j][0][v];
          m1o[s][i][j][1][v] = m2o[s][i][j][1][v];
        }
      }
    }
  }
}


/* M1 = M1 + M2 */
void QPHIX_M_peq_M(QPHIX_ColorMatrix *m1, QPHIX_ColorMatrix *m2)
{
  int s, i, j, v;
  QSU3M *m1e, *m1o, *m2e, *m2o;

  m1e = (QSU3M *)m1->even;
  m1o = (QSU3M *)m1->odd;
  m2e = (QSU3M *)m2->even;
  m2o = (QSU3M *)m2->odd;

  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<qphix_fused_sites_on_node; s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1e[s][i][j][0][v] += m2e[s][i][j][0][v];
          m1e[s][i][j][1][v] += m2e[s][i][j][1][v];
          m1o[s][i][j][0][v] += m2o[s][i][j][0][v];
          m1o[s][i][j][1][v] += m2o[s][i][j][1][v];
        }
      }
    }
  }
}

/* M1 = r * M2 */
void QPHIX_M_eq_r_times_M(QPHIX_ColorMatrix *m1, QPHIX_Real a, QPHIX_ColorMatrix *m2)
{
  int s, i, j, v;
  QSU3M *m1e, *m1o, *m2e, *m2o;

  m1e = (QSU3M *)m1->even;
  m1o = (QSU3M *)m1->odd;
  m2e = (QSU3M *)m2->even;
  m2o = (QSU3M *)m2->odd;

  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<qphix_fused_sites_on_node; s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1e[s][i][j][0][v] = a * m2e[s][i][j][0][v];
          m1e[s][i][j][1][v] = a * m2e[s][i][j][1][v];
          m1o[s][i][j][0][v] = a * m2o[s][i][j][0][v];
          m1o[s][i][j][1][v] = a * m2o[s][i][j][1][v];
        }
      }
    }
  }
}

/* M1 = M1 + (r * M2) */
void QPHIX_M_peq_r_times_M(QPHIX_ColorMatrix *m1, QPHIX_Real a, QPHIX_ColorMatrix *m2)
{
  int s, i, j, v;
  QSU3M *m1e, *m1o, *m2e, *m2o;

  m1e = (QSU3M *)m1->even;
  m1o = (QSU3M *)m1->odd;
  m2e = (QSU3M *)m2->even;
  m2o = (QSU3M *)m2->odd;

  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<qphix_fused_sites_on_node; s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1e[s][i][j][0][v] += a * m2e[s][i][j][0][v];
          m1e[s][i][j][1][v] += a * m2e[s][i][j][1][v];
          m1o[s][i][j][0][v] += a * m2o[s][i][j][0][v];
          m1o[s][i][j][1][v] += a * m2o[s][i][j][1][v];
        }
      }
    }
  }
}

/* M1 = M1 + (M2 x M3) */
void QPHIX_M_peq_M_times_M(QPHIX_ColorMatrix *m1, QPHIX_ColorMatrix *m2, QPHIX_ColorMatrix *m3)
{
  int s, i, j, v, k;
  QSU3M *m1e, *m1o, *m2e, *m2o, *m3e, *m3o;

  m1e = (QSU3M *)m1->even;
  m1o = (QSU3M *)m1->odd;
  m2e = (QSU3M *)m2->even;
  m2o = (QSU3M *)m2->odd;
  m3o = (QSU3M *)m3->odd;
  m3e = (QSU3M *)m3->even;

  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);
  __assume_aligned(m3e, VECLENBYTES);
  __assume_aligned(m3o, VECLENBYTES);

#pragma omp parallel for private(i, j, v, k)
  for(s=0; s<qphix_fused_sites_on_node; s++)
  {
    for(i=0; i<3; i++)
    {
      for(k=0; k<3; k++)
      {
        for(j=0; j<3; j++)
        {
#pragma omp simd
          for(v=0; v<VECLEN; v++)
          {
            m1e[s][i][j][0][v] += (m2e[s][i][k][0][v] * m3e[s][k][j][0][v]) -
              (m2e[s][i][k][1][v] * m3e[s][k][j][1][v]);
            m1e[s][i][j][1][v] += (m2e[s][i][k][0][v] * m3e[s][k][j][1][v]) +
              (m2e[s][i][k][1][v] * m3e[s][k][j][0][v]);
            m1o[s][i][j][0][v] += (m2o[s][i][k][0][v] * m3o[s][k][j][0][v]) -
              (m2o[s][i][k][1][v] * m3o[s][k][j][1][v]);
            m1o[s][i][j][1][v] += (m2o[s][i][k][0][v] * m3o[s][k][j][1][v]) +
              (m2o[s][i][k][1][v] * m3o[s][k][j][0][v]);
          }
        }
      }
    }
  }
}

/* M1 = M2 x M3 */
void QPHIX_M_eq_M_times_M(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2)
{
  int s, i, j, v,k;
  QSU3M *m1e, *m1o, *m2e, *m2o, *m3e, *m3o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  m3e = (QSU3M *)M3->even;
  m3o = (QSU3M *)M3->odd;
  QPHIX_M_eq_zero(M3);
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);
  __assume_aligned(m3e, VECLENBYTES);
  __assume_aligned(m3o, VECLENBYTES);

#pragma omp parallel for private(i,j,v,k)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(k=0; k<3; k++)
      {
        for(j=0; j<3; j++)
        {
#pragma omp simd
          for(v=0; v<VECLEN; v++)
          {
            m3e[s][i][j][0][v] +=  m1e[s][i][k][0][v] * m2e[s][k][j][0][v] - m1e[s][i][k][1][v] * m2e[s][k][j][1][v];
            m3e[s][i][j][1][v] +=  m1e[s][i][k][0][v] * m2e[s][k][j][1][v] + m1e[s][i][k][1][v] * m2e[s][k][j][0][v];
            m3o[s][i][j][0][v] +=  m1o[s][i][k][0][v] * m2o[s][k][j][0][v] - m1o[s][i][k][1][v] * m2o[s][k][j][1][v];
            m3o[s][i][j][1][v] +=  m1o[s][i][k][0][v] * m2o[s][k][j][1][v] + m1o[s][i][k][1][v] * m2o[s][k][j][0][v];

          }
        }
      }
    }
  }
}


/* M1 = M1 + (M2 x M3a) */
void QPHIX_M_peq_M_times_Ma(QPHIX_ColorMatrix *m1, QPHIX_ColorMatrix *m2, QPHIX_ColorMatrix *m3)
{
  int s, i, j, v, k;
  QSU3M *m1e, *m1o, *m2e, *m2o, *m3e, *m3o;

  m1e = (QSU3M *)m1->even;
  m1o = (QSU3M *)m1->odd;
  m2e = (QSU3M *)m2->even;
  m2o = (QSU3M *)m2->odd;
  m3o = (QSU3M *)m3->odd;
  m3e = (QSU3M *)m3->even;

  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);
  __assume_aligned(m3e, VECLENBYTES);
  __assume_aligned(m3o, VECLENBYTES);

#pragma omp parallel for private(i, j, v, k)
  for(s=0; s<qphix_fused_sites_on_node; s++)
  {
    for(i=0; i<3; i++)
    {
      for(k=0; k<3; k++)
      {
        for(j=0; j<3; j++)
        {
#pragma omp simd
          for(v=0; v<VECLEN; v++)
          {
            m1e[s][i][j][0][v] += (m2e[s][i][k][0][v] * m3e[s][j][k][0][v]) +
              (m2e[s][i][k][1][v] * m3e[s][j][k][1][v]);
            m1e[s][i][j][1][v] += (m2e[s][i][k][1][v] * m3e[s][j][k][0][v]) -
              (m2e[s][i][k][0][v] * m3e[s][j][k][1][v]);
            m1o[s][i][j][0][v] += (m2o[s][i][k][0][v] * m3o[s][j][k][0][v]) +
              (m2o[s][i][k][1][v] * m3o[s][j][k][1][v]);
            m1o[s][i][j][1][v] += (m2o[s][i][k][1][v] * m3o[s][j][k][0][v]) -
              (m2o[s][i][k][0][v] * m3o[s][j][k][1][v]);
          }
        }
      }
    }
  }
}

/* M1 = M2 x adjoint(M3) */
void QPHIX_M_eq_M_times_Ma(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2)
{
  int s, i, j, v,k;
  QSU3M *m1e, *m1o, *m2e, *m2o, *m3e, *m3o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  m3e = (QSU3M *)M3->even;
  m3o = (QSU3M *)M3->odd;
  QPHIX_M_eq_zero(M3);
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);
  __assume_aligned(m3e, VECLENBYTES);
  __assume_aligned(m3o, VECLENBYTES);

#pragma omp parallel for private(i,j,v,k)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(k=0; k<3; k++)
      {
        for(j=0; j<3; j++)
        {
#pragma omp simd
          for(v=0; v<VECLEN; v++)
          {
            m3e[s][i][j][0][v] +=  m1e[s][i][k][0][v] * m2e[s][j][k][0][v] + m1e[s][i][k][1][v] * m2e[s][j][k][1][v];
            m3e[s][i][j][1][v] +=  m1e[s][i][k][1][v] * m2e[s][j][k][0][v] - m1e[s][i][k][0][v] * m2e[s][j][k][1][v];
            m3o[s][i][j][0][v] +=  m1o[s][i][k][0][v] * m2o[s][j][k][0][v] + m1o[s][i][k][1][v] * m2o[s][j][k][1][v];
            m3o[s][i][j][1][v] +=  m1o[s][i][k][1][v] * m2o[s][j][k][0][v] - m1o[s][i][k][0][v] * m2o[s][j][k][1][v];

          }
        }
      }
    }
  }
}


/* M1 += adjoint(M2) */
void QPHIX_M_peq_Ma(QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2)
{
  int s, i, j, v;
  QSU3M *m1e, *m1o, *m2e, *m2o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1e[s][i][j][0][v] += m2e[s][j][i][0][v];
          m1e[s][i][j][1][v] += -m2e[s][j][i][1][v];
          m1o[s][i][j][0][v] += m2o[s][j][i][0][v];
          m1o[s][i][j][1][v] += -m2o[s][j][i][1][v];
        }
      }
    }
  }
}

/* M1 = adjoint(M2) */
void QPHIX_M_eq_Ma(QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2)
{
  int s, i, j, v;
  QSU3M *m1e, *m1o, *m2e, *m2o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1e[s][i][j][0][v] = m2e[s][j][i][0][v];
          m1e[s][i][j][1][v] = -m2e[s][j][i][1][v];
          m1o[s][i][j][0][v] = m2o[s][j][i][0][v];
          m1o[s][i][j][1][v] = -m2o[s][j][i][1][v];
        }
      }
    }
  }
}
/* M3 = adjoint(M1) * M2 */
void QPHIX_M_eq_Ma_times_M(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2)
{
  int s, i, j, v,k;
  QSU3M *m1e, *m1o, *m2e, *m2o, *m3e, *m3o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  m3e = (QSU3M *)M3->even;
  m3o = (QSU3M *)M3->odd;
  QPHIX_M_eq_zero(M3);
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);
  __assume_aligned(m3e, VECLENBYTES);
  __assume_aligned(m3o, VECLENBYTES);

#pragma omp parallel for private(i,j,v,k)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(k=0; k<3; k++)
      {
        for(j=0; j<3; j++)
        {
#pragma omp simd
          for(v=0; v<VECLEN; v++)
          {
            m3e[s][i][j][0][v] +=  m1e[s][k][i][0][v] * m2e[s][k][j][0][v] + m1e[s][k][i][1][v] * m2e[s][k][j][1][v];
            m3e[s][i][j][1][v] +=  m1e[s][k][i][0][v] * m2e[s][k][j][1][v] - m1e[s][k][i][1][v] * m2e[s][k][j][0][v];
            m3o[s][i][j][0][v] +=  m1o[s][k][i][0][v] * m2o[s][k][j][0][v] + m1o[s][k][i][1][v] * m2o[s][k][j][1][v];
            m3o[s][i][j][1][v] +=  m1o[s][k][i][0][v] * m2o[s][k][j][1][v] - m1o[s][k][i][1][v] * m2o[s][k][j][0][v];

          }
        }
      }
    }
  }
}

/* M3 = r *M1 + M2 */
void QPHIX_M_eq_r_times_M_plus_M(QPHIX_ColorMatrix *M3, QPHIX_Real r, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2)
{
  register int i,j,k,s,v;
  QSU3M *m1e, *m1o, *m2e, *m2o, *m3e, *m3o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  m3e = (QSU3M *)M3->even;
  m3o = (QSU3M *)M3->odd;
  QPHIX_M_eq_zero(M3);
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);
  __assume_aligned(m3e, VECLENBYTES);
  __assume_aligned(m3o, VECLENBYTES);

#pragma omp parallel for private(i,j,v,k)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m3e[s][i][j][0][v] = r * m1e[s][i][j][0][v] + m2e[s][i][j][0][v];
          m3e[s][i][j][1][v] = r * m1e[s][i][j][1][v] + m2e[s][i][j][1][v];
          m3o[s][i][j][0][v] = r * m1o[s][i][j][0][v] + m2o[s][i][j][0][v];
          m3o[s][i][j][1][v] = r * m1o[s][i][j][1][v] + m2o[s][i][j][1][v];
        }
      }
    }
  }
}

/* M3 = adjoint(M1) * adjoint(M2) */
void QPHIX_M_eq_Ma_times_Ma(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2)
{
  int s, i, j, k, v;
  QSU3M *m1e, *m1o, *m2e, *m2o, *m3e, *m3o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  m3e = (QSU3M *)M3->even;
  m3o = (QSU3M *)M3->odd;
  QPHIX_M_eq_zero(M3);
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);
  __assume_aligned(m3e, VECLENBYTES);
  __assume_aligned(m3o, VECLENBYTES);

#pragma omp parallel for private(i,j,v,k)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(k=0; k<3; k++)
      {
        for(j=0; j<3; j++)
        {
#pragma omp simd
          for(v=0; v<VECLEN; v++)
          {
            m3e[s][i][j][0][v] +=  m1e[s][k][i][0][v] * m2e[s][j][k][0][v] - m1e[s][k][i][1][v] * m2e[s][j][k][1][v];
            m3e[s][i][j][1][v] +=  m1e[s][k][i][0][v] * -m2e[s][j][k][1][v] - m1e[s][k][i][1][v] * m2e[s][j][k][0][v];
            m3o[s][i][j][0][v] +=  m1o[s][k][i][0][v] * m2o[s][j][k][0][v] - m1o[s][k][i][1][v] * m2o[s][j][k][1][v];
            m3o[s][i][j][1][v] +=  m1o[s][k][i][0][v] * -m2o[s][j][k][1][v] - m1o[s][k][i][1][v] * m2o[s][j][k][0][v];

          }
        }
      }
    }
  }
}

/* M1 -= M2 */
void QPHIX_M_meq_M(QPHIX_ColorMatrix *M1, QPHIX_ColorMatrix *M2)
{
  register int i,j,s,v;
  QSU3M *m1e, *m1o, *m2e, *m2o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);

#pragma omp parallel for private(i,j,v)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(j=0; j<3; j++)
      {
#pragma omp simd
        for(v=0; v<VECLEN; v++)
        {
          m1e[s][i][j][0][v] -= m2e[s][i][j][0][v];
          m1e[s][i][j][1][v] -= m2e[s][i][j][1][v];
          m1o[s][i][j][0][v] -= m2o[s][i][j][0][v];
          m1o[s][i][j][1][v] -= m2o[s][i][j][1][v];
        }
      }
    }
  }
}

/* M3 = M3 - (M1 * M2) */ 
void QPHIX_M_meq_M_times_M(QPHIX_ColorMatrix *M3, QPHIX_ColorMatrix *M2, QPHIX_ColorMatrix *M1)
{
  int s, i, j, v,k;
  QSU3M *m1e, *m1o, *m2e, *m2o, *m3e, *m3o;
  m1e = (QSU3M *)M1->even;
  m1o = (QSU3M *)M1->odd;
  m2e = (QSU3M *)M2->even;
  m2o = (QSU3M *)M2->odd;
  m3e = (QSU3M *)M3->even;
  m3o = (QSU3M *)M3->odd;
  __assume_aligned(m1e, VECLENBYTES);
  __assume_aligned(m1o, VECLENBYTES);
  __assume_aligned(m2e, VECLENBYTES);
  __assume_aligned(m2o, VECLENBYTES);
  __assume_aligned(m3e, VECLENBYTES);
  __assume_aligned(m3o, VECLENBYTES);

#pragma omp parallel for private(i,j,v,k)
  for(s=0; s<(qphix_even_sites_on_node/VECLEN); s++)
  {
    for(i=0; i<3; i++)
    {
      for(k=0; k<3; k++)
      {
        for(j=0; j<3; j++)
        {
          for(v=0; v<VECLEN; v++)
          {
            m3e[s][i][j][0][v] -=  m2e[s][i][k][0][v] * m1e[s][k][j][0][v] - m2e[s][i][k][1][v] * m1e[s][k][j][1][v];
            m3e[s][i][j][1][v] -=  m2e[s][i][k][0][v] * m1e[s][k][j][1][v] + m2e[s][i][k][1][v] * m1e[s][k][j][0][v];
            m3o[s][i][j][0][v] -=  m2o[s][i][k][0][v] * m1o[s][k][j][0][v] - m2o[s][i][k][1][v] * m1o[s][k][j][1][v];
            m3o[s][i][j][1][v] -=  m2o[s][i][k][0][v] * m1o[s][k][j][1][v] + m2o[s][i][k][1][v] * m1o[s][k][j][0][v];

          }
        }
      }
    }
  }
}

/* operations on exposed matrices */
void QSU3_M_eq_M_times_Ma ( qphix_su3_matrix *restrict r, qphix_su3_matrix *restrict a, qphix_su3_matrix *restrict b)
{
  for(int i_c=0; i_c<3; i_c++) {
    for(int j_c=0; j_c<3; j_c++) {
      QPHIX_Complex x;
      QPHIX_c_eq_r(x,0.);
      for(int k_c=0; k_c<3; k_c++) {
        QPHIX_c_peq_c_times_ca(x, QPHIX_elem_M(*a,i_c,k_c), QPHIX_elem_M(*b,j_c,k_c));
      }
      QPHIX_c_eq_c(QPHIX_elem_M(*r,i_c,j_c),x);
    }
  }
}

void QSU3_M_eq_Ma(qphix_su3_matrix *restrict r, qphix_su3_matrix *restrict a)
{
#ifdef HAVE_XLC  
#pragma disjoint(*r,*a)
  __alignx(16,r);
  __alignx(16,a);
#endif  
  for(int i_c=0; i_c<3; i_c++) {
    for(int j_c=0; j_c<3; j_c++) {
      QPHIX_c_eq_ca(QPHIX_elem_M(*r,i_c,j_c),QPHIX_elem_M(*a,j_c,i_c));
    }
  }
}

void QSU3_M_eq_zero ( qphix_su3_matrix *restrict r)
{
#ifdef HAVE_XLC
  __alignx(16,r);
#endif
  for(int i_c=0; i_c<3; i_c++) {
    for(int j_c=0; j_c<3; j_c++) {
      QPHIX_c_eq_r(QPHIX_elem_M(*r,i_c,j_c),0.);
    }
  }
}

void QSU3_M_eq_M_times_M ( qphix_su3_matrix *restrict r, qphix_su3_matrix *restrict a, qphix_su3_matrix *restrict b)
{
#ifdef HAVE_XLC
#pragma disjoint(*r,*a,*b)
  __alignx(16,r);
  __alignx(16,a);
  __alignx(16,b);
#endif
  for(int i_c=0; i_c<3; i_c++) {
    for(int j_c=0; j_c<3; j_c++) {
      QPHIX_Complex x;
      QPHIX_c_eq_r(x,0.);
      for(int k_c=0; k_c<3; k_c++) {
        QPHIX_c_peq_c_times_c(x, QPHIX_elem_M(*a,i_c,k_c), QPHIX_elem_M(*b,k_c,j_c));
      }
      QPHIX_c_eq_c(QPHIX_elem_M(*r,i_c,j_c),x);
    }
  }
}

void QSU3_M_eq_Ma_times_M ( qphix_su3_matrix *restrict r, qphix_su3_matrix *restrict a, qphix_su3_matrix *restrict b)
{
#ifdef HAVE_XLC
#pragma disjoint(*r,*a,*b)
  __alignx(16,r);
  __alignx(16,a);
  __alignx(16,b);
#endif
  for(int i_c=0; i_c<3; i_c++) {
    for(int j_c=0; j_c<3; j_c++) {
      QPHIX_Complex x;
      QPHIX_c_eq_r(x,0.);
      for(int k_c=0; k_c<3; k_c++) {
        QPHIX_c_peq_ca_times_c(x, QPHIX_elem_M(*a,k_c,i_c), QPHIX_elem_M(*b,k_c,j_c));
      }
      QPHIX_c_eq_c(QPHIX_elem_M(*r,i_c,j_c),x);
    }
  }
}

void QSU3_M_eq_r_times_M ( qphix_su3_matrix *restrict r, QPHIX_Real *restrict a, qphix_su3_matrix *restrict b)
{
#ifdef HAVE_XLC
#pragma disjoint(*r,*a,*b)
  __alignx(16,r);
  __alignx(16,a);
  __alignx(16,b);
#endif
  for(int i_c=0; i_c<3; i_c++) {
    for(int j_c=0; j_c<3; j_c++) {
      QPHIX_c_eq_r_times_c(QPHIX_elem_M(*r,i_c,j_c), *a, QPHIX_elem_M(*b,i_c,j_c));
    }
  }
}

void QSU3_M_eq_r_times_M_plus_M ( qphix_su3_matrix *restrict r, QPHIX_Real *restrict a, qphix_su3_matrix *restrict b, qphix_su3_matrix *restrict c)
{
#ifdef HAVE_XLC
#pragma disjoint(*r,*a,*b,*c)
  __alignx(16,r);
  __alignx(16,a);
  __alignx(16,b);
  __alignx(16,c);
#endif
  for(int i_c=0; i_c<3; i_c++) {
    for(int j_c=0; j_c<3; j_c++) {
      QPHIX_c_eq_r_times_c_plus_c(QPHIX_elem_M(*r,i_c,j_c), *a, QPHIX_elem_M(*b,i_c,j_c), QPHIX_elem_M(*c,i_c,j_c));
    }
  }
}

void QSU3_R_eq_re_trace_M ( QPHIX_Real *restrict r, qphix_su3_matrix *restrict a)
{
#ifdef HAVE_XLC
#pragma disjoint(*r,*a)
  __alignx(16,r);
  __alignx(16,a);
#endif
  QPHIX_Real x;
  x = 0.;
  for(int i_c=0; i_c<3; i_c++) {
    QPHIX_r_peq_Re_c(x,QPHIX_elem_M(*a,i_c,i_c));
  }
  *r = x;
}


