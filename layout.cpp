#include <omp.h>
#include <stdio.h>
#include <sys/time.h>
#include "qphix_internal.h"
#include "ks_globals.h"
#include "ff_boundary.h"
//#define COMM
//#define COMM2
//#define COMM3
#define COMM4

extern "C"
{
#include "layout.h"

extern site *lattice;
extern int qphix_sites_on_node;
extern int qphix_even_sites_on_node;
extern int qphix_fused_sites_on_node;

static const int nGX = (VECLEN < 2  ? 1 : 2);
static const int nGY = (VECLEN < 4  ? 1 : 2);
static const int nGZ = (VECLEN < 8  ? 1 : 2);
static const int nGT = (VECLEN < 16 ? 1 : 2);

#ifdef FF_DEBUG
void QPHIX_checksum_M(QPHIX_ColorMatrix *m, int verbosity)
{
  int i;
  double *d;
  SU3_Matrix *data;
  double sum[15] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  data =  (SU3_Matrix *) _mm_malloc(qphix_sites_on_node * sizeof(SU3_Matrix), VECLENBYTES);
  
  QPHIX_layout_to_su3m(data, m);
  d = (double *) (data);
  
  for(i=0; i<18; i++)
  {
    sum[0]  += d[i];                  //sum of all real & imag
    if(i%2 == 0)
    {
      sum[1]  += d[i];                //sum of all real
    if(i/6 == 0)     sum[3]  += d[i]; //sum of real row1
    if(i/6 == 1)     sum[4]  += d[i]; //sum of real row2
    if(i/6 == 2)     sum[5]  += d[i]; //sum of real row3
    if(i%6 == 0)     sum[9]  += d[i]; //sum of real col1
    if((i-2)%6 == 0) sum[10] += d[i]; //sum of real col2
    if((i-4)%6 == 0) sum[11] += d[i]; //sum of real col3
    }
    else
    {
      sum[2]  += d[i];                  //sum of all imag
      if(i/6 == 0)     sum[6]  += d[i]; //sum of imag row1
      if(i/6 == 1)     sum[7]  += d[i]; //sum of imag row2
      if(i/6 == 2)     sum[8]  += d[i]; //sum of imag row3
      if((i-1)%6 == 0) sum[12] += d[i]; //sum of imag col1
      if((i-3)%6 == 0) sum[13] += d[i]; //sum of imag col2
      if((i-5)%6 == 0) sum[14] += d[i]; //sum of imag col3
    }
  }
  
  if (verbosity == 0)
  {
    printf("Real matrix:\n");
    for(i=0; i<18; i=i+2)
    {
      if(i==6 || i==12) printf("\n");
      printf("%e  ", *(d+i));
    }
    printf("\n");
    
    printf("Imag matrix:\n");
    for(i=1; i<18; i=i+2)
    {
      if(i==7 || i==13) printf("\n");
      printf("%e  ", *(d+i));
    }
    printf("\n");
  }
  printf("Checksum:\n");
  if(verbosity == 1)
  {
    printf("{%e}\n", sum[0]);
  }
  if(verbosity == 2)
  {
    printf("{%e, %e, %e}\n", sum[0], sum[1], sum[2]);
  }
  if(verbosity == 3)
  {
    printf("{{%e, %e, %e},\n", sum[0], sum[1], sum[2]);
    printf(" {%e, %e, %e},\n", sum[3], sum[4], sum[5]);
    printf(" {%e, %e, %e},\n", sum[6], sum[7], sum[8]);
    printf(" {%e, %e, %e},\n", sum[9], sum[10], sum[11]);
    printf(" {%e, %e, %e}}\n", sum[12], sum[13], sum[14]);
  }
  
  _mm_free(data);
}
#endif


void QPHIX_get_milc_nbr_crd(int x, int y, int z, int t, int dir, int *xp,
                            int *yp, int *zp, int *tp)
{
  /*
  x, y, z, t         => coordinates of site
  dir                => direction (e.g. XUP)
  *xp, *yp, *zp, *tp => pointers to coordinates of neighbor
  */
  
#ifdef DEBUG
  assert((dir >= XUP) && (dir <= X3DOWN));
  assert((x >= 0) && (x < Nx));
  assert((y >= 0) && (y < Ny));
  assert((z >= 0) && (z < Nz));
  assert((t >= 0) && (t < Nt));
#endif
  
  *xp = x; *yp = y; *zp = z; *tp = t;
  
  switch(dir)
  {
    case XUP   : *xp = (x+1)%Nx; break;
    case XDOWN : *xp = (x+Nx-1)%Nx; break;
    case YUP   : *yp = (y+1)%Ny; break;
    case YDOWN : *yp = (y+Ny-1)%Ny; break;
    case ZUP   : *zp = (z+1)%Nz; break;
    case ZDOWN : *zp = (z+Nz-1)%Nz; break;
    case TUP   : *tp = (t+1)%Nt; break;
    case TDOWN : *tp = (t+Nt-1)%Nt; break;
    case X3UP  : *xp = (x+3)%Nx; break;
    case X3DOWN: *xp = (x+4*Nx-3)%Nx; break;
    case Y3UP  : *yp = (y+3)%Ny; break;
    case Y3DOWN: *yp = (y+4*Ny-3)%Ny; break;
    case Z3UP  : *zp = (z+3)%Nz; break;
    case Z3DOWN: *zp = (z+4*Nz-3)%Nz; break;
    case T3UP  : *tp = (t+3)%Nt; break;
    case T3DOWN: *tp = (t+4*Nt-3)%Nt; break;
    default    : printf("BOTCH: bad direction\n"); exit(1);
  }
}

void QPHIX_get_nbr_index(int ind, int v, int p, int dir, int *nind, int *nv,
                         int *np)
{
  /*
  Given a qphix array index and vector index, find the index and vector
  offset of my neighbor in the specified direction.
  Args:
  ind    - index into qphix fused array
  v      - vector offset into qphix fused array
  p      - parity, even or odd site
  dim    - dimension in which you want to find neighbor
           will be in the range [XUP,TUP] or [X3UP,T3UP]
  fb     - forwards or backwards in that dimension
  nind   - Neighbor index into qphix fused array
  nv     - Neighbor's vector offset
  np     - Neighbor's parity
  */
  
#ifdef DEBUG
  assert(ind <= qphix_fused_sites_on_node);
  assert((v >=0) && (v <= VECLEN));
  assert((dir >= XUP) && (dir <= X3DOWN));
#endif
  
  int tind = ind, tv = v, tnx, tdir;
  int x, y, z, t;
  int x1, y1, z1, t1;
  int x2, y2, z2, t2;
  int nx, ny, nz, nt;
  
    /* ind = x2 + Vxh*y2 + Vxh*Vy*z2 + Vxh*Vy*Vz*t2 */
    /* v = x1 + nGX*y1 + nGX*nGY*z1 + nGX*nGY*nGZ*t1 */

//    printf("Pxy %d Pxyz %d qphix_fused_sites_on_node %d\n",Pxy,Pxyz,qphix_fused_sites_on_node);
    x2 = ind % Vxh;
    x1 = v % nGX;
    x  = x1*Vxh + x2;

    //t1 = ((ind - (ind % Vxh)) / Vxh) % Vy
    //t2 = ((v   - (v   % nGx)) / nGx) % nGy
    //y  = t1 + (t2 * Vy)
    ind -= x2;    /* ind = Vxh*y2 + Vxh*Vy*z2 + Vxh*Vy*Vz*t2 */
    ind /= Vxh;   /* ind = y2 + Vy*z2 + Vy*Vz*z2 */
    y2 = ind % Vy;
    v -= x1;      /* v = nGX*y1 + nGX*nGY*z1 + nGX*nGY*nGZ*t1 */
    v /= nGX;     /* v = y1 + nGY*z1 + nGY*nGZ*t1 */
    y1 = v % nGY;
    y  = y1*Vy + y2;

    //t1 = ((ind - (ind % Vy )) / Vy) % Vz
    //t2 = ((v   - (v   % nGy)) / nGy) % nGz
    //z  = t1 + (t2 * Vz)
    ind -= y2;    /* ind = Vy*z2 + Vy*Vz*t2 */
    ind /= Vy;    /* ind = z2 + Vz*t2 */
    z2 = ind % Vz;
    v -= y1;      /* v = nGY*z1 + nGY*nGZ*t1 */
    v /= nGY;     /* v = z1 + nGZ*t1 */
    z1 = v % nGZ;
    z  = z1*Vz + z2;

    ind -= z2;    /* ind = Vz*t2 */
    t2 = ind / Vz;
    v -= z1;      /* v = nGZ*t1 */
    t1 = v / nGZ;
    t  = t1*Vy + t2;

    /*
     * (x,y,z,t) is the sub-lattice coords.
     * Convert it to global lattice coords now.
     * p determines whether we are looking at the even or the odd array
     *
     * Even lattice
     * E o E o E o E o
     * o E o E o E o E
     * E o E o E o E o
     * o E o E o E o E <- x must be 2*x+1 when (y+z+t) is odd
     * E o E o E o E o <- x must be 2*x   when (y+z+t) is even
     *
     * Odd lattice
     * e O e O e O e O
     * O e O e O e O e
     * e O e O e O e O
     * O e O e O e O e <- x must be 2*x   when (y+z+t) is odd
     * e O e O e O e O <- x must be 2*x+1 when (y+z+t) is even
     */
    if (p == QPHIX_EVEN) {
      x = 2*x + ((y+z+t) % 2);
    } else if (p == QPHIX_ODD) {
      x = 2*x + ((y+z+t+1) % 2);
    }

    QPHIX_get_milc_nbr_crd(x,y,z,t,dir,&nx,&ny,&nz,&nt);
    tnx = nx;

    if(((nx+ny+nz+nt) % 2) == 0) { /* neighbor is even */
      *np = QPHIX_EVEN;
      nx = ((ny+nz+nt) % 2) ? (nx-1)/2 : (nx/2);
    } else { /* neighbor is odd */
      *np = QPHIX_ODD;
      nx = ((ny+nz+nt) % 2) ? (nx/2) : (nx-1)/2;
    }

    x1 = nx / Vxh;
    x2 = nx % Vxh;
    y1 = ny / Vy;
    y2 = ny % Vy;
    z1 = nz / Vz;
    z2 = nz % Vz;
    t1 = nt / Vt;
    t2 = nt % Vt;

    *nind = t2 * Pxyz + z2 * Pxy + y2 * Vxh + x2;
    *nv = ((t1 * nGZ + z1) * nGY + y1) * nGX + x1;
}
  
void QPHIX_V_eq_sV(QPHIX_ColorVector *v2, QPHIX_ColorVector *v1, int dim,
                   int fb)
{
#ifdef DEBUG
  assert((v1 != NULL) && (v1->even != NULL) && (v1->odd != NULL));
  assert((v2 != NULL) && (v2->even != NULL) && (v2->odd != NULL));
  assert(((dim >= XUP)  && (dim <= TUP))  ||
        ((dim >= X3UP) && (dim <= T3UP)));
  assert((fb == FORWARDS) || (fb == BACKWARDS));
#endif
  QSU3V *v1e, *v1o, *v2e, *v2o, *vs;
  int i, v, p, el;
  int ni, nv, np;
  int dir;
// fbv is used to indicate forward or backward to QPhiX function calls. In MILC this is +1 or -1. QPhiX expects 0 or 1.
  int fbv = (fb == 1)? 0:1;
  
  v1e = (QSU3V *)v1->even;
  v1o = (QSU3V *)v1->odd;
  v2e = (QSU3V *)v2->even;
  v2o = (QSU3V *)v2->odd;
  
  dir = dim;
  if (fb == BACKWARDS)
  {
    if (dim <= TUP) dir = OPP_DIR(dim);
    else dir = OPP_3_DIR(dim);
  }
  
  p = QPHIX_EVEN;
// the ifdef COMM, 2, 3 and 4 are defined for debug convenience. these correspond to the progressive communication in 4 directions.
#ifdef COMM
  if (dir == X3UP){
    ff_pack_and_send_boundaries_v3(0, v1o, 1, fbv, (dim-8));
  }
  else if((dir == XUP)  || (dir == XDOWN))
    ff_pack_and_send_boundaries_v(0, v1o, 1, fbv, dim);
#endif
#ifdef COMM2
  if (dir == X3UP || dir == Y3UP)
    ff_pack_and_send_boundaries_v3(0, v1o, 1, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN))
    ff_pack_and_send_boundaries_v(0, v1o, 1, fbv, dim);
#endif
#ifdef COMM3
  if (dir == X3UP || dir == Y3UP || dir == Z3UP)
    ff_pack_and_send_boundaries_v3(0, v1o, 1, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN))
    ff_pack_and_send_boundaries_v(0, v1o, 1, fbv, dim);
#endif
#ifdef COMM4
  if (dir == X3UP || dir == Y3UP || dir == Z3UP || dir == T3UP)
    ff_pack_and_send_boundaries_v3(0, v1o, 1, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN) || (dir == TUP) || (dir == TDOWN))
    ff_pack_and_send_boundaries_v(0, v1o, 1, fbv, dim);
#endif

#pragma omp parallel for private(v, ni, nv, np, vs, el)
  for(i=0; i<qphix_fused_sites_on_node; i++)
  {
    for(v=0; v<VECLEN; v++)
    {
      QPHIX_get_nbr_index(i, v, p, dir, &ni, &nv, &np);
      vs = (np == QPHIX_EVEN) ? v1e : v1o;
      for(el=0;el<3;++el)
      {
        v2e[i][el][0][v] = vs[ni][el][0][nv];
        v2e[i][el][1][v] = vs[ni][el][1][nv];
      }
    }
  }

#ifdef COMM
  if (dir == X3UP)
    ff_recv_and_unpack_boundaries_v3(0, v2e, 1, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN))
    ff_recv_and_unpack_boundaries_v(0, v2e, 1, fbv, dim);
#endif
#ifdef COMM2
  if (dir == X3UP || dir == Y3UP)
    ff_recv_and_unpack_boundaries_v3(0, v2e, 1, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN))
    ff_recv_and_unpack_boundaries_v(0, v2e, 1, fbv, dim);
#endif
#ifdef COMM3
  if (dir == X3UP || dir == Y3UP || dir == Z3UP)
    ff_recv_and_unpack_boundaries_v3(0, v2e, 1, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN))
    ff_recv_and_unpack_boundaries_v(0, v2e, 1, fbv, dim);
#endif
#ifdef COMM4
  if (dir == X3UP || dir == Y3UP || dir == Z3UP || dir == T3UP)
    ff_recv_and_unpack_boundaries_v3(0, v2e, 1, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN) || (dir == TUP) || (dir == TDOWN)){
    ff_recv_and_unpack_boundaries_v(0, v2e, 1, fbv, dim);
    }
#endif
  
  p = QPHIX_ODD;
#ifdef COMM
  if (dir == X3UP)
    ff_pack_and_send_boundaries_v3(0, v1e, 0, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN))
    ff_pack_and_send_boundaries_v(0, v1e, 0, fbv, dim);
#endif
#ifdef COMM2
  if (dir == X3UP || dir == Y3UP)
    ff_pack_and_send_boundaries_v3(0, v1e, 0, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN))
    ff_pack_and_send_boundaries_v(0, v1e, 0, fbv, dim);
#endif
#ifdef COMM3
  if (dir == X3UP || dir == Y3UP || dir == Z3UP)
    ff_pack_and_send_boundaries_v3(0, v1e, 0, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN))
    ff_pack_and_send_boundaries_v(0, v1e, 0, fbv, dim);
#endif
#ifdef COMM4
  if (dir == X3UP || dir == Y3UP || dir == Z3UP || dir == T3UP)
    ff_pack_and_send_boundaries_v3(0, v1e, 0, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN) || (dir == TUP) || (dir == TDOWN))
    ff_pack_and_send_boundaries_v(0, v1e, 0, fbv, dim);
#endif

#pragma omp parallel for private(v, ni, nv, np, vs, el)
  for(i=0; i<qphix_fused_sites_on_node; i++)
  {
    for(v=0; v<VECLEN; v++)
    {
      QPHIX_get_nbr_index(i, v, p, dir, &ni, &nv, &np);
      vs = (np == QPHIX_EVEN) ? v1e : v1o;
      for(el=0;el<3;++el)
      {
        v2o[i][el][0][v] = vs[ni][el][0][nv];
        v2o[i][el][1][v] = vs[ni][el][1][nv];
      }
    }
  }

#ifdef COMM
  if (dir == X3UP)
    ff_recv_and_unpack_boundaries_v3(0, v2o, 0, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN))
    ff_recv_and_unpack_boundaries_v(0, v2o, 0, fbv, dim);
#endif
#ifdef COMM2
  if (dir == X3UP || dir == Y3UP)
    ff_recv_and_unpack_boundaries_v3(0, v2o, 0, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN))
    ff_recv_and_unpack_boundaries_v(0, v2o, 0, fbv, dim);
#endif
#ifdef COMM3
  if (dir == X3UP || dir == Y3UP || dir == Z3UP)
    ff_recv_and_unpack_boundaries_v3(0, v2o, 0, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN))
    ff_recv_and_unpack_boundaries_v(0, v2o, 0, fbv, dim);
#endif
#ifdef COMM4
  if (dir == X3UP || dir == Y3UP || dir == Z3UP || dir == T3UP)
    ff_recv_and_unpack_boundaries_v3(0, v2o, 0, fbv, (dim-8));
  else if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN) || (dir == TUP) || (dir == TDOWN))
    ff_recv_and_unpack_boundaries_v(0, v2o, 0, fbv, dim);
#endif

}

void QPHIX_M_eq_sM(QPHIX_ColorMatrix *m2, QPHIX_ColorMatrix *m1, int dim, int fb)
{
#ifdef DEBUG
  assert((m1 != NULL) && (m1->even != NULL) && (m1->odd != NULL));
  assert((m2 != NULL) && (m2->even != NULL) && (m2->odd != NULL));
  assert(((dim >= XUP)  && (dim <= TUP))  ||
      ((dim >= X3UP) && (dim <= T3UP)));
  assert((fb == FORWARDS) || (fb == BACKWARDS));
#endif

  QSU3M *m1e, *m1o, *m2e, *m2o, *ms;
  int i, v, p, row, col;
  int ni, nv, np;
  int dir;
  int fbv = (fb == 1)? 0:1;
#ifdef FF_SHIFT_PROFILE  
  static double shift_total=0.0;
  struct timeval shift_start, shift_end;

  gettimeofday(&shift_start, NULL);
#endif  

  m1e = (QSU3M *)m1->even;
  m1o = (QSU3M *)m1->odd;
  m2e = (QSU3M *)m2->even;
  m2o = (QSU3M *)m2->odd;

  dir = dim;
  if (fb == BACKWARDS)
  {
    if (dim <= TUP) dir = OPP_DIR((dim));
    else dir = OPP_3_DIR(dim);
  }

  p = QPHIX_EVEN;
#ifdef COMM4
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN) || (dir == TUP) || (dir == TDOWN))
  ff_pack_and_send_boundaries(0, m1o, 1, fbv, dim);
#endif
#ifdef COMM3
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP)  || (dir == YDOWN) || (dir == ZUP)  || (dir == ZDOWN))
  ff_pack_and_send_boundaries(0, m1o, 1, fbv, dim);
#endif
#ifdef COMM2
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP)  || (dir == YDOWN))
  ff_pack_and_send_boundaries(0, m1o, 1, fbv, dim);
#endif
#ifdef COMM
  if((dir == XUP)  || (dir == XDOWN))
  ff_pack_and_send_boundaries(0, m1o, 1, fbv, dim);
#endif

#pragma omp parallel for private(v,ni,nv,np,ms,row,col)
  for(i=0;i<qphix_fused_sites_on_node;++i)
  {
    for(v=0;v<VECLEN;++v)
    {
      QPHIX_get_nbr_index(i,v,p,dir,&ni,&nv,&np);
      ms = (np == QPHIX_EVEN) ? m1e : m1o;
      for(row=0;row<3;++row)
      {
        for(col=0;col<3;++col)
        {
          m2e[i][row][col][0][v] = ms[ni][row][col][0][nv];
          m2e[i][row][col][1][v] = ms[ni][row][col][1][nv];
        }
      }
    }
  }

#ifdef COMM4
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN) || (dir == TUP) || (dir == TDOWN))
  ff_recv_and_unpack_boundaries(0, m2e, 1, fbv, dim);
#endif
#ifdef COMM3
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP)  || (dir == YDOWN) || (dir == ZUP)  || (dir == ZDOWN))
  ff_recv_and_unpack_boundaries(0, m2e, 1, fbv, dim);
#endif
#ifdef COMM2
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP)  || (dir == YDOWN))
  ff_recv_and_unpack_boundaries(0, m2e, 1, fbv, dim);
#endif
#ifdef COMM
  if((dir == XUP)  || (dir == XDOWN))
  ff_recv_and_unpack_boundaries(0, m2e, 1, fbv, dim);
#endif
  
  p = QPHIX_ODD;
#ifdef COMM4
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN) || (dir == TUP) || (dir == TDOWN))
  ff_pack_and_send_boundaries(0, m1e, 0, fbv, dim);
#endif
#ifdef COMM3
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP)  || (dir == YDOWN) || (dir == ZUP)  || (dir == ZDOWN))
  ff_pack_and_send_boundaries(0, m1e, 0, fbv, dim);
#endif
#ifdef COMM2
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP)  || (dir == YDOWN))
  ff_pack_and_send_boundaries(0, m1e, 0, fbv, dim);
#endif
#ifdef COMM
  if((dir == XUP)  || (dir == XDOWN))
  ff_pack_and_send_boundaries(0, m1e, 0, fbv, dim);
#endif

#pragma omp parallel for private(v,ni,nv,np,ms,row,col)
  for(i=0;i<qphix_fused_sites_on_node;++i)
  {
    for(v=0;v<VECLEN;++v)
    {
      QPHIX_get_nbr_index(i,v,p,dir,&ni,&nv,&np);
      ms = (np == QPHIX_EVEN) ? m1e : m1o;
      for(row=0;row<3;++row)
      {
        for(col=0;col<3;++col)
        {
          m2o[i][row][col][0][v] = ms[ni][row][col][0][nv];
          m2o[i][row][col][1][v] = ms[ni][row][col][1][nv];
        }
      }
    }
  }

#ifdef COMM4
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP) || (dir == YDOWN) || (dir == ZUP) || (dir == ZDOWN) || (dir == TUP) || (dir == TDOWN))
  ff_recv_and_unpack_boundaries(0, m2o, 0, fbv, dim);
#endif
#ifdef COMM3
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP)  || (dir == YDOWN) || (dir == ZUP)  || (dir == ZDOWN))
  ff_recv_and_unpack_boundaries(0, m2o, 0, fbv, dim);
#endif
#ifdef COMM2
  if((dir == XUP)  || (dir == XDOWN) || (dir == YUP)  || (dir == YDOWN))
  ff_recv_and_unpack_boundaries(0, m2o, 0, fbv, dim);
#endif
#ifdef COMM
  if((dir == XUP)  || (dir == XDOWN))
  ff_recv_and_unpack_boundaries(0, m2o, 0, fbv, dim);
#endif

#ifdef FF_SHIFT_PROFILE  
  gettimeofday(&shift_end, NULL);
  shift_total += (shift_end.tv_sec + (double) shift_end.tv_usec/1000000) -
                 (shift_start.tv_sec + (double) shift_start.tv_usec/1000000);
  printf("cumulative shift time = %f\n", shift_total);
#endif  
}

}
