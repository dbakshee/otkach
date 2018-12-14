#ifndef SUSHKOV_H
#define SUSHKOV_H

//================== BEGIN COPY FOR MINT ============================
//#define real float
#define real double
struct Params
{
    real t_ins, t_split_gate, t_w, t1_AlGaAs, t2_AlGaAs, t_GaAs, t2_GaAs;
    real xmax, ymax, ymin, xmin,hx, hy, hz, hz_gate, hz_ins, hzw, hzb;
    int  ixmin, ixmax, iymin, iymax, iz_ins1, iz_ins2, izc, iz1, iz2, izb;
    real Vtg, Vsg, eps_gaalas, eps_gaas, eps_ins;
    int NX, NY, NZ;
};

//extern struct Params params; // defined and initialized in main()

#if defined(_OPENMP)
#include <omp.h>
#define OMP(sth) __pragma(omp sth)
#else /**/
#define OMP(sth) /*none*/
#endif /**/

#if defined __INTEL_COMPILER
#define SIMD __pragma(simd)
#else /**/
#define SIMD                    /*nothing */
#endif /**/

#if !defined(M_PI)
#define M_PI 3.1415926535897932384626433832795
#endif /**/

#define ONEOVER3 (1.0/3.0)
//================== END COPY FOR MINT ============================

//================== HOW DO WE ADDRESS ARRAYS =====================
// Use column-major Fortran style
// REAL U0(0:II-1,0:JJ-1,0:KK-1) --> u0(i,j,k) = u0[i+II*j+II*JJ*k]
#if 1
//#define aimp1(i,j)   aimp1[((i)-ixmin) + ((j)-iymin)*(ixmax-ixmin+1)]
//#define aimp2(i,j)   aimp2[((i)-ixmin) + ((j)-iymin)*(ixmax-ixmin+1)]
#define Ez0(i,j)   Ez0[((i)-ixmin) + ((j)-iymin)*(ixmax-ixmin+1) ]
#define aimp1(i,j)   aimp1[((i)-ixmin) + ((j)-iymin)*(ixmax-ixmin+1) ]
#define Ez1(i,j)   Ez1[((i)-ixmin) + ((j)-iymin)*(ixmax-ixmin+1) ]
#define uold(i,j,k) uold[((i)-ixmin) + ((j)-iymin)*(ixmax-ixmin+1) + (k)*(ixmax-ixmin+1)*(iymax-iymin+1)]
#define u0(i,j,k) u0[((i)-ixmin) + ((j)-iymin)*(ixmax-ixmin+1) + (k)*(ixmax-ixmin+1)*(iymax-iymin+1)]
#define u1(i,j,k) u1[((i)-ixmin) + ((j)-iymin)*(ixmax-ixmin+1) + (k)*(ixmax-ixmin+1)*(iymax-iymin+1)]
#else
inline real& at3(int i,int j,int k,real *a,int ixmin,int ixmax,int iymin,int iymax, int izmax)
{
    int ofs = (i) + (j)*(ixmax-ixmin+1) + (k)*(ixmax-ixmin+1)*(iymax-iymin+1);
    int hi = (ixmax-ixmin+1)*(iymax-iymin+1)*(izmax);
    if (ofs < 0 || ofs >= hi
       || i<0 || i >= ixmax-ixmin+1
       || j<0 || j >= iymax-iymin+1
       || k<0 || k >= izmax )
    {
        printf("oops in at3: ofs=%i but hi=%i, ijk=%i,%i,%i\n",ofs,hi,i,j,k);
        exit(1);
    }
    return a[ofs];
}
inline real& at2(int i,int j,real *a,int ixmin,int ixmax,int iymin,int iymax)
{
    int ofs = (i) + (j)*(ixmax-ixmin+1);
    int hi = (ixmax-ixmin+1)*(iymax-iymin+1);
    if (ofs < 0 || ofs >= hi
       || i<0 || i >= ixmax-ixmin+1
       || j<0 || j >= iymax-iymin+1)
    {
        printf("oops in at2: ofs=%i but hi=%i, ij=%i,%i\n",ofs,hi,i,j);
        exit(1);
    }
    return a[ofs];
}
inline int& at2i(int i,int j,int *a,int ixmin,int ixmax,int iymin,int iymax)
{
    int ofs = (i) + (j)*(ixmax-ixmin+1);
    int hi = (ixmax-ixmin+1)*(iymax-iymin+1);
    if (ofs < 0 || ofs >= hi
       || i<0 || i >= ixmax-ixmin+1
       || j<0 || j >= iymax-iymin+1)
    {
        printf("oops in at2i: ofs=%i but hi=%i, ij=%i,%i\n",ofs,hi,i,j);
        exit(1);
    }
    return a[ofs];
}
//#define aimp1(i,j)    at2(i,j,aimp1,ixmin,ixmax,iymin,iymax)
//#define aimp2(i,j)    at2(i,j,aimp2,ixmin,ixmax,iymin,iymax)
#define Ez0(i,j)      at2(i,j,Ez0,ixmin,ixmax,iymin,iymax)
#define aimp1(i,j)   at2(i,j,aimp1,ixmin,ixmax,iymin,iymax)
#define Ez1(i,j)      at2(i,j,Ez0,ixmin,ixmax,iymin,iymax)
#define uold(i,j,k)     at3(i,j,k,u0,ixmin,ixmax,iymin,iymax,izb+1)
#define u0(i,j,k)     at3(i,j,k,u0,ixmin,ixmax,iymin,iymax,izb+1)
#define u1(i,j,k)     at3(i,j,k,u1,ixmin,ixmax,iymin,iymax,izb+1)
#endif


//================== FUNCTION DECLARATIONS ========================

#if defined(__cplusplus)
extern "C"
#endif
int compute_cuda(real ***u0,    /* U[n3][n2][n1] */
                 real ***u1,
                 int nrounds,   /* one round is U->Unew, and swap */
                 double *nflops, /* number of floating point operations */
                 struct Params *params);

int compute_host(real *u0,  /* U(n1,n2,n3) */
                 real *u1,
                 real *Ez0, 
                 real *aimp1, 
                 int nrounds,   /* one round is U->Unew, and swap */
                 double *nflops, /* number of floating point operations */
                 struct Params *params);

#endif
