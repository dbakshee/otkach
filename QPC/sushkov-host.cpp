// -*- c -*-
// Trying 3d
#include <stdio.h>
//#include <malloc.h>
#include <math.h>
#include <xmmintrin.h>
#define _rdtsc __rdtsc
#include "ringnew2.h"
#include <vector>
#include <cstdlib>

#if defined(_WIN32)
#define ZI "%Ii"
#else
#define ZI "%zi"
#endif
/*
#include <stdio.h>
//#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <xmmintrin.h>
#include "ringnew2.h"

#if defined(_WIN32)
#define ZI "%Ii"
#else
#define ZI "%zi"
#endif
*/
inline int splitGate(real xx, real yy)
{
  real Lsg=600;
  real xg=0.5*Lsg+50;
  real yg=200;//150;
  real x0=0.5*Lsg;
  real y0=yg+50;
  real x=abs(xx);
  real y=abs(yy);
  real r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
  int res=0;
  if(x<=xg&&y>=y0) {res=1;}
  else
    {
    if(y<y0&&y>=yg&&x<=x0) res=1;
    else if(y<y0&&y>=yg&&x>x0&&x<xg&&r<=50) res=1;
    }
  return res;
}
inline int topGate(real xx, real yy, real t_ins)
{
  real Lsg=600;
  real xg=0.5*Lsg+50+t_ins-5;
  real yg=200-t_ins+5;
  real x0=0.5*Lsg;
  real y0=yg+50;
  real x=abs(xx);
  real y=abs(yy);
  real r = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
  int res=0;
  if(x>xg) {res=1;}
  if(y<yg&&x<=xg&&r>70) res=1;
  res=1;//!!!!
  return res;
}


int compute_host(real *u0,      /* U(n1,n2,n3) */
                 real *u1,
                 real *Ez0,
                 real *aimp1,
                 int nrounds,   /* one round is U->Unew and swap */
                 double *nflops, /* number of floating point operations */
                 struct Params *params)
{
    real hx = params->hx;
    real hy = params->hy;
    real hz = params->hz;
    real hzw = params->hzw;
    real hzb = params->hzb;
    real hz_ins = params->hz_ins;
    real hz_gate = params->hz_gate;

    real eps_gaalas = params->eps_gaalas;
    real eps_gaas = params->eps_gaas;
    real eps_ins = params->eps_ins;

    real t_ins = params->t_ins;
    int izb = params->izb;
    int iz1 = params->iz1;
    int iz2 = params->iz2;
    int izc = params->izc;
    int iz_ins1 = params->iz_ins1;
    int iz_ins2 = params->iz_ins2;

    real ymin = params->ymin;
    real xmin = params->xmin;
    int iymin = params->iymin;
    int iymax = params->iymax;
    int ixmin = params->ixmin;
    int ixmax = params->ixmax;
    real Vtg = params->Vtg;
    real Vsg = params->Vsg;
    real dNs = 0;//-1.875;// 10^11 cm^-2
    real AA1 = 4.*M_PI*4.803*4.803*1e-4/1.6;
    real efermi = 0;
    real ebandoffset = 0.152;//0.7*0.46;//0.275*2./3.;//0.24;
    real me = 0.067;
//    real mh = 0.3;
    real dw = params->t_w;//16.;//20;//5;
    real E_1=(M_PI/dw)*(M_PI/dw)*0.0381/me;//0.02;//0.01;
    real delta=0.26;
    real Eg=1.5;
    real fn2D = 10.*AA1/eps_gaas;
    real A4 = 324.*me/0.07766*AA1/eps_gaas;
    real LL = 1000;
    int NX = params->NX;
    int NY = params->NY;
    int NZ = params->NZ;
    real hxhx = 1. / (hx * hx);
    real hyhy = 1. / (hy * hy);
    real hzhz = 1. / (hz * hz);
////////////////////////////////////////
// Для рандомизированных примесей
//       Число примесей $N_S$ на площади $S$~nm$^2$ при концентрации
//       $n_s\cdot10^{11}$~см$^{-2}$ есть $$N_S=0.001 n_s S$$
//        parameter (S_spacer=hy*(iymaxspacer-iyminspacer)*(ixmaxspacer-ixminspacer)*hx)
//        parameter (iS_spacer=(iymaxspacer-iyminspacer+1)*(ixmaxspacer-ixminspacer+1))!???
//        parameter (nimp1=0.001*dNdSi1surf*S_spacer) !total number of single charge impurities in layer
//!       Перенормированная плотность заряда в примесях
//        parameter (ASi1=nimp*AA1/eps_gaas)
//        parameter (ASi1=-dNdSi1*4*pi*(4.803)**2*1e-4/1.6/eps_gaalas)
//        parameter (Asi1imp=ASi1*iS_spacer/nimp1)
//////////////////////////////////////
    real nimp=-7.5;//3;//1;//2;//3;//3;//10^11cm^{-2}=10^{-3}nm^{-2}
    real ns=-7.5;//0;//0;//3.43-nimp;//1.77;//0.65;//73;//525;
    real Aimp = 0;//1000*AA1/eps_gaas/hz/hx/hy;// вместо Asi1
    real s=1./(hx*hy);
    real Utg=-Vtg+delta;//Eg/2;//delta+0.05;//+0.01;//Eg/2;
    real Usg=-Vsg+delta;//Eg/2;//delta;//Eg/2;
    real c_sh_U=0;//6;//5.84853;// шунтирующая емкость поверхностных состояний по отношению к U_new()-Delta
    real hxhyhz = 1. / (2. * hzhz + 2. * hxhx + 2. * hyhy);
    real hzhzw = 1. / (hzw*hzw);
    real hxhyhzw = 1. / (2.*hzhzw+2.*hxhx+2.*hyhy);
    real f3 = A4*hxhyhzw*M_PI*0.5/dw;
    long long dt[15];
    for (int n = 0; n < nrounds; ++n)
    {
        //boundary conditions: vertical bounds
        dt[0] = _rdtsc();
        OMP(parallel for)
        for (int j = iymin; j <= iymax; ++j)
        {
#pragma ivdep
#pragma simd
             for (int i = ixmin; i <= ixmax; ++i)
             {
              Ez0(i,j) = u0(i,j,int(iz1+0.5*dw/hzw))+E_1; //  !!!!!!!!!!!!!!!!
             }
        }

// 0000
        OMP(parallel for)
        for (int j = iymin; j <= iymax; ++j)
        {
#pragma ivdep
#pragma simd
            for (int i = ixmin; i <= ixmax; ++i)
            {
                u0(i, j, 0) = Utg;
                u0(i,j,izb)=u0(i, j,izb-1);//-hzb*delta/LL;//(0.5*Eg)/LL;  ?????????????????
            }
        }
        dt[0] = _rdtsc()-dt[0];
        dt[1] = _rdtsc();
        int kk;
        hzhz = 1. / (hz_gate * hz_gate);
        hxhyhz = 1. / (2. * hzhz + 2. * hxhx + 2. * hyhy);

//1111
        for (int k = 1; k < iz_ins2; ++k)
        {
                u0(ixmin, iymin, k) = ((u0(ixmin, iymin, k - 1) + u0(ixmin, iymin, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymin, k) + u0(ixmin + 1, iymin, k)) * hxhx
                                       +      (u0(ixmin, iymax - 1, k) + u0(ixmin, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmin, iymax, k) = ((u0(ixmin, iymax, k - 1) + u0(ixmin, iymax, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymax, k) + u0(ixmin + 1, iymax, k)) * hxhx
                                       +      (u0(ixmin, iymax - 1, k) + u0(ixmin, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmax, iymax, k) = ((u0(ixmax, iymax, k - 1) + u0(ixmax, iymax, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymax, k) + u0(ixmin + 1, iymax, k)) * hxhx
                                       +      (u0(ixmax, iymax - 1, k) + u0(ixmax, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmax, iymin, k) = ((u0(ixmax, iymin, k - 1) + u0(ixmax, iymin, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymin, k) + u0(ixmin + 1, iymin, k)) * hxhx
                                       +      (u0(ixmax, iymax - 1, k) + u0(ixmax, iymin + 1, k)) * hyhy) * hxhyhz;
                u1(ixmax, iymin, k) = u0(ixmax, iymin, k);
                u1(ixmax, iymax, k) = u0(ixmax, iymax, k);
                u1(ixmin, iymin, k) = u0(ixmin, iymin, k);
                u1(ixmin, iymax, k) = u0(ixmin, iymax, k);
        }

                       kk=iz_ins2;
                       {
                       int i=ixmin;
                       int j=iymin;
                       real x=hx*i+xmin;
                       real y=hy*j+ymin;
                       real dU1=-1000*s*0.25*(aimp1(ixmin,iymin)+aimp1(ixmin,iymax)+aimp1(ixmax,iymax)+aimp1(ixmax,iymin))*(0.1*1.602/8.85);
                      if(!splitGate(x,y))
                        {
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hz + u0(i, j, kk - 1)*eps_ins/hz_gate)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hz + u0(i, j, kk - 2)*eps_ins/hz_gate)-2.*dU1)/
                                               (3.*(eps_gaas/hz+eps_ins/hz_gate));

                        }
                       else u0(i, j, kk) = Usg;
                       u1(ixmin, j, kk) = u0(i, j, kk);

                       u0(ixmax, j, kk) = u0(i, j, kk);
                       u1(ixmax, j, kk) = u1(i, j, kk);

                       u0(ixmax, iymax, kk) = u0(i, j, kk);
                       u1(ixmax, iymax, kk) = u1(i, j, kk);

                       u0(ixmin, iymax, kk) = u0(i, j, kk);
                       u1(ixmin, iymax, kk) = u1(i, j, kk);
                       }
    hzhz = 1. / (hz * hz);
    hxhyhz = 1. / (2. * hzhz + 2. * hxhx + 2. * hyhy);

        for (int k = iz_ins2+1; k < iz1; ++k)
        {
                u0(ixmin, iymin, k) = ((u0(ixmin, iymin, k - 1) + u0(ixmin, iymin, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymin, k) + u0(ixmin + 1, iymin, k)) * hxhx
                                       +      (u0(ixmin, iymax - 1, k) + u0(ixmin, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmin, iymax, k) = ((u0(ixmin, iymax, k - 1) + u0(ixmin, iymax, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymax, k) + u0(ixmin + 1, iymax, k)) * hxhx
                                       +      (u0(ixmin, iymax - 1, k) + u0(ixmin, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmax, iymax, k) = ((u0(ixmax, iymax, k - 1) + u0(ixmax, iymax, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymax, k) + u0(ixmin + 1, iymax, k)) * hxhx
                                       +      (u0(ixmax, iymax - 1, k) + u0(ixmax, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmax, iymin, k) = ((u0(ixmax, iymin, k - 1) + u0(ixmax, iymin, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymin, k) + u0(ixmin + 1, iymin, k)) * hxhx
                                       +      (u0(ixmax, iymax - 1, k) + u0(ixmax, iymin + 1, k)) * hyhy) * hxhyhz;
                u1(ixmax, iymin, k) = u0(ixmax, iymin, k);
                u1(ixmax, iymax, k) = u0(ixmax, iymax, k);
                u1(ixmin, iymin, k) = u0(ixmin, iymin, k);
                u1(ixmin, iymax, k) = u0(ixmin, iymax, k);
        }

//сшивка k=iz1
        kk=iz1;
        int j=iymax;
        int i = ixmin;
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzw + u0(i, j, kk - 1)*eps_gaalas/hz)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzw + u0(i, j, kk - 2)*eps_gaalas/hz))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hz));
        i = ixmax;
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzw + u0(i, j, kk - 1)*eps_gaalas/hz)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzw + u0(i, j, kk - 2)*eps_gaalas/hz))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hz));
        j=iymin;
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzw + u0(i, j, kk - 1)*eps_gaalas/hz)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzw + u0(i, j, kk - 2)*eps_gaalas/hz))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hz));
        i = ixmin;
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzw + u0(i, j, kk - 1)*eps_gaalas/hz)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzw + u0(i, j, kk - 2)*eps_gaalas/hz))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hz));
        u1(ixmax, iymin, kk) = u0(ixmax, iymin, kk);
        u1(ixmin, iymin, kk) = u0(ixmin, iymin, kk);
        u1(ixmax, iymax, kk) = u0(ixmax, iymax, kk);
        u1(ixmin, iymax, kk) = u0(ixmin, iymax, kk);
//--------------------
        for (int k = iz1 + 1; k < iz2; ++k)
        {
            real zz = (k-iz1)*hzw;
            real z1 = f3*sin(M_PI*zz/dw);
                u0(ixmin, iymin, k) = ((u0(ixmin, iymin, k - 1) + u0(ixmin, iymin, k + 1)) * hzhzw
                                       +      (u0(ixmax - 1, iymin, k) + u0(ixmin + 1, iymin, k)) * hxhx
                                       +      (u0(ixmin, iymax - 1, k) + u0(ixmin, iymin + 1, k)) * hyhy) * hxhyhzw;
                       real dE = efermi-Ez0(ixmin, iymin);
                       if(dE>0) u0(ixmin, iymin,k) = u0(ixmin, iymin,k) + z1 * dE;

                u1(ixmin, iymin, k) = u0(ixmin, iymin, k);
                u0(ixmin, iymax, k) = ((u0(ixmin, iymax, k - 1) + u0(ixmin, iymax, k + 1)) * hzhzw
                                       +      (u0(ixmax - 1, iymax, k) + u0(ixmin + 1, iymax, k)) * hxhx
                                       +      (u0(ixmin, iymax - 1, k) + u0(ixmin, iymin + 1, k)) * hyhy) * hxhyhzw;
                       dE = efermi-Ez0(ixmin, iymax);
                       if(dE>0) u0(ixmin, iymax,k) = u0(ixmin, iymax,k) + z1 * dE;
                u1(ixmin, iymax, k) = u0(ixmin, iymax, k);
                u0(ixmax, iymax, k) = ((u0(ixmax, iymax, k - 1) + u0(ixmax, iymax, k + 1)) * hzhzw
                                       +      (u0(ixmax - 1, iymax, k) + u0(ixmin + 1, iymax, k)) * hxhx
                                       +      (u0(ixmax, iymax - 1, k) + u0(ixmax, iymin + 1, k)) * hyhy) * hxhyhzw;
                       dE = efermi-Ez0(ixmax, iymax);
                       if(dE>0) u0(ixmax, iymax,k) = u0(ixmax, iymax,k) + z1 * dE;
                u1(ixmax, iymax, k) = u0(ixmax, iymax, k);
                u0(ixmax, iymin, k) = ((u0(ixmax, iymin, k - 1) + u0(ixmax, iymin, k + 1)) * hzhzw
                                       +      (u0(ixmax - 1, iymin, k) + u0(ixmin + 1, iymin, k)) * hxhx
                                       +      (u0(ixmax, iymax - 1, k) + u0(ixmax, iymin + 1, k)) * hyhy) * hxhyhzw;
                       dE = efermi-Ez0(ixmax, iymin);
                       if(dE>0) u0(ixmax, iymin,k) = u0(ixmax, iymin,k) + z1 * dE;
                u1(ixmax, iymin, k) = u0(ixmax, iymin, k);
        }
//----------------------
//сшивка k=iz2
        kk=iz2;
        j = iymin;
        i = ixmin;
        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaalas/hzb + u0(i, j, kk - 1)*eps_gaas/hzw)
                                       -      (u0(i, j, kk + 2)*eps_gaalas/hzb + u0(i, j, kk - 2)*eps_gaas/hzw))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hzb));
        j = iymax;
        i = ixmin;
        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaalas/hzb + u0(i, j, kk - 1)*eps_gaas/hzw)
                                       -      (u0(i, j, kk + 2)*eps_gaalas/hzb + u0(i, j, kk - 2)*eps_gaas/hzw))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hzb));
        j = iymax;
        i = ixmax;
        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaalas/hzb + u0(i, j, kk - 1)*eps_gaas/hzw)
                                       -      (u0(i, j, kk + 2)*eps_gaalas/hzb + u0(i, j, kk - 2)*eps_gaas/hzw))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hzb));
        j = iymin;
        i = ixmax;
        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaalas/hzb + u0(i, j, kk - 1)*eps_gaas/hzw)
                                       -      (u0(i, j, kk + 2)*eps_gaalas/hzb + u0(i, j, kk - 2)*eps_gaas/hzw))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hzb));
        u1(ixmax, iymin, kk) = u0(ixmax, iymin, kk);
        u1(ixmin, iymin, kk) = u0(ixmin, iymin, kk);
        u1(ixmax, iymax, kk) = u0(ixmax, iymax, kk);
        u1(ixmin, iymax, kk) = u0(ixmin, iymax, kk);

        for (int k = iz2 + 1; k < izb; ++k)
        {
                u0(ixmin, iymin, k) = ((u0(ixmin, iymin, k - 1) + u0(ixmin, iymin, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymin, k) + u0(ixmin + 1, iymin, k)) * hxhx
                                       +      (u0(ixmin, iymax - 1, k) + u0(ixmin, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmin, iymax, k) = ((u0(ixmin, iymax, k - 1) + u0(ixmin, iymax, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymax, k) + u0(ixmin + 1, iymax, k)) * hxhx
                                       +      (u0(ixmin, iymax - 1, k) + u0(ixmin, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmax, iymax, k) = ((u0(ixmax, iymax, k - 1) + u0(ixmax, iymax, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymax, k) + u0(ixmin + 1, iymax, k)) * hxhx
                                       +      (u0(ixmax, iymax - 1, k) + u0(ixmax, iymin + 1, k)) * hyhy) * hxhyhz;
                u0(ixmax, iymin, k) = ((u0(ixmax, iymin, k - 1) + u0(ixmax, iymin, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, iymin, k) + u0(ixmin + 1, iymin, k)) * hxhx
                                       +      (u0(ixmax, iymax - 1, k) + u0(ixmax, iymin + 1, k)) * hyhy) * hxhyhz;
        u1(ixmax, iymin, k) = u0(ixmax, iymin, k);
        u1(ixmin, iymin, k) = u0(ixmin, iymin, k);
        u1(ixmax, iymax, k) = u0(ixmax, iymax, k);
        u1(ixmin, iymax, k) = u0(ixmin, iymax, k);
        }

//2222
//        OMP(parallel for)
        for (int k = 1; k < iz_ins2; ++k)
        {
//#pragma ivdep
//#pragma simd
            for (int j = iymin+1; j < iymax; ++j)
            {
                u0(ixmin, j, k) = ((u0(ixmin, j, k - 1) + u0(ixmin, j, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, j, k) + u0(ixmin + 1, j, k)) * hxhx
                                       +      (u0(ixmin, j - 1, k) + u0(ixmin, j + 1, k)) * hyhy) * hxhyhz;

                u1(ixmax, j, k) = u0(ixmin, j, k);
                u1(ixmax, j, k) = u0(ixmin, j, k);
                u1(ixmin, j, k) = u0(ixmin, j, k);
            }
        }
//        OMP(parallel for)
        for (int k = 1; k < iz_ins2; ++k)
        {
//#pragma ivdep
//#pragma simd
            for (int i = ixmin+1; i < ixmax; ++i)
            {
                u0(i, iymin, k) = ((u0(i, iymin, k - 1) + u0(i, iymin, k + 1)) * hzhz
                                       +      (u0(i - 1, iymin, k) + u0(i + 1, iymin, k)) * hxhx
                                       +      (u0(i, iymax - 1, k) + u0(i, iymin + 1, k)) * hyhy) * hxhyhz;

                u0(i, iymax, k) = u0(i, iymin, k);
                u1(i, iymin, k) = u0(i, iymin, k);
                u1(i, iymax, k) = u0(i, iymin, k);
            }
        }
//!!!!!!!!!!!!!!!!!!!!!!! раскомментировала
        kk=iz_ins2;
        for (int j = iymin+1; j < iymax; ++j)
        {
                       int i=ixmin;
                       real x=hx*ixmin+xmin;
                       real y=hy*j+ymin;
//                       real dU1=ns*(0.1*1.602/8.85);
                       real dU1=-1000*s*(aimp1(ixmin,j)+aimp1(ixmax,j))*(0.1*1.602/8.85);
//                       real dU1=-1000*s*0.5*(aimp1(ixmin,j)+aimp1(ixmax,j))*(0.1*1.602/8.85);
                      if(!splitGate(x,y))
                        {
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hz + u0(i, j, kk - 1)*eps_ins/hz_gate)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hz + u0(i, j, kk - 2)*eps_ins/hz_gate)-2.*dU1)/
                                               (3.*(eps_gaas/hz+eps_ins/hz_gate));

                        }
                       else u0(i, j, kk) = Usg;
                       u1(ixmin, j, kk) = u0(ixmin, j, kk);
                       u0(ixmax, j, kk) = u0(ixmin, j, kk);
                       u1(ixmax, j, kk) = u0(ixmax, j, kk);
          }
        for (int i = ixmin+1; i < ixmax; ++i)
        {
                       int j=iymin;
                       real x=hx*i+xmin;
                       real y=hy*j+ymin;
//                       real dU1=ns*(0.1*1.602/8.85);
                       real dU1=-1000*s*(aimp1(i,iymin)+aimp1(i,iymax))*(0.1*1.602/8.85);
                       if(!splitGate(x,y))
                        {
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hz + u0(i, j, kk - 1)*eps_ins/hz_gate)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hz + u0(i, j, kk - 2)*eps_ins/hz_gate)-2.*dU1)/
                                               (3.*(eps_gaas/hz+eps_ins/hz_gate));
                        }
                       else u0(i, j, kk) = Usg;
                       u1(i, j, kk) = u0(i, j, kk);
                       u0(i, iymax, kk) = u0(i, j, kk);
                       u1(i, iymax, kk) = u0(i, j, kk);
         }
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//        OMP(parallel for)
        for (int k = iz_ins2+1; k < iz1; ++k)
        {
//#pragma ivdep
//#pragma simd
            for (int i = ixmin+1; i < ixmax; ++i)
            {
                u0(i, iymin, k) = ((u0(i, iymin, k - 1) + u0(i, iymin, k + 1)) * hzhz
                                       +      (u0(i - 1, iymin, k) + u0(i + 1, iymin, k)) * hxhx
                                       +      (u0(i, iymax - 1, k) + u0(i, iymin + 1, k)) * hyhy) * hxhyhz;

                u0(i, iymax, k) = u0(i, iymin, k);
                u1(i, iymin, k) = u0(i, iymin, k);
                u1(i, iymax, k) = u0(i, iymin, k);

            }
        }

//        OMP(parallel for)
        for (int k = iz_ins2+1; k < iz1; ++k)
        {
//#pragma ivdep
//#pragma simd
            for (int j = iymin+1; j < iymax; ++j)
            {
                u0(ixmin, j, k) = ((u0(ixmin, j, k - 1) + u0(ixmin, j, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, j, k) + u0(ixmin + 1, j, k)) * hxhx
                                       +      (u0(ixmin, j - 1, k) + u0(ixmin, j + 1, k)) * hyhy) * hxhyhz;
                u0(ixmax, j, k) = u0(ixmin, j, k);
                u1(ixmax, j, k) = u0(ixmin, j, k);
                u1(ixmin, j, k) = u0(ixmin, j, k);
        }
        }
//сшивка k=iz1   +
        kk=iz1;
          for (int j=iymax-1; j>iymin; --j)
        {
            int i = ixmin;
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzw + u0(i, j, kk - 1)*eps_gaalas/hz)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzw + u0(i, j, kk - 2)*eps_gaalas/hz))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hz));
                        u0(ixmax, j, kk) = u0(ixmin, j, kk);
                        u0(ixmax, j, kk) = u0(ixmin, j, kk);
                        u1(ixmin, j, kk) = u0(ixmin, j, kk);
        }
          for (int i=ixmax-1; i>ixmin; --i)
        {
            int j = iymin;
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzw + u0(i, j, kk - 1)*eps_gaalas/hz)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzw + u0(i, j, kk - 2)*eps_gaalas/hz))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hz));
                        u0(i, iymax, kk) = u0(i, iymin, kk);
                        u1(i, iymax, kk) = u0(i, iymin, kk);
                        u1(i, iymin, kk) = u0(i, iymin, kk);
        }
/////////////////////////
//        OMP(parallel for)
        for (int k = iz1 + 1; k < iz2; ++k)
        {
            real zz = (k-iz1)*hzw;
            real z1 = f3*sin(M_PI*zz/dw);
//#pragma ivdep
//#pragma simd
            for (int j = iymin+1; j < iymax; ++j)
            {
                u0(ixmin, j, k) = ((u0(ixmin, j, k - 1) + u0(ixmin, j, k + 1)) * hzhzw
                                       +      (u0(ixmax - 1, j, k) + u0(ixmin + 1, j, k)) * hxhx
                                       +      (u0(ixmin, j - 1, k) + u0(ixmin, j + 1, k)) * hyhy) * hxhyhzw;
                       real dE = efermi-Ez0(ixmin,j);
                       if(dE>0) u0(ixmin,j,k) = u0(ixmin,j,k) + z1 * dE;
                       u1(ixmin, j, k) = u0(ixmin, j, k);
                u0(ixmax, j, k) = ((u0(ixmax, j, k - 1) + u0(ixmax, j, k + 1)) * hzhzw
                                       +      (u0(ixmax - 1, j, k) + u0(ixmin + 1, j, k)) * hxhx
                                       +      (u0(ixmax, j - 1, k) + u0(ixmax, j + 1, k)) * hyhy) * hxhyhzw;
                       dE = efermi-Ez0(ixmax,j);
                       if(dE>0) u0(ixmax,j,k) = u0(ixmax,j,k) + z1 * dE;
                u1(ixmax, j, k) = u0(ixmax, j, k);
            }
        }
//        OMP(parallel for)
        for (int k = iz1 + 1; k < iz2; ++k)
        {
            real zz = (k-iz1)*hzw;
            real z1 = f3*sin(M_PI*zz/dw);
//#pragma ivdep
//#pragma simd
            for (int i = ixmin+1; i < ixmax; ++i)
            {
                u0(i, iymin, k) = ((u0(i, iymin, k - 1) + u0(i, iymin, k + 1)) * hzhzw
                                       +      (u0(i - 1, iymin, k) + u0(i + 1, iymin, k)) * hxhx
                                       +      (u0(i, iymax - 1, k) + u0(i, iymin + 1, k)) * hyhy) * hxhyhzw;
                       real dE = efermi-Ez0(i,iymin);
                       if(dE>0) u0(i,iymin,k) = u0(i,iymin,k) + z1 * dE;
                u1(i, iymin, k) = u0(i, iymin, k);
                u0(i, iymax, k) = ((u0(i, iymax, k - 1) + u0(i, iymax, k + 1)) * hzhzw
                                       +      (u0(i - 1, iymax, k) + u0(i + 1, iymax, k)) * hxhx
                                       +      (u0(i, iymax - 1, k) + u0(i, iymin + 1, k)) * hyhy) * hxhyhzw;
                       dE = efermi-Ez0(i,iymax);
                       if(dE>0) u0(i,iymax,k) = u0(i,iymax,k) + z1 * dE;
                u1(i, iymax, k) = u0(i, iymax, k);
            }
        }

//////////////////////////
//сшивка k=iz2   +
        kk=iz2;
          for (int j=iymax-1; j>iymin; --j)
        {
            int i = ixmin;
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzb + u0(i, j, kk - 1)*eps_gaas/hzw)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzb + u0(i, j, kk - 2)*eps_gaas/hzw))/
                                               (3.*(eps_gaas/hzw+eps_gaas/hzb));
                        u0(ixmax, j, kk) = u0(ixmin, j, kk);
                        u1(ixmax, j, kk) = u0(ixmin, j, kk);
                        u1(ixmin, j, kk) = u0(ixmin, j, kk);
        }
        for (int i = ixmin+1; i < ixmax; ++i)
        {
            int j = iymin;
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzb + u0(i, j, kk - 1)*eps_gaas/hzw)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzb + u0(i, j, kk - 2)*eps_gaas/hzw))/
                                               (3.*(eps_gaas/hzw+eps_gaas/hzb));
                        u0(i, iymax, kk) = u0(i, iymin, kk);
                        u1(i, iymax, kk) = u0(i, iymin, kk);
                        u1(i, iymin, kk) = u0(i, iymin, kk);
        }
//        OMP(parallel for)
        for (int k = iz2 + 1; k < izb; ++k)
        {
//#pragma ivdep
//#pragma simd
            for (int j = iymin+1; j < iymax; ++j)
            {
                u0(ixmin, j, k) = ((u0(ixmin, j, k - 1) + u0(ixmin, j, k + 1)) * hzhz
                                       +      (u0(ixmax - 1, j, k) + u0(ixmin + 1, j, k)) * hxhx
                                       +      (u0(ixmin, j - 1, k) + u0(ixmin, j + 1, k)) * hyhy) * hxhyhz;

                u1(ixmax, j, k) = u0(ixmin, j, k);
                u1(ixmin, j, k) = u0(ixmin, j, k);
                u0(ixmax, j, k) = u0(ixmin, j, k);
            }
        }
//        OMP(parallel for)
        for (int k = iz2 + 1; k < izb; ++k)
        {
//#pragma ivdep
//#pragma simd
            for (int i = ixmin+1; i < ixmax; ++i)
            {
                u0(i, iymin, k) = ((u0(i, iymin, k - 1) + u0(i, iymin, k + 1)) * hzhz
                                       +      (u0(i - 1, iymin, k) + u0(i + 1, iymin, k)) * hxhx
                                       +      (u0(i, iymax - 1, k) + u0(i, iymin + 1, k)) * hyhy) * hxhyhz;

                u0(i, iymax, k) = u0(i, iymin, k);
                u1(i, iymin, k) = u0(i, iymin, k);
                u1(i, iymax, k) = u0(i, iymin, k);
            }
        }
// j=iymin j=iymax
        dt[1] = _rdtsc()-dt[1];

        dt[2] = _rdtsc();
//----first layer--insolator
    hzhz = 1. / (hz_ins * hz_ins);
    hxhyhz = 1. / (2. * hzhz + 2. * hxhx + 2. * hyhy);
        for (int k = 1; k < iz_ins1; ++k)
        {
//            OMP(parallel for)
            for (int j = iymin + 1; j < iymax; ++j)
            {
                real y=hy*j+ymin;
//#pragma ivdep
//#pragma simd
                for (int i = ixmin + 1; i < ixmax; ++i)
                {
                    real x=hx*i+xmin;
                    if(topGate(x,y,t_ins)) {u1(i, j, k) = Utg; u0(i, j, k) = Utg;}
                    else
                    {
                        u1(i, j, k) = ((u0(i, j, k - 1) + u0(i, j, k + 1)) * hzhz
                                       +      (u0(i - 1, j, k) + u0(i + 1, j, k)) * hxhx
                                       +      (u0(i, j - 1, k) + u0(i, j + 1, k)) * hyhy) * hxhyhz;
                     }
                }
            }
        }

//----first layer--insolator
    hzhz = 1. / (hz_ins * hz_ins);
    hxhyhz = 1. / (2. * hzhz + 2. * hxhx + 2. * hyhy);
        for (int k = iz_ins1; k < iz_ins2; ++k)
        {
//            OMP(parallel for)
            for (int j = iymin + 1; j < iymax; ++j)
            {
                real y=hy*j+ymin;
//#pragma ivdep
//#pragma simd
                for (int i = ixmin + 1; i < ixmax; ++i)
                {
                    real x=hx*i+xmin;
                    if(topGate(x,y,t_ins)&&k<=iz_ins2-iz_ins1) {u1(i, j, k) = Utg; u0(i, j, k) = Utg;}
                    else
                    {
                        u1(i, j, k) = ((u0(i, j, k - 1) + u0(i, j, k + 1)) * hzhz
                                       +      (u0(i - 1, j, k) + u0(i + 1, j, k)) * hxhx
                                       +      (u0(i, j - 1, k) + u0(i, j + 1, k)) * hyhy) * hxhyhz;
                    }
                    if(splitGate(x,y)) {u1(i, j, k) = Usg; u0(i, j, k) = Usg;}
//                 u1(i, j, kk) = u0(i, j, kk);
                }
            }
        }
        dt[2] = _rdtsc()-dt[2];
// correct derivatives for between insolator and GaAs layers
        dt[3] = _rdtsc();
            kk=iz_ins2;
//            OMP(parallel for)
            for (int j = 1; j <iymax; ++j)
            {
                real y=hy*j+ymin;
//#pragma ivdep
//#pragma simd
                for (int i = 1; i < ixmax; ++i)
                {
                       real x=hx*i+xmin;
//                       real dU1=ns*(0.1*1.602/8.85);
                       real dU1=-1000*s*aimp1(i,j)*(0.1*1.602/8.85);
//                       real dU1=(1000*s*aimp1(i,j)+ns)*(0.1*1.602/8.85);
                      if(!splitGate(x,y))
                        {
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hz + u0(i, j, kk - 1)*eps_ins/hz_gate)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hz + u0(i, j, kk - 2)*eps_ins/hz_gate)-2.*dU1)/
                                               (3.*(eps_gaas/hz+eps_ins/hz_gate));
                        }
                        else u0(i, j, kk) = Usg;
                        u1(i, j, kk) = u0(i, j, kk);
                }
            }
// calculate U_new for cap layer
    hzhz = 1. / (hz * hz);
    hxhyhz = 1. / (2. * hzhz + 2. * hxhx + 2. * hyhy);
        for (int k = iz_ins2+1; k < izc; ++k)
        {
            OMP(parallel for)
            for (int j = iymin + 1; j < iymax; ++j)
            {
#pragma ivdep
#pragma simd
                for (int i = ixmin + 1; i < ixmax; ++i)
                {
                        u1(i, j, k) = ((u0(i, j, k - 1) + u0(i, j, k + 1)) * hzhz
                                       +      (u0(i - 1, j, k) + u0(i + 1, j, k)) * hxhx
                                       +      (u0(i, j - 1, k) + u0(i, j + 1, k)) * hyhy) * hxhyhz;
                }
            }
        }
        dt[3] = _rdtsc()-dt[3];
        dt[4] = _rdtsc();
// correct derivatives for between cap and spacer AlGaAs layers
        kk=izc;
        OMP(parallel for)
        for (int j=0; j<=iymax; ++j)
        {
#pragma ivdep
#pragma simd
            for (int i = ixmin+1; i < ixmax; ++i)
            {
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaalas/hz + u0(i, j, kk - 1)*eps_gaas/hz)
                                       -      (u0(i, j, kk + 2)*eps_gaalas/hz + u0(i, j, kk - 2)*eps_gaas/hz))/
                                               (3.*(eps_gaas/hz+eps_gaalas/hz));
                        u1(i, j, kk) = u0(i, j, kk);
            }
        }

// calculate U_new for AlGaAs layer
        hzhz = 1. / (hz*hz);
        hxhyhz = 1. / (2.*hzhz+2.*hxhx+2.*hyhy);
        for (int k = izc+1; k < iz1; ++k)
        {
            OMP(parallel for)
            for (int j = iymin + 1; j < iymax; ++j)
            {
#pragma ivdep
#pragma simd
                for (int i = ixmin + 1; i < ixmax; ++i)
                {
                        u1(i, j, k) = ((u0(i, j, k - 1) + u0(i, j, k + 1)) * hzhz
                                       +      (u0(i - 1, j, k) + u0(i + 1, j, k)) * hxhx
                                       +      (u0(i, j - 1, k) + u0(i, j + 1, k)) * hyhy) * hxhyhz;
                }
            }
        }
        dt[4] = _rdtsc()-dt[4];

// correct derivatives for between spacer and 2deg layers
        dt[5] = _rdtsc();
        kk=iz1;
        OMP(parallel for)
        for (int j=iymin+1; j<iymax; ++j)
        {
#pragma ivdep
#pragma simd
            for (int i = ixmin+1; i < ixmax; ++i)
            {
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaas/hzw + u0(i, j, kk - 1)*eps_gaalas/hz)
                                       -      (u0(i, j, kk + 2)*eps_gaas/hzw + u0(i, j, kk - 2)*eps_gaalas/hz))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hz));
                        u1(i, j, kk) = u0(i, j, kk);
            }
        }

        // 2deg layer
        for (int k = iz1 + 1; k < iz2; ++k)
        {
            real zz = (k-iz1)*hzw;
            real z1 = f3*sin(M_PI*zz/dw);
            OMP(parallel for)
            for (int j = iymin + 1; j < iymax; ++j)
            {
#pragma ivdep
#pragma simd
                for (int i = ixmin + 1; i < ixmax; ++i)
                {
                    u1(i, j, k) = ((u0(i, j, k - 1) + u0(i, j, k + 1)) * hzhzw
                                   +      (u0(i - 1, j, k) + u0(i + 1, j, k)) * hxhx
                                   +      (u0(i, j - 1, k) + u0(i, j + 1, k)) * hyhy) * hxhyhzw;
                       real dE = efermi-Ez0(i,j);
                       if(dE>0) u1(i,j,k) = u1(i,j,k) + z1 * dE;
                }
            }
        }
        kk=iz2;
        OMP(parallel for)
        for (int j=iymin+1; j<iymax; ++j)
        {
#pragma ivdep
#pragma simd
            for (int i = ixmin+1; i < ixmax; ++i)
            {
                        u0(i, j, kk) = (4.*(u0(i, j, kk + 1)*eps_gaalas/hzb + u0(i, j, kk - 1)*eps_gaas/hzw)
                                       -      (u0(i, j, kk + 2)*eps_gaalas/hzb + u0(i, j, kk - 2)*eps_gaas/hzw))/
                                               (3.*(eps_gaas/hzw+eps_gaalas/hzb));
                        u1(i, j, kk) = u0(i, j, kk);
            }
        }
        dt[5] = _rdtsc()-dt[5];
        // bulk layer
        dt[6] = _rdtsc();
        hzhz = 1. / (hzb*hzb);
        hxhyhz = 1. / (2.*hzhz+2.*hxhx+2.*hyhy);
        for (int k = iz2 + 1; k < izb; ++k)
        {
            OMP(parallel for)
            for (int j = iymin + 1; j < iymax; ++j)
            {
#pragma ivdep
#pragma simd
                for (int i = ixmin + 1; i < ixmax; ++i)
                {
                    u1(i, j, k) = ((u0(i, j, k - 1) + u0(i, j, k + 1)) * hzhz
                                   +      (u0(i - 1, j, k) + u0(i + 1, j, k)) * hxhx
                                   +      (u0(i, j - 1, k) + u0(i, j + 1, k)) * hyhy) * hxhyhz;
                }
             }
         }
        dt[6] = _rdtsc()-dt[6];
        { real *tmp = u0; u0 = u1; u1 = tmp; }
        //for (int i=0;i<10;++i) printf(ZI" ",dt[i]); printf("\n"); fflush(0);
    }
    *nflops = 8.0 * NX*NY*NZ * nrounds;
    return 0;
}
