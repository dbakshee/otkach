// -*- c -*-
// By Olga Tkachenko.
// Main program.
// structure 713
#include <stdio.h>
#include <malloc.h>
#include <xmmintrin.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#define malloc(sz) _mm_malloc(sz,64)
#define free(p) _mm_free(p)

#include "ringnew2.h"
#include "mytime.h"
#include <omp.h>

#if defined(_WIN32)
#define ZI "%Ii"
#else
#define ZI "%zi"
#endif

/*
#include <stdio.h>
#include <malloc.h>
#include <xmmintrin.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#define malloc(sz) _mm_malloc(sz,64)
#define free(p) _mm_free(p)

#include "ringnew2.h"
#include "mytime.h"
#include <omp.h>

#if defined(_WIN32)
#define ZI "%Ii"
#else
#define ZI "%zi"
#endif
#define ZZ "/home2/ssd1/otkach"
*/
void setimp(int nimpur, real *aimp1, int ixa, int ixb, int iya, int iyb)
{
    int ixmin = ixa;
    int iymin = iya;
    int ixmax = ixb;
    //first clear array
    for (int j=iya; j<=iyb; ++j)
	for (int i=ixa; i<=ixb; ++i)
	    {
//		aimp1(i,j) = 0.001*7.5*25;//0.0;
		aimp1(i,j) = 0;//0.001*7.5*25;//0.0;
	    }
    //fill the array

    int k = 0;
    while ( k < nimpur )
	{
	    int i = int(0.5 + (double)rand()/RAND_MAX * (ixb - ixa) + ixa);
	    if (i > ixb) continue;
	    int j = int(0.5 + (double)rand()/RAND_MAX * (iyb - iya) + iya);
	    if (j > iyb) continue;
	    aimp1(i,j) += 1.;
	    ++k;
	}
}

void writeimp(const char *fname, real *aimp1, int ixa, int ixb, int iya, int iyb)
{
    int ixmin = ixa;
    int iymin = iya;
    int ixmax = ixb;
    FILE *o = fopen(fname, "wt");
    assert(o);
    for (int j=iya; j<=iyb; ++j)
	for (int i=ixa; i<=ixb; ++i)
	    {
		fprintf(o,"%i %i %lg\n", i, j, aimp1(i,j));
	    }
    fclose(o);
}
void readimp(const char *fname, real *aimp1, int Np, int ixa, int ixb, int iya, int iyb)
{
    int ixmin = ixa;
    int iymin = iya;
    int ixmax = ixb;
    FILE *o = fopen(fname, "rt");
    assert(o);
    int lineno = 0;
    for (int j=iya; j<=iyb; ++j){
	for (int i=ixa; i<=ixb; ++i){
		aimp1(i,j) = 0.0;}}
	for (int n=0; n<Np; ++n)
	    {
            int ix,iy;
            assert(2 == fscanf(o, "%i %i\n", &ix, &iy));
            if (ix != ix||iy != iy)
            {
               printf("OOpps on line %i i have found a NaN\n",lineno);
               exit(1);
            }
            lineno++;
            aimp1(ix,iy)+=1;
//            printf("line=%i %i %i %lg\n",lineno,ix,iy,aimp1(ix,iy));
	    }
    fclose(o);
    printf("Reading1 is done!\n");
}

// real *u, ***v;
// #define u(x,y,z)  u[....]
// v = map3d(u, ... )
//
// u is addressed thus: u(x,y,z)
// v is addressed thus: v[z][y][x]
real ***
make_map3d (real * u, int X, int Y, int Z)
{
  //real *u = (real*  )malloc(sizeof(real  ) * Z*Y*X );
  real **t = (real **) malloc (sizeof (real *) * Z * Y);
  real ***v = (real ***) malloc (sizeof (real **) * Z);
  for (int z = 0; z < Z; ++z)
  {
    v[z] = &t[z * Y];
    for (int y = 0; y < Y; ++y)
    {
      v[z][y] = &u[y * X + z * Y * X];
    }
  }
  return v;
}

void
free_map3d (real *** v)
{
  free (v[0]);
  free (v);
}

real **
make_map2d (real * u, int X, int Y)
{
  //real *u = (real*  )malloc(sizeof(real  ) * Y*X );
  real **v = (real **) malloc (sizeof (real *) * Y);
  for (int y = 0; y < Y; ++y)
  {
    v[y] = &u[0 + y * X];
  }
  return v;
}

void
free_map2d (real ** v)
{
  free (v);
}

#define HBEAT 5000.0		/* milliseconds */
void
write3_to_file (const char *fname, real * u0, int i, int j, struct Params *params)
{
  real hx = params->hx;
  real hy = params->hy;
  real hz = params->hz;
  real hz_ins = params->hz_ins;
  real hz_gate = params->hz_gate;
  real hzw = params->hzw;
  real hzb = params->hzb;
  int iz_ins1 = params->iz_ins1;
  int iz_ins2 = params->iz_ins2;
  int izb = params->izb;
  int izc = params->izc;
  int iz1 = params->iz1;
  int iz2 = params->iz2;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int iymin = params->iymin;
  int iymax = params->iymax;

  FILE *o = fopen (fname, "wt");
  //  assert(o);
  fprintf (o, "Ueff along z at x_i=%.8g y_j=%.8g\n", i * hx, j * hy);
  real z = 0;
  for (int k = 0; k < iz_ins1 + 1; ++k)
  {
    z = hz_ins * k;
    fprintf (o, "%.8g %.8g\n", z, u0 (i, j, k));
  }
  real z1 = z;
  for (int k = iz_ins1 + 1; k <= iz_ins2; ++k)
  {
    z = z1 + hz_gate * (k - iz_ins1);
    fprintf (o, "%.8g %.8g\n", z, u0 (i, j, k));
  }
  z1 = z;
  for (int k = iz_ins2 + 1; k <= izc; ++k)
  {
    z = z1 + hz * (k - iz_ins2);
    fprintf (o, "%.8g %.8g\n", z, u0 (i, j, k));
  }
  z1 = z;
  for (int k = izc; k <= iz1; ++k)
  {
    z = z1 + hz * (k - izc);
    fprintf (o, "%.8g %.8g\n", z, u0 (i, j, k)+0.152);
  }
  z1 = z;
  for (int k = iz1; k <= iz2; ++k)
  {
    z = z1 + hzw * (k - iz1);
    fprintf (o, "%.8g %.8g\n", z, u0 (i, j, k));
  }
  z1 = z;
  for (int k = iz2; k <= izb; ++k)
  {
    z = z1 + hzb * (k - iz2);
    fprintf (o, "%.8g %.8g\n", z, u0 (i, j, k)+0.152);
  }
  fclose (o);
}
/*
void
write2_to_file (const char *fname, real * Ez0, struct Params *params)
{
  FILE *o = fopen (fname, "wt");
  assert (o);
  int iymin = params->iymin;	//24;//iymin;//28;
  int iymax = params->iymax;	//109;//iymax;//105;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int NX = params->NX;
  int NY = params->NY;
  real hx = params->hx;
  real hy = params->hy;
  int iy0 = 1;			//3*50;//10*50;//10*17;//15*17;//68
  fprintf (o, "%3i  %.8g %.8g\n", ixmax - ixmin - 2 * iy0 + 1, 0., hx * (ixmax - ixmin - 2 * iy0));
  fprintf (o, "%3i  %.8g %.8g\n", iymax - iymin - 2 * iy0 + 1, 0., hy * (iymax - iymin - 2 * iy0));

  double n_av=0;
  double us=0;
  double coef = 324.*0.067/0.07766;
  int isum = 0;
  for (int i = ixmin+1; i < ixmax; ++i)
  {
    for (int j = iymin+1; j < iymax; ++j)
    {
      double ueff = -Ez0 (i, j);
      double  n=0;
      us = us+ueff;
      if(ueff>0) n=ueff*coef;
      n_av=n_av+n;
      isum = isum+1;
//      fprintf (o, "%.9lg\n", n);
//      fprintf (o, "%lg\n", -ueff);
    }
  }
      n_av=n_av/isum;
      us=us/isum;
  double rms = 0;
  for (int i = ixmin+1; i < ixmax; ++i)
  {
    for (int j = iymin+1; j < iymax; ++j)
    {
      double ueff = -Ez0 (i, j);
      double du   = ueff-us;
      rms = rms+du*du;
    }
   }
   rms = sqrt(rms/isum);
   fprintf(o,"%lg %lg %lg\n",n_av, us, rms);
  for (int j = iymin; j < iymax+1; ++j)
  {
    for (int i = ixmin; i < ixmax+1; ++i)
    {
      fprintf(o,"%lg ", Ez0 (i, j));
    }
        fprintf(o,"\n");
   }
  fclose (o);
}
*/

void
write2_to_file (const char *fname, real * Ez0, struct Params *params)
{
  FILE *o = fopen (fname, "wt");
  assert (o);
  int iymin = params->iymin;	//24;//iymin;//28;
  int iymax = params->iymax;	//109;//iymax;//105;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int NX = params->NX;
  int NY = params->NY;
  real hx = params->hx;
  real hy = params->hy;
  int iy0 = 1;			//3*50;//10*50;//10*17;//15*17;//68
  fprintf (o, "%3i  %.8g %.8g\n", ixmax - ixmin - 2 * iy0 + 1, 0., hx * (ixmax - ixmin - 2 * iy0));
  fprintf (o, "%3i  %.8g %.8g\n", iymax - iymin - 2 * iy0 + 1, 0., hy * (iymax - iymin - 2 * iy0));

  double n_av=0;
  double us=0;
  double coef = 324.*0.067/0.07766;
  int isum = 0;
  for (int i = ixmin+20; i <= ixmax-20; ++i)
  {
    for (int j = iymin+20; j < iymax-20; ++j)
    {
      double ueff = -Ez0 (i, j);
      double  n=0;
      us = us+ueff;
      if(ueff>0) n=ueff*coef;
      n_av=n_av+n;
      isum = isum+1;
//      fprintf (o, "%.9lg\n", n);
//      fprintf (o, "%lg\n", -ueff);
    }
  }
      n_av=n_av/isum;
      us=us/isum;
  double rms = 0;
  for (int i = ixmin+20; i < ixmax-20; ++i)
  {
    for (int j = iymin+20; j < iymax-20; ++j)
    {
      double ueff = -Ez0 (i, j);
      double du   = ueff-us;
      rms = rms+du*du;
    }
   }
   rms = sqrt(rms/isum);
   fprintf(o,"%lg %lg %lg\n",n_av, us, rms);
  for (int j = iymin; j <= iymax; ++j)
//  for (int j = iymin+20; j < iymax-20; ++j)
  {
//    for (int i = ixmin+20; i < ixmax-20; ++i)
    for (int i = ixmin; i <= ixmax; ++i)
    {
      fprintf(o,"%lg ", Ez0 (i, j));
    }
        fprintf(o,"\n");
   }
  fclose (o);
}
void
write2grad_to_file (const char *fname, real * Ez0, struct Params *params)
{
  FILE *o = fopen (fname, "wt");
  assert (o);
  int iymin = params->iymin;	//24;//iymin;//28;
  int iymax = params->iymax;	//109;//iymax;//105;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  real hx = params->hx;
  real hy = params->hy;
//      fprintf (o, "%lg %lg\n", hx, hy);

  int iy0 = 0;			//3*50;//10*50;//10*17;//15*17;//68
//  fprintf (o, "%3i  %.8g %.8g\n", ixmax - ixmin - 2 * iy0 + 1, 0., hx * (ixmax - ixmin - 2 * iy0));
//  fprintf (o, "%3i  %.8g %.8g\n", iymax - iymin - 2 * iy0 + 1, 0., hy * (iymax - iymin - 2 * iy0));
//  for (int i = ixmin + 1; i <= ixmax - 1; ++i)
//  {
  for (int j = iymin; j <= iymax; ++j)
    {  real y=hx*j;
  for (int i = ixmin; i <= ixmax; ++i)
  {
       real x=hx*i;
//            real ueff = Ez0(i, j);
//          fprintf(o,"%lg %lg %lg\n", hx*i, hy*j, ueff);
      real F, Fx,Fy, Uxp1, Uxm1, Uyp1, Uym1;
      if(i!=0&&j!=0&&i!=ixmax&&j!=iymax)
       {
       Fx=-0.5*1000*(Ez0 (i+1,  j)-Ez0 (i-1,  j))/hx;
       Fy=-0.5*1000*(Ez0 (i,  j+1)-Ez0 (i,  j-1))/hy;
       }
      else {
      if(i==0) Uxm1=Ez0(ixmax-1,j);
      if(j==0) Uym1=Ez0(i,ixmax-1);
      if(i==ixmax) Uxp1=Ez0(ixmin+1,j);
      if(j==iymax) Uyp1=Ez0(i,iymin+1);
       Fx=-0.5*(Uxp1-Uxm1)/hx;
       Fy=-0.5*(Uyp1-Uym1)/hy;
      }
       F = 1000*sqrt(Fx*Fx+Fy*Fy);
//          fprintf(o,"%lg %lg %lg\n", hx*(i-2*iy0), hy*(j-iy0), ueff);
//      fprintf (o, "%16.7g", F);
      fprintf (o, "%15.4g %15.4g %15.7g %15.7g\n", x,y,Fx,Fy);
    }
//      fprintf(o,"\n");
  }
  fclose (o);
}

void
write_Contour (const char *fname, real * Ez0, struct Params *params)
{
  FILE *o = fopen (fname, "wt");
  assert (o);
  int iymin = params->iymin;	//10;
  int iymax = params->iymax;	//134-10;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  for (int j = iymin; j <= iymax; ++j)
  {
    for (int i = ixmin; i <= ixmax; ++i)
    {
      real ueff = Ez0 (i, j);
      if (ueff < 0)
	fprintf (o, "%i %i\n", i, j);
    }
  }
  fclose (o);
}

void
writed_to_file (const char *fname, real * Ez0, int ix0, struct Params *params)
{
  real hx = params->hx;
  real hy = params->hy;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int iymin = params->iymin;
  int iymax = params->iymax;

  FILE *o = fopen (fname, "wt");
  assert (o);
  int ixs=10;
  int iys = 1;
  fprintf (o, "Ueff along diagonal starting at x_i=%.8g\n", hx * ixs);
  for (int i = 0; i <= ixmax - ixs; ++i)
  {
    if (i + ixs <= ixmax && i + iys <= iymax)
    {
      real ueff = Ez0 (i + ixs, i + iys);
      fprintf (o, "%.8g %.8g\n", 2 * hx * i, ueff);
    }
  }
  fclose (o);
}

void
writeEzofy_to_file(const char *fname, real * u0, int i, struct Params *params)
{
  real hx = params->hx;
  real hy = params->hy;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int iymin = params->iymin;
  int iymax = params->iymax;
  real hzw = params->hzw;
  int iz1 = params->iz1;
  int iz2 = params->iz2;

  FILE *o = fopen (fname, "wt");
  assert (o);
  double w=(iz2-iz1)*hzw;
    fprintf (o, "%.8g \n", w);
  fprintf (o, "Ez along y at x_i=%.8g\n", i * hx);
  for (int j = iymin; j <= iymax; ++j)
  {
    double Ez = (u0 (i, j, iz1) - u0 (i, j, iz2))/w;
    fprintf (o, "%.8g %.8g\n", hy * j, Ez);
  }
  fclose (o);
}

void
writex_to_file (const char *fname, real * Ez0, int i, struct Params *params)
{
  real hx = params->hx;
  real hy = params->hy;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int iymin = params->iymin;
  int iymax = params->iymax;
  real xmin = params->xmin;
  real ymin = params->ymin;

  FILE *o = fopen (fname, "wt");
  assert (o);
//  fprintf (o, "Ueff along y at x_i=%.8g\n", i * hx);
  double n=0;
  double coef = 324.*0.067/0.07766;
  for (int j = iymin; j <= iymax; ++j)
  {
    real ueff = Ez0 (i, j);
    if(ueff>0) n=ueff*coef;
    else n=0;
//    fprintf (o, "%.8g %.8g\n", hy * j, n);
    fprintf (o, "%.8g %.8g\n", ymin+hy * j, ueff);
  }
  fclose (o);
}

void
writey_to_file (const char *fname, real * Ez0, int i, struct Params *params)
{
  real hx = params->hx;
  real hy = params->hy;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int iymin = params->iymin;
  int iymax = params->iymax;
  real xmin = params->xmin;
  real ymin = params->ymin;

  FILE *o = fopen (fname, "wt");
  assert (o);
//  fprintf (o, "Ueff along x at y_i=%.8g\n", i * hy);
  double n=0;
  double coef = 324.*0.067/0.07766;   //542.364
  for (int j = ixmin; j <= ixmax; ++j)
  {
    real ueff = Ez0 (j, i);
    if(ueff>0) n=ueff*coef;
    else n=0;
    fprintf (o, "%.8g %.8g\n", xmin+hx * j, ueff);
//    fprintf (o, "%.8g %.8g\n", hx * j, n);
  }
  fclose (o);
}

void
read_dt (const char *fname, real * u0, struct Params *params)
{
//return;
  FILE *f = fopen (fname, "rb");
//    FILE *f = fopen("u_sushkov.dt","rb");//fopen(fname,"rb");
  if (!f)
    return;
  printf ("Reading dt file %s\n", fname);
//      printf("Reading dt file %s\n","u_sushkov.dt");
  size_t bytes;
  assert (1 == fread (&bytes, sizeof (bytes), 1, f));

  size_t iter;
  assert (1 == fread (&iter, sizeof (iter), 1, f));

  size_t n = (bytes - sizeof (iter)) / sizeof (real);
//    size_t n = (bytes-sizeof(iter))/sizeof(float);//!!!!
//    float *x = (float*)malloc(sizeof(float)*n);//!!!!
  real *x = (real *) malloc (sizeof (real) * n);
  assert (x);
//    size_t nread = fread( x, sizeof(float), n, f );
  size_t nread = fread (x, sizeof (real), n, f);
  if (n != nread)
  {
    printf ("OOps: wanted " ZI " numbers, but read only " ZI "\n", n, nread);
    exit (1);
  }
  fclose (f);

  int izb = params->izb;
  int iymin = params->iymin;
  int iymax = params->iymax;
  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int NX = params->NX;
  int NY = params->NY;
  int NZ = params->NZ;
  if (n == NX * NY * NZ)
  {
    printf ("file has been read: iter=" ZI " n=" ZI "\n", iter, n);
    OMP (parallel for)
      for (int k = 0; k <= izb; ++k)
      {
	for (int j = iymin; j <= iymax; ++j)
	{
	  for (int i = ixmin; i <= ixmax; ++i)
	  {
	    int ijk = ((i) - ixmin) + ((j) - iymin) * (ixmax - ixmin + 1) + (k) * (ixmax - ixmin + 1) * (iymax - iymin + 1);
	    u0 (i, j, k) = x[ijk];
	  }
	}
      }
  }
  free (x);
}

static void initParams (struct Params *);
void
damn_compute (const char *dtfilename)
{
  struct Params params0;
  initParams (&params0);

  struct Params *params = &params0;

  int ixmin = params->ixmin;
  int ixmax = params->ixmax;
  int iymin = params->iymin;
  int iymax = params->iymax;
  int iz1 = params->iz1;
  int iz2 = params->iz2;
  int izc = params->izc;
  int iz_ins1 = params->iz_ins1;
  int iz_ins2 = params->iz_ins2;
  int izb = params->izb;
  int NX = params->NX;
  int NY = params->NY;
  int NZ = params->NZ;
  real hx = params->hx;
  real hy = params->hy;
  printf ("allocate memory for %ix%ix%i mesh\n", NX, NY, NZ);
  real *u0 = (real *) malloc (sizeof (real) * NX * NY * NZ);
  real *u1 = (real *) malloc (sizeof (real) * NX * NY * NZ);
  real *Ez0 = (real *) malloc (sizeof (real) * NX * NY);
  real *aimp1 = (real *) malloc (sizeof (real) * NX * NY);
  real *Ez1 = (real *) malloc (sizeof (real) * NX * NY);
  printf ("u0=%p u1=%p\n", u0, u1);

  printf ("initial approximation\n");
  memset (u0, 0, sizeof (real) * NX * NY * NZ);
  memset (u1, 0, sizeof (real) * NX * NY * NZ);
  memset (Ez0, 0, sizeof (real) * NX * NY);
  memset (aimp1, 0, sizeof (real) * NX * NY);
//    if (dtfilename)
  if (1)
  {
      read_dt("pot_n_imp7i5_500Kn.dt",u0,params); //  uniform 50K

      double us=0;
      for (int j = iymin; j <= iymax; ++j)
      {
	for (int i = ixmin; i <= ixmax; ++i)
	{
	  Ez0(i,j)=u0 (i, j, iz1+8 )+0.0219;
          us = us+Ez0(i,j);
	}
      }
      us=us/NX/NY;
      printf ("us=%13.6g\n", us);
      write2_to_file ("ueff_500K.dat", Ez0, &params0);
      int ix0 = 0;
      writex_to_file ("Uofy_500K_c.dat", Ez0, (NX+1)/2, &params0);
      writex_to_file ("Uofy_500K_0.dat", Ez0, ix0, &params0);
      writex_to_file ("Uofy_500K_xmax.dat", Ez0, ixmax, &params0);
      write3_to_file ("Uofz_uni_2deg_500K.dat", u1, 5, 100, &params0);
      write3_to_file ("Uofz_uni_gate_500K.dat", u1, 100, 180, &params0);
      int iy0 = 0;
      writey_to_file ("Uofx_500K_c.dat", Ez0, (NY+1)/2, &params0);
      writey_to_file ("Uofx_500K_0.dat", Ez0, iy0, &params0);
      writey_to_file ("Uofx_500K_ymax.dat", Ez0, iymax, &params0);

    //  exit(0);
  }
  else
  {
    {
      for (int j = iymin; j <= iymax; ++j)
      {
	for (int i = ixmin; i <= ixmax; ++i)
	{
	  Ez1 (i, j) = u0 (i, j, iz2);
	}
      }
    }
    exit (0);
  }

  // And now we create splitGate - array of 0/1 such that 1 means
  // 'gate' and 0 means 'hole'
  //    int *splitGate = create_splitGate(....);


#ifdef _OPENMP
  printf ("computing with simple iterations on %i threads...\n", omp_get_max_threads ());
#else
  printf ("computing with simple iterations ON 1 THREAD!...\n");
#endif

  /* Create 3d-map u0_3 for u0 (this is for MINT only) */
  real ***u0_3 = make_map3d (u0, NX, NY, NZ);
  real ***u1_3 = make_map3d (u1, NX, NY, NZ);
  real **Ez0_2 = make_map2d (Ez0, NX, NY);

  mytime_t myt0 = mytime ();	/* myt0 is the time when we started computation */
  int itn = 1;
  int itlim = 1000;
  int hbeat = 3;//50; /* we will adjust hbeat to report every HBEAT milliseconds */
  int imp1seed=37;//11;//17;//11;//71117;//11;//71117;
  srand(imp1seed);
  real nimp=7.5;//3;//1;//3;//3;//1;
  int nimpur=45114;//0.001*nimp*(ixmax-ixmin)*(iymax-iymin)*hx*hy;
//  setimp(abs(nimpur),aimp1,ixmin,ixmax,iymin,iymax);
//impn7i5_500K.dat
//impn7i5_500K_2deg.dat
  readimp("imp_500K.dat", aimp1, nimpur, ixmin, ixmax, iymin, iymax);
//  readimp("impn7i5_500K_2deg.dat", aimp1, nimpur, ixmin, ixmax, iymin, iymax);
//  int nimpur=4080;//0.001*nimp*(ixmax-ixmin)*(iymax-iymin)*hx*hy;
//  readimp("impn1_5K_new.dat", aimp1, nimpur, ixmin, ixmax, iymin, iymax);
//  writeimp("imp_500Ka.dat", aimp1, ixmin, ixmax, iymin, iymax);
  for (int it = 0; it < 4001000;)
  {
    /* Compute a block of hbeat iterations */
    mytime_t myt1 = mytime ();

    double nflops;
//#if defined(USE_CUDA)
//    compute_cuda (u0_3, u1_3, hbeat, &nflops, params);
//#else
    compute_host (u0, u1, Ez0, aimp1, hbeat, &nflops, params);
//#endif

    mytime_t myt2 = mytime ();

    /* Compute various things to report */

    /* Report status and adjust hbeat */
    double ms = mytime_ms (myt1, myt2);
    double s0 = mytime_s (myt0, myt2);
    real speriter = ms * 1e-3 / hbeat;
    real gf = nflops / ms * 1e-6;

    it += hbeat;
    hbeat *= HBEAT / ms;
    hbeat = (hbeat & (~1));
    if (hbeat < 2)
      hbeat = 100;//50;		//2;

    if (1)			/* write some potential profiles */
    {
      int ix0 = (NX+1)/2;
      writex_to_file ("Uofy_500K.dat", Ez0, ix0, &params0);
      writex_to_file ("Uofy_0_500K.dat", Ez0, 0, &params0);
      write3_to_file ("Uofz_2deg_500K.dat", u1, 5, 100, &params0);
      write3_to_file ("Uofz_gate_500K.dat", u1, 100, 180, &params0);
      int iy0 = (NY+1)/2;
      writey_to_file ("Uofx_500K.dat", Ez0, iy0, &params0);
      writey_to_file ("Uofx_0_500K.dat", Ez0, 0, &params0);
    }

    /* look for max difference between u0 and u1 */
    real error = 0;
    int ier = 0;
    int jer = 0;
    int ker = 0;
    {
      OMP (parallel for)
	for (int k = 1; k < izb; ++k)
	{
	  for (int j = iymin + 1; j < iymax; ++j)
	  {
	    for (int i = ixmin + 1; i < ixmax; ++i)
	    {
	      real errtest = u0 (i, j, k) - u1 (i, j, k);
	      if (fabs (errtest) > fabs (error))
	      {
		OMP (critical)
		{
		  error = errtest;
//          if (fabs(errtest) > fabs(error)) error = errtest;
//  ier=i; jer=j; ker=k;
		}
	      }
	    }
	  }
	}
     }
      printf (" %5i(%3i) %4.0fs err=%13.6g  %10.5gs/iter %10.3gGF\n", it, hbeat, s0, error, speriter, gf);
      fflush(0);


    if (it>itlim)			/* save breakpoint file dt */
    {
      write2_to_file ("ueff_500K.dat", Ez0, &params0);
      itn=itn+1;
      itlim = 1000*itn;
      size_t iter = it;
      FILE *f = fopen ("pot_n_imp7i5_500Kpar.dt", "wb");
      assert (f);
      size_t bytes = sizeof (iter) + sizeof (real) * NX * NY * NZ;
      fwrite (&bytes, sizeof (bytes), 1, f);
      fwrite (&iter, sizeof (iter), 1, f);
      fwrite (u0, sizeof (real), NX * NY * NZ, f);
      fwrite (&bytes, sizeof (bytes), 1, f);
      fclose (f);
    }
  }

  free (u0);
  free (u1);

  free_map3d (u0_3);
  free_map3d (u1_3);
}

int
main (int argc, char *argv[])
{
  /* process command line arguments */
  if (argc < 2)
  {
    printf ("Please say something like: sushkov.exe something\n");
    exit (1);
  }
#if defined(__MIC__)
  printf("COMPILED FOR MIC\n");
#endif
  printf ("Will use %i openmp threads\n", omp_get_max_threads ());
  damn_compute (argv[1]);	/* if argv[1] names a real file, it will be read as a dt file */
}


void
initParams (struct Params *p)
{
  //sizes and layer thicknesses in nm

  p->t_ins = 25;//30;		//60.;//25.;//30;//50.;          // space between top gate and split gate
  p->t_split_gate = 40;//30;		//50.;//20.;//10.;      //thickness of split gate
  p->t_GaAs = 10;		//50;//50;//50;//50;//100;// 50;//30.;         //distance between split gate and 2DEG = AlGaAs layer thickness
  p->t1_AlGaAs = 27;		//50;//50;//50;//50;//100;// 50;//30.;         //distance between split gate and 2DEG = AlGaAs layer thickness
  p->t_w = 16;		//50;//50;//50;//50;//100;// 50;//30.;         //distance between split gate and 2DEG = AlGaAs layer thickness
  p->t2_AlGaAs = 100;		//50;//50;//50;//50;//100;// 50;//30.;         //distance between split gate and 2DEG = AlGaAs layer thickness
  p->t2_GaAs = 100.;		//100.;//70;//30.;            // GaAs layer thickness
  p->Vsg = -0.5;//-0.55;//0.45;//-0.2;//-0.6;//-0.75;//0.1;//-0.7;//-0.7;//-0.4;//-0.2;//-0.362;//-0.661;//-0.525;//-0.853+0.0525;//-0.72-0.2+0.0675;//-0.46;//-0.50-3*0.0579;//0.181;//3*0.06;//-0.6;//-0.9;//-1.3;//-1.1;//-1.;//-1.1;//-1.1;//-1.0;//-1.2;			//+5;
  p->Vtg = 1.06;//-1.418;//-1.4+3*0.05;//-1.22;//+1;//-1.1;//+0.1;//2.5;//1.4;//1.3;//1.20;//1;4.;			//5.;//-5.;
//0-- -0.55;
//0a-- -.525 definition
  //grid parameters (nm)
  p->hz = 1.;			// step along z
  p->hzb = 1.;			// step along z
  p->hzw = 1;//0.5;			// step along z in well
  p->hz_gate  = 1.;		//2.;  // step along around gate
  p->hz_ins = 1.;		//2.;  // step along z in insulator
//    xmax=2038.85309382;//2000.;
//    ymax=2001.12505034;//2000.;

  p->xmax= 1500;//2038.85309382/2;//1200.;
  p->ymax= 1200;//2001.12505034/2;//1200.;
  p->xmin= -1500;//-2038.85309382/2;//-1200.;
  p->ymin= -1200;//-2001.12505034/2;//-1200.;
//  p->xmax= 1500;//2038.85309382/2;//1200.;
//  p->ymax= 1200;//2001.12505034/2;//1200.;
//  p->xmin= -1500;//-2038.85309382/2;//-1200.;
//  p->ymin= -1200;//-2001.12505034/2;//-1200.;
  p->ixmin = 0;
  p->iymin = 0;
  p->ixmax = 1200;//int((p->xmax-p->xmin)/p->hx);
  p->iymax = 960;//int((p->ymax-p->ymin)/p->hy);
//  p->ixmax = 600;//int((p->xmax-p->xmin)/p->hx);
//  p->iymax = 480;//int((p->ymax-p->ymin)/p->hy);
  p->hx =  2.5;//3000./p->ixmax;//5.0844;//4.;//2.;			//step along x
  p->hy =  2.5;//2400./p->iymax;//4.99033;
//  p->hx =  3000./p->ixmax;//5.0844;//4.;//2.;			//step along x
//  p->hy =  2400./p->iymax;//4.99033;
//  p->hx =  2038.85309382/(p->ixmax+1);//5.0844;//4.;//2.;			//step along x
//  p->hy =  2001.12505034/(p->iymax+1);//4.99033;
  p->iz_ins1 = (int) (p->t_ins/ p->hz_ins);
  p->iz_ins2 = p->iz_ins1 + (int) (p->t_split_gate / p->hz_gate);	//46;
  p->izc = p->iz_ins2 + (int) (p->t_GaAs / p->hz_gate);	//50;
  p->iz1 = p->izc + (int) (p->t1_AlGaAs / p->hz);	//60;
  p->iz2 = p->iz1 + (int) (p->t_w / p->hzw);	//60;
  p->izb = p->iz2 + (int) (p->t2_AlGaAs / p->hzb);	//70;

  p->eps_gaalas = 12.1;		//12.2;        /* GaAlAs dielectric constant */
  p->eps_gaas = 13;//13.2;		// GaAs dielectric constant
  p->eps_ins = 8;//4.35;		//8.;//3.;          // insulator dielectric constant

  p->NX = p->ixmax - p->ixmin + 1;
  p->NY = p->iymax - p->iymin + 1;
  p->NZ = p->izb + 1;

};
