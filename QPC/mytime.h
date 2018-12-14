/* -*- c++ -*- */
#ifndef MYTIME_H_
#define MYTIME_H_

#if !defined(INLINE)
#if defined(__cplusplus)
#define INLINE inline
#elif defined(__GNUC__)
#define INLINE __inline__
#elif defined(__INTEL_COMPILER)
#define INLINE __inline
#else
#define INLINE static
#endif
#endif

//-------------------------------------------------------------------
#if defined(_WIN32)
#define WINDOWS_LEAN_AND_MEAN
#include <windows.h>
typedef long long int mytime_t;
INLINE mytime_t mytime()
{
  LARGE_INTEGER t;
  QueryPerformanceCounter( &t );
  return (mytime_t)(t.QuadPart);
}
static double mytime_invhz()
{
  static double invhz;
  if (!invhz)
  {
    LARGE_INTEGER t;
    QueryPerformanceFrequency(&t);
    invhz = 1.0 / t.QuadPart;
  }
  return invhz;
}
INLINE double mytime_s(mytime_t t0,mytime_t t1)
{
  return mytime_invhz() * (t1-t0);
}
INLINE double mytime_ms(mytime_t t0,mytime_t t1)
{
  return 1e3 * mytime_invhz() * (t1-t0);
}
INLINE double mytime_us(mytime_t t0,mytime_t t1)
{
  return 1e6 * mytime_invhz() * (t1-t0);
}
INLINE double mytime_ns(mytime_t t0,mytime_t t1)
{
  return 1e9 * mytime_invhz() * (t1-t0);
}
#endif /* _WIN32 */
//-------------------------------------------------------------------
#if defined(__linux__)
#include <time.h>
typedef long long int mytime_t;
INLINE mytime_t mytime()
{
    struct timespec tp;
    clock_gettime(CLOCK_MONOTONIC, &tp);
    return (mytime_t)(1e9*tp.tv_sec + tp.tv_nsec);
}
INLINE double mytime_s(mytime_t t0,mytime_t t1)
{
    return 1e-9*(t1-t0);
}
INLINE double mytime_ms(mytime_t t0,mytime_t t1)
{
    return 1e-6*(t1-t0);
}
INLINE double mytime_us(mytime_t t0,mytime_t t1)
{
    return 1e-3*(t1-t0);
}
INLINE double mytime_ns(mytime_t t0,mytime_t t1)
{
    return (t1-t0);
}
#endif /* __linux__ */
//-------------------------------------------------------------------

#endif /* MYTIME_H_ */
