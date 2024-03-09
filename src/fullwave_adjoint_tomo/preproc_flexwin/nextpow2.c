#include <limits.h>
#include <stdio.h>

int  nextpow2(int numsamp)
{
  int npts_max;
  
  if (numsamp > INT_MAX) {
    fprintf(stderr,"Too long time series\n");
    return 0;
  }
  npts_max = 1;
  while( 1 )
    {
      npts_max *= 2;
      if (npts_max >= numsamp) return npts_max;
    }
}

