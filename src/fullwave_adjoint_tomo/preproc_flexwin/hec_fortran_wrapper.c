#include "hec.h"

// fortran wraper 

// int convolve_(const float *data, int *ndata,
// 	      const float *syn,  int *nsyn,
// 	      float *conv)
// {
//   return convolve(data,*ndata,syn,*nsyn,conv);
// }
// int CONVOLVE_(const float *data, int *ndata,
// 	      const float *syn,  int *nsyn,
// 	      float *conv)
// {
//   return convolve(data,*ndata,syn,*nsyn,conv);
// }


// //
// int xcorr_(const float *data, int *pndata,
//  	   const float *syn, int *pnsyn,
// 	   float *corr, int *nshift)
// {
  
//   return xcorr(data,*pndata,syn,*pnsyn,corr,nshift);
  
// }
// int XCORR_(const float *data, int *pndata,
//  	   const float *syn, int *pnsyn,
// 	   float *corr, int *nshift)
// {
  
//   return xcorr(data,*pndata,syn,*pnsyn,corr,nshift);
  
// }

// //
// // return the index of zero-shift in fortran
// int crosscorr_(const float *data, int *pndata,
// 	   const float *syn, int *pnsyn,
// 	   float *corr)
// {
  
//   return crosscorr(data,*pndata,syn,*pnsyn,corr)+1;

// }
// int CROSSCORR_(const float *data, int *pndata,
// 	   const float *syn, int *pnsyn,
// 	   float *corr)
// {
  
//   return crosscorr(data,*pndata,syn,*pnsyn,corr)+1;

// }

// //
// int autocorr_(const float *data, int *pndata,
// 	   float *corr)
// {
 
//   return autocorr(data,*pndata,corr)+1;

// }
// int AUTOCORR_(const float *data, int *pndata,
// 	   float *corr)
// {
 
//   return autocorr(data,*pndata,corr)+1;

// }

// //
// int autocorr2_(const float *data, int *pndata,
// 	   float *corr)
// {
 
//   return autocorr2(data,*pndata,corr)+1;

// }
// int AUTOCORR2_(const float *data, int *pndata,
// 	   float *corr)
// {
 
//   return  autocorr2(data,*pndata,corr)+1;

// }

//
int envelope_(const int *numsamp,const float *f,float *fenv)

{

  return envelope(*numsamp,f,fenv);

}
int ENVELOPE_(const int *numsamp,const float *f,float *fenv)

{

  return envelope(*numsamp,f,fenv);

}

//
// int hilbert_(const int *numsamp,float *f,float *fhilb)

// {

//   return hilbert(*numsamp,f,fhilb);

// }
// int HILBERT_(const int *numsamp,float *f,float *fhilb)

// {

//   return hilbert(*numsamp,f,fhilb);

// }
