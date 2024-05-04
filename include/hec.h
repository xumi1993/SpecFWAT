#ifndef _hec_h
#define _hec_h

#ifndef _dcomplex_h
#define _dcomplex_h
#endif

// spectrum analysis
int  nextpow2(int numsamp);

int envelope(int numsamp,const float *f,float *fenv);

int hilbert(int numsamp,const float *f,float *fhilb);

void detrend(float *x, int n);

//convolution
int convolve(const float *data, int ndata,
	     const float *syn,  int nsyn,
	     float *conv);

//cross-correlation
int xcorr(const float *data, int ndata,
          const float *syn, int nsyn,
          float *corr, int *nshift);

int crosscorr(const float *data, int ndata,
              const float *syn,  int nsyn,
              float *corr);
int autocorr(const float *data, int ndata,
             float *corr);
int autocorr2(const float *data, int ndata,
	      float *corr);

//fortran wrapper

int convolve_c(const float *data, int *ndata,
              const float *syn,  int *nsyn,
	       float *conv);
int CONVOLVE_C(const float *data, int *ndata,
              const float *syn,  int *nsyn,
	       float *conv);
int xcorr_c(const float *data, int *pndata,
           const float *syn, int *pnsyn,
	    float *corr, int *nshift);
int XCORR_C(const float *data, int *pndata,
           const float *syn, int *pnsyn,
	    float *corr, int *nshift);
int crosscorr_c(const float *data, int *pndata,
           const float *syn, int *pnsyn,
		float *corr);
int CROSSCORR_C(const float *data, int *pndata,
           const float *syn, int *pnsyn,
		float *corr);
int autocorr_c(const float *data, int *pndata,
	       float *corr);
int AUTOCORR_C(const float *data, int *pndata,
	       float *corr);
int autocorr2_c(const float *data, int *pndata,
		float *corr);
int AUTOCORR2_C(const float *data, int *pndata,
		float *corr);

int envelope_c(const int *numsamp,const float *f,float *fenv);
int ENVELOPE_C(const int *numsamp,const float *f,float *fenv);
int hilbert_c(const int *numsamp,float *f,float *fhilb);
int HILBERT_C(const int *numsamp,float *f,float *fhilb);


#endif
