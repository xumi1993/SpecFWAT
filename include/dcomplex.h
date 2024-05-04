/* double-precision complex arithmetic subroutines */
#ifndef _dcomplex_h
#define _dcomplex_h

typedef struct { double re; double im;}  dcomplex;


dcomplex dcmult(dcomplex a,dcomplex b);
dcomplex dconj(dcomplex a);
dcomplex dcdiv(dcomplex a,dcomplex b);
dcomplex dcsum(dcomplex a,dcomplex b);
dcomplex dcdiff(dcomplex a,dcomplex b);
dcomplex dscamult(dcomplex z,double c);
dcomplex dscadiv(dcomplex z,double c);
double dmodu(dcomplex z);
double dcphase(dcomplex *c);

dcomplex dsca_add(double k0, dcomplex k);
dcomplex dcexpon(dcomplex c);
dcomplex dccos(dcomplex c);
dcomplex dcsin(dcomplex c);

// fft prototype
void cfft(dcomplex *x,int n,int isign);


#endif
