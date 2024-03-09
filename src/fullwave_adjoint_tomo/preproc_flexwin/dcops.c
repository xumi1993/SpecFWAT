
#include <math.h>
#include "dcomplex.h"

/*---------------------------------------------------------------------*/
dcomplex dcmult(a,b)		            /* complex multiplication  */
	dcomplex a,b;
{
	dcomplex c;

		c.re = (a.re * b.re) - (a.im * b.im);
			c.im = (a.re * b.im) + (a.im * b.re);
				return(c);

} 

/*---------------------------------------------------------------------*/
dcomplex dconj(a)		               /* complex conjugation  */
	dcomplex a;
{
	dcomplex abar;

		abar.re = a.re;
			abar.im = -a.im;
				return(abar);

}

/*---------------------------------------------------------------------*/
dcomplex dcdiv(a,b)		                  /* complex division  */
	dcomplex a,b;

{
	dcomplex c,temp;
	dcomplex dconj(),dcmult(),dscadiv();
	   
		temp = dconj(b);
			c = dcmult(a,temp);
				temp.re = b.re*b.re + b.im*b.im;
					temp.im = 0.e0;
						c = dscadiv(c,temp.re);
							return(c);

} 

/*---------------------------------------------------------------------*/
dcomplex dcsum(a,b)	                	  /* complex addition  */
	dcomplex a,b;
{
	dcomplex c;

		c.re = a.re + b.re;
			c.im = a.im + b.im;
				return(c);

}  
 
/*---------------------------------------------------------------------*/
dcomplex dcdiff(a,b)		               /* complex subtraction  */
	dcomplex a,b;
{
	dcomplex c;

		c.re = a.re - b.re;
			c.im = a.im - b.im;
				return(c);

}  

/*---------------------------------------------------------------------*/
dcomplex dscamult(z,c)	 /* scalar multiplication of a complex number  */
	dcomplex z;
	double c;
{
	dcomplex ans;

		ans.re = c * z.re;
			ans.im = c * z.im;
				return(ans);

}

/*---------------------------------------------------------------------*/
dcomplex dscadiv(z,c)          /* scalar division of a complex number  */
	dcomplex z;
	double c;
{
	dcomplex ans;

		ans.re =  z.re / c ;
			ans.im =  z.im / c ;
				return(ans);

} 


/*---------------------------------------------------------------------*/
double dmodu(z)	                       /* modulus of a complex number  */
	dcomplex z;
{
	double ans;

		ans = (z.re * z.re) + (z.im * z.im);
			ans = sqrt((double) ans);
				return(ans);

}  

/*---------------------------------------------------------------------*/
double dcphase(c)                        /* phase of a complex number  */
	dcomplex *c;
{
	double phase;
	float pi = 3.1415927;

	     if(c->re == 0.)
		            if(c->im == 0.)
				    	 phase = 0.;
	            else if(c->im < 0.)
			    	 phase = -pi/2.;
	            else
			    	 phase = pi/2.;
	          else 
			         if(c->im == 0. && c->re < 0.)
					 	 phase = pi;
		       else
			              {
					             phase = (c->im) / (c->re);
						            phase = atan(phase);
							           }
		       if(c->re < 0. && c->im <0.)
			              phase = phase + pi;
		          
		            return(phase);
} 

/*---------------------------------------------------------------------*/
dcomplex dsca_add(k0,k)           /* sum of scalar and complex number  */
	double k0;
	dcomplex k;
{
	dcomplex c;

		c.re = k.re + k0;
			c.im = k.im;
				return(c);
} 

/*---------------------------------------------------------------------*/
dcomplex dcexpon(c)                /* exponential of a complex number  */
	dcomplex c;
{
	dcomplex expc;
	expc.re = exp(c.re)*cos(c.im);
	expc.im = exp(c.re)*sin(c.im);
	return (expc);
}

/*---------------------------------------------------------------------*/
dcomplex dccos(c)                          /* cos of a complex number  */
	dcomplex c;
{
	double fac1, fac2, fac3, fac4;
	dcomplex dccos;
	fac1 = cos(c.re);
	fac2 = sin(c.re);
	fac3 = exp(c.im);
	fac4 = exp(-c.im);
	dccos.re = (fac1*fac3 + fac1*fac4)/2.;
	dccos.im = (fac2*fac3 - fac2*fac4)/2.;
	return(dccos);
}

/*---------------------------------------------------------------------*/
dcomplex dcsin(c)                          /* sin of a complex number  */
	dcomplex c;
{
	double fac1, fac2, fac3, fac4;
	dcomplex dcsin;
	fac1 = sin(c.re);
	fac2 = cos(c.re);
	fac3 = exp(c.im);
	fac4 = exp(-c.im);
	dcsin.re = (fac1*fac3 + fac1*fac4)/2.;
	dcsin.im = (fac2*fac3 - fac2*fac4)/2.;
	return(dcsin);
}

