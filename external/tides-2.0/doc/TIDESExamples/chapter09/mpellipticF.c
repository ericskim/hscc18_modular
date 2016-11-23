#include "mpfr.h"
#include "mp_tides.h"
#include "ellipMP.h"

void ellipticF(mpfr_t ellipf,  mpfr_t phi, mpfr_t k);


int main()
{
	int dig = 30;
	set_precision_digits(dig);
	int i, npoints = 5;
	mpfr_t phi, par, dphi, ppar, ellipF;
	mpfr_init2(phi,   TIDES_PREC);
	mpfr_init2(dphi,  TIDES_PREC);
	mpfr_init2(par,   TIDES_PREC);
	mpfr_init2(ppar,  TIDES_PREC);
	mpfr_init2(ellipF,TIDES_PREC);
	
	mpfr_set_str(dphi,"1.57079632679489661923132169164", 10, TIDES_RND);
	mpfr_div_si(dphi, dphi, npoints, TIDES_RND);

	printf("\nElliptic integral of first kind: \n");
	mpfr_set_str(ppar,"0.1", 10, TIDES_RND);

	for(i = 0; i <= npoints; i++) {
		mpfr_mul_si(phi, dphi, i, TIDES_RND);
		mpfr_mul_si(par, ppar, i, TIDES_RND);
		ellipticF(ellipF, phi, par);
		mpfr_printf("F(%d Pi/2, %.2Rf) = %.29Re\n", i, par, ellipF);
	}
	return 0; 
}

void ellipticF(mpfr_t ellipf,  mpfr_t phi, mpfr_t k)
{

	int pphi,pk, prec, dprec;
	pphi = (int)mpfr_get_prec (phi);
	pk = (int)mpfr_get_prec (phi);
	if(pphi > pk) prec = pk;
	else prec = pphi;

	mpfr_t  x, par, tini, tol;
	mpfr_init2(x,	  prec);
	mpfr_init2(par,	  prec);
	mpfr_init2(tini,  prec);
	mpfr_init2(tol,   prec);
	
	mpfr_set_str(x,   "0", 10, TIDES_RND);
	mpfr_set_str(tini,"0", 10, TIDES_RND);
	mpfr_set(par, k, TIDES_RND);
	
	dprec = floor(prec/3.3219);
	mpfr_set_si(tol, -dprec, TIDES_RND);
	mpfr_exp10(tol,  tol, TIDES_RND);
	
	mp_tides_delta(ellipMP, NULL, 1, 1, 0, &x, &par, tini, phi, 1, tol, tol, NULL, NULL);
	mpfr_set(ellipf, x, TIDES_RND);

}

