/****************************************************************************
	Driver file of the mp_tides program
	This file has been created by MathTIDES (2.00) January 27, 2011, 18:43

	Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
	Grupo de Mecanica Espacial
	University of Zaragoza
	SPAIN

	http://gme.unizar.es/software/tides
	Contact: <tides@unizar.es>

	This file is part of TIDES.

	TIDES is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	TIDES is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with TIDES.  If not, see <http://www.gnu.org/licenses/>.

*****************************************************************************/

#include "mpfr.h"
#include "mp_tides.h"
#include "mpfr_lorenz.h"

int main() {

	int i;

/* --- SET PRECISION  ------------ */
	set_precision_digits(50);

/* --- PARAMETERS  --------------- */
	int npar = 3;
	mpfr_t p[npar];
	for(i=0; i<npar; i++) mpfr_init2(p[i], TIDES_PREC);
	mpfr_set_str(p[0], "10.", 10, TIDES_RND); 
	mpfr_set_str(p[1], "28.", 10, TIDES_RND); 
	mpfr_set_str(p[2], "2.6666666666666666666666666666666666666666666666667", 10, TIDES_RND); 

/* --- VARIABLES   --------------- */
	int nvar = 3;
	mpfr_t v[nvar];
	for(i=0; i<nvar; i++) mpfr_init2(v[i], TIDES_PREC);
	mpfr_set_str(v[0], "-13.763610682134200525014401054361653864100864854092", 10, TIDES_RND); 
	mpfr_set_str(v[1], "-19.578751942451795538838041446009558866114240053428", 10, TIDES_RND); 
	mpfr_set_str(v[2], "27.", 10, TIDES_RND); 

/* --- NUMBER OF FUNCTIONS   ----- */
	int nfun = 0;

/* --- TOLERANCES  --------------- */
	mpfr_t tolrel, tolabs;
	mpfr_init2(tolrel, TIDES_PREC); 
	mpfr_init2(tolabs, TIDES_PREC); 
	mpfr_set_str(tolrel, "1.e-49", 10, TIDES_RND);
	mpfr_set_str(tolabs, "1.e-49", 10, TIDES_RND);

/* --- INTEGRATION POINTS   ------ */
	mpfr_t tini, dt, tend; 
	mpfr_init2(tini, TIDES_PREC); 
	mpfr_init2(dt, TIDES_PREC); 
	mpfr_init2(tend, TIDES_PREC); 
	mpfr_set_str(tini, "0", 10, TIDES_RND);
	mpfr_set_str(tend, "1.5586522107161747275678702092126960705284805489972", 10, TIDES_RND);
	int  nipt  = 1;
	mpfr_set(dt,tend, TIDES_RND); 
	mpfr_sub(dt,dt,tini, TIDES_RND); 
	mpfr_div_si(dt,dt,nipt, TIDES_RND); 

/* --- INTEGRATOR  --------------- */
	mpfr_t x[nvar];
	for (i=0; i<nvar; i++) {
		mpfrts_init(&x[i]);
		mpfrts_set(&x[i], v[i]);
	}
	
//	mpfrts_set_info_taylor();
	mp_tides_delta(mpfr_lorenz, NULL, nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, NULL, NULL);

	mpfr_t aux, error, errm;
	mpfrts_init(&aux);
	mpfrts_init(&error);
	mpfrts_init(&errm);
	mpfr_set_str(errm, "1.e-45", 10, TIDES_RND);
	
	mpfrts_set_str (&aux, "0."); mpfrts_set_str (&error, "-1.");
	for (i=0; i<nvar; i++) {
		mpfrts_sub (&aux, v[i],x[i]);
		mpfrts_abs (&aux, aux);
		/* mpfrts_div (&aux, aux, x[i]);*/
		if (mpfrts_greater (aux, error))
			mpfrts_set (&error, aux);
	}
	
	if (mpfrts_greater (errm, error)) {
		return 0;
	} else {
		mpfrts_write ("Error in test_mpfr_lorenz", error);
		return 1;
	}
}


