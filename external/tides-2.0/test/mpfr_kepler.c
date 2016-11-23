/****************************************************************************
	This file has been created by MathTIDES (2.00) April 16, 2011, 20:06

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

#include "mp_tides.h"
#include "mpfr_kepler.h"


long  mpfr_kepler(iteration_data *itd, mpfr_t t, mpfr_t v[], mpfr_t p[], int ORDER, mpfr_t *cvfd)
{

	int i;
	static int   VARIABLES        = 4;
	static int   PARAMETERS       = 1;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 10;
	static int   POS_FUNCTIONS[1] = {0};

	initialize_mp_case();

	int NCONST = 2;
	mpfr_t ct[NCONST];
	for(i = 0; i < NCONST ; i++ ) { 
		mpfrts_init(&ct[i]);
	}
	mpfrts_set_i(&ct[0],-3);
	mpfrts_div_i(&ct[0],ct[0],2);
	mpfrts_set_i(&ct[1],-1);

	for(i=0;  i<=ORDER; i++) {
		mpfrts_var_t(itd, var[3],var[1], i);
		mpfrts_var_t(itd, var[4],var[2], i);
		mpfrts_var_t(itd, link[8],var[3], i);
		mpfrts_var_t(itd, link[9],var[4], i);
		mpfrts_mul_t_cc(itd, ct[1],var[1],link[0],i);
		mpfrts_mul_t_cc(itd, ct[1],var[2],link[1],i);
		mpfrts_mul_t(itd, var[1],var[1],link[2],i);
		mpfrts_mul_t(itd, var[2],var[2],link[3],i);
		mpfrts_add_t(itd, link[2],link[3],link[4],i);
		mpfrts_mul_t(itd, link[0],par[0],link[5],i);
		mpfrts_mul_t(itd, link[1],par[0],link[6],i);
		mpfrts_pow_t_cc(itd, link[4],ct[0],link[7],i);
		mpfrts_mul_t(itd, link[5],link[7],link[8],i);
		mpfrts_mul_t(itd, link[6],link[7],link[9],i);
	}

	write_mp_solution();
	clear_vpl();
	clear_cts();

	return NUM_COLUMNS;
}

