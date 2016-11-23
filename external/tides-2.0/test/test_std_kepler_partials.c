/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (2.00) January 28, 2011, 11:35

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

#include "dp_tides.h"
#include "std_kepler_partials.h"

int   std_kepler_partials_PDData[]  = {2, 1, 2, 6, 7, 0, 1, 3, 5, 8, 12, 15, 15, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 15, 0, 0, 1, 0, 2, 0, 1, 3, 0, 2, 1, 4, 0, 2, 5, 15, 0, 1, 0, 2, 0, 3, 1, 0, 4, 1, 2, 0, 5, 2, 0, 7, 0, 1, 2, 3, 5, 7, 9, 9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 0, 0, 0, 0, 1, 0, 2, 0, 2, 9, 0, 1, 2, 3, 1, 4, 1, 5, 2, 0, 0, 1, 0, 0, 1, 2, 0, 1, 1, 0, 2};


int main() {

/* --- PARAMETERS  --------------- */
	int npar = 1;
	double p[npar];
	p[0] = 1. ; 

/* --- VARIABLES   --------------- */
	int nvar = 4;
	double v[nvar];
//	v[0] = 0.30000000000000004 ; 
	v[1] = 0.0 ; 
	v[2] = 0.0 ; 
//	v[3] = 2.3804761428476167 ; 
	v[0] = .99;
	v[3] = 1.010050503787816;
/* --- NUMBER OF FUNCTIONS   ----- */
	int nfun = 0;

/* --- TOLERANCES  --------------- */
	double tolrel = 1.e-16 ;
	double tolabs = 1.e-16 ;

/* --- INTEGRATION POINTS   ------ */
	double tini = 0.0;
	double tend = 62.83185307179586;
	int    nipt = 1;
	double dt = (tend - tini)/nipt ;


/* --- OUTPUT  ------------------- */
	dp_data_matrix std_kepler_partials_DataMatrix;

/* --- INTEGRATOR  --------------- */
	dp_tides_delta(std_kepler_partials, std_kepler_partials_PDData, 
				   nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, 
				   &std_kepler_partials_DataMatrix, NULL);
	
	double sol[25];
	sol[0] =  6.2831853071795862e+01;
	sol[1] =  9.8999999999999866e-01;
	sol[2] =  5.0195958500864890e-14;
	sol[3] =  -5.0404125317982107e-14;
	sol[4] =  1.0100505037878160e+00;
	sol[5] =  1.0000000000095639e+00;
	sol[6] =  -1.9425572344378682e+02;
	sol[7] =  1.9622771409947413e+02;
	sol[8] =  9.9475983006414026e-12;
	sol[9] =  4.7169950237457847e-14;
	sol[10] =  9.9999999999767608e-01;
	sol[11] =  2.3402391136073675e-12;
	sol[12] =  -8.3183460120039854e-14;
	sol[13] =  -3.7739059996668650e+04;
	sol[14] =  -5.9856370248138046e+02;
	sol[15] =  -1.8819921564279412e+02;
	sol[16] =  -3.8503390466781886e+04;
	sol[17] =  -5.2428106300794752e-10;
	sol[18] =  1.1388969767267554e-09;
	sol[19] =  -1.1789609288825886e-09;
	sol[20] =  1.9820981222129046e+02;
	sol[21] =  4.1552539187250659e-11;
	sol[22] =  -1.9621790246971091e+02;
	sol[23] =  1.9820981222293122e+02;
	sol[24] =  -6.0005334034940461e-12;
	
	

	double aux = 0.;
	double error = -1.;
	int i;
	for (i=0; i<25; i++) {
		aux = fabs(sol[i]-std_kepler_partials_DataMatrix.data[1][i]);
		if (aux > error) error = aux;
	}
	
	
	if (error < 1.e-8) {
		return 0;
	} else {
		printf("Error in test_std_kepler_partials: %25.17e\n",error);
		return 1;
	}
}


