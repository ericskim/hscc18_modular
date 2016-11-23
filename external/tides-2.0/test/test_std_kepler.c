/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (2.00) January 26, 2011, 18:16

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
#include "std_kepler.h"

int main() {

/* --- PARAMETERS  --------------- */
	int npar = 1;
	double p[npar];
	p[0] = 1. ; 

/* --- VARIABLES   --------------- */
	int nvar = 4;
	double v[nvar];
	v[0] = 0.30000000000000004 ; 
	v[1] = 0.0 ; 
	v[2] = 0.0 ; 
	v[3] = 2.3804761428476167 ; 

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


/* --- INTEGRATOR  --------------- */
	double x[nvar];
	int i;
	for (i=0; i<nvar; i++) x[i] = v[i];
	
//	set_info_taylor();


	dp_tides_delta(std_kepler, NULL, nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, NULL, NULL);

	double aux = 0.;
	double error = -1.;
	
	for (i=0; i<nvar; i++) {
		aux = fabs(v[i]-x[i]);
		if (aux > error) error = aux;
	}
	
	
	if (error < 1.e-10) {
		return 0;
	} else {
		printf("Error in test_std_kepler: %25.17e\n",error);
		return 1;
	}
}


