/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (2.00) January 5, 2011, 11:55

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
#include "SatJ2.h"

double energy(double *v, double *p)
{
	double r, cener, pener, j2ener, ener;
	r = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	cener =  (v[3]*v[3]+v[4]*v[4]+v[5]*v[5])/2.;
	pener = -p[0]/r;
	j2ener = v[2]/r; 
	j2ener = (3*j2ener*j2ener -1)/2.;
	j2ener = p[0]*p[1]*p[1]*p[2]*j2ener/(r*r);
	ener = cener+pener+j2ener;
	return ener;
}

int main() {

/* --- PARAMETERS  --------------- */
	int npar = 3;
	double p[npar];
	p[0] = 0.005530428042714393 ; 
	p[1] = 1. ; 
	p[2] = 0.0010826266835531513 ; 

/* --- VARIABLES   --------------- */
	int nvar = 6;
	double v[nvar];
	v[0] = 1.3 ; 
	v[1] = 0.0 ; 
	v[2] = 0.0 ; 
	v[3] = 0.0 ; 
	v[4] = 0.06423314045257492 ; 
	v[5] = 0.011326035717425298 ; 

/* --- NUMBER OF FUNCTIONS   ----- */
	int nfun = 2;

/* --- TOLERANCES  --------------- */
	double tolrel = 1.e-16 ;
	double tolabs = 1.e-16 ;

/* --- INTEGRATION POINTS   ------ */
	double tini = 0.0;
	double dt   = 25.;
	double tend = 125.;
	int    nipt = (int) floor ((tend-tini)/dt);

/* --- OUTPUT  ------------------- */
	dp_data_matrix datj2;

/* --- INTEGRATOR  --------------- */
	dp_tides_delta(SatJ2, NULL, nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, &datj2, NULL);
	
/* --- Written manually  --------- */
	int i,j;
	double var[6], ener;
	for(i = 0 ; i <= nipt; i++) {
		for(j = 0; j <6; j++)  var[j] = datj2.data[i][j+1];
		ener = energy(var,p);
		printf("T+V = %.15le, H1 = %.15le, H1 - H2 = %.15le\n", datj2.data[i][7], datj2.data[i][8], datj2.data[i][8]-ener);
	}

/* --- END  ---------------------- */
	return 0;
}



