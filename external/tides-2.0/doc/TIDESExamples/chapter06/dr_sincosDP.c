/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (2.00) October 27, 2015, 14:11

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
#include "sincosDP.h"

int main() {

/* --- PARAMETERS  --------------- */
	int npar = 0;

/* --- VARIABLES   --------------- */
	int nvar = 2;
	double v[nvar];
	v[0] = 0.0 ; 
	v[1] = 1. ; 

/* --- NUMBER OF FUNCTIONS   ----- */
	int nfun = 0;

/* --- TOLERANCES  --------------- */
	double tolrel = 1.e-16 ;
	double tolabs = 1.e-16 ;

/* --- INTEGRATION POINTS   ------ */
	double tini = 0.0;
	double tend = 1.570796326794897;
	int    nipt = 1;
	double dt = (tend - tini)/nipt ;

/* --- OUTPUT  ------------------- */
	FILE* fd = stdout;

/* --- INTEGRATOR  --------------- */
	dp_tides_delta(sincosDP, NULL, nvar, npar, nfun, v, NULL, tini, dt, nipt, tolrel, tolabs, NULL, fd);

/* --- END  ---------------------- */
	return 0;
}


