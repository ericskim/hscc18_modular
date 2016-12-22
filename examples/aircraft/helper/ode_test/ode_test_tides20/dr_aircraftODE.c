/****************************************************************************
	Driver file of the dp_tides program
	This file has been created by MathTIDES (2.00) January 20, 2016, 9:21

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
#include "aircraftODE.h"

int main() {

/* --- PARAMETERS  --------------- */
	int npar = 2;
	double p[npar];
	p[0] = ***** ; 
	p[1] = ***** ; 

/* --- VARIABLES   --------------- */
	int nvar = 3;
	double v[nvar];
	v[0] = 1. ; 
	v[1] = 1. ; 
	v[2] = 1. ; 

/* --- NUMBER OF FUNCTIONS   ----- */
	int nfun = 0;

/* --- TOLERANCES  --------------- */
	double tolrel = 1.e-16 ;
	double tolabs = 1.e-16 ;

/* --- INTEGRATION POINTS   ------ */
	double tini = 0.0;
	double tend = 0.25;
	int    nipt = 1;
	double dt = (tend - tini)/nipt ;

/* --- OUTPUT  ------------------- */
	FILE* fd = stdout;

/* --- INTEGRATOR  --------------- */
	dp_tides_delta(aircraftODE, NULL, nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, NULL, fd);

/* --- END  ---------------------- */
	return 0;
}


