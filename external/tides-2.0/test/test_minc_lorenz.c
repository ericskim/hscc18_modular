/****************************************************************************
	Driver file of the minc_tides program
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

#include "minc_tides.h"

int main() {

	int  i, VARS, PARS; 
	VARS = 3;
	PARS = 3;
	double tolrel, tolabs, tini, tend, dt; 
	double v[VARS], p[PARS]; 


/************************************************************/
/************************************************************/
/*      INITIAL CONDITIONS, INTEGRATION TIMES, TOLERANCES    */
/************************************************************/
/************************************************************/

/* --- PARAMETERS VALUE --- */
	p[0] = 10. ; 
	p[1] = 28. ; 
	p[2] = 2.666666666666667 ; 

/* --- INITIAL VALUES --- */
	v[0] = -13.7636106821342 ; 
	v[1] = -19.5787519424518 ; 
	v[2] = 27. ; 

/* --- INITIAL INTEGRATION POINT --- */
	tini = 0.0 ;

/* --- ENDPOINT OF INTEGRATION   --- */
	tend = 1.558652210716175 ;

/* --- DELTA t FOR DENSE OUTPUT  --- */
	dt   = 1.558652210716175 ;

/* --- REQUIRED TOLERANCES --- */
	tolrel = 1.e-16 ;
	tolabs = 1.e-16 ;


/***********************************************************/
/***********************************************************/
/*       CALL THE INTEGRATOR                               */
/***********************************************************/
/***********************************************************/
	double x[VARS];
	for (i=0; i<VARS; i++) x[i] = v[i];
	minc_tides(v,VARS,p,PARS,tini,tend,dt,tolrel,tolabs);

/***********************************************************/
/***********************************************************/
/*       SHOW FINAL POINT ON THE SCREEN                    */
/*       SHOW ACCEPTED AND REJECTED STEPS ON THE SCREEN    */
/***********************************************************/
/***********************************************************/

	double aux = 0.;
	double error = -1.;
	
	for (i=0; i<VARS; i++) {
		aux = fabs(v[i]-x[i]);
		if (aux > error) error = aux;
	}
	
	
	if (error < 1.e-10) {
		return 0;
	} else {
		printf("Error in test_minc_lorenz: %25.17e\n",error);
		return 1;
	}
	
	return 0;
}




