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
#include "std_lorenz_partials.h"

int   std_lorenz_partials_PDData[]  = {2, 4, 6, 6, 7, 0, 1, 3, 5, 8, 12, 15, 15, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 15, 0, 0, 1, 0, 2, 0, 1, 3, 0, 2, 1, 4, 0, 2, 5, 15, 0, 1, 0, 2, 0, 3, 1, 0, 4, 1, 2, 0, 5, 2, 0, 7, 0, 1, 2, 3, 5, 7, 9, 9, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 0, 0, 0, 0, 1, 0, 2, 0, 2, 9, 0, 1, 2, 3, 1, 4, 1, 5, 2, 0, 0, 1, 0, 0, 1, 2, 0, 1, 1, 0, 2};


int main() {

/* --- PARAMETERS  --------------- */
	int npar = 3;
	double p[npar];
	p[0] = 10. ; 
	p[1] = 28. ; 
	p[2] = 2.666666666666667 ; 

/* --- VARIABLES   --------------- */
	int nvar = 3;
	double v[nvar];
	v[0] = -13.7636106821342 ; 
	v[1] = -19.5787519424518 ; 
	v[2] = 27. ; 

/* --- NUMBER OF FUNCTIONS   ----- */
	int nfun = 0;

/* --- TOLERANCES  --------------- */
	double tolrel = 1.e-16 ;
	double tolabs = 1.e-16 ;

/* --- INTEGRATION POINTS   ------ */
	double tini = 0.0;
	double tend = 1.558652210716175;
	int    nipt = 1;
	double dt = (tend - tini)/nipt ;

/* --- OUTPUT  ------------------- */
	dp_data_matrix std_lorenz_partials_DataMatrix;

/* --- INTEGRATOR  --------------- */
	double sol[19];
	sol[0] =  1.55865221071617510e+00;
	sol[1] = -1.37636106821342086e+01;
	sol[2] = -1.95787519424518131e+01;
	sol[3] = 2.69999999999999964e+01;
	sol[4] = -9.12802195703595109e-01;
	sol[5] = -1.30905810340089968e-01;
	sol[6] = 1.67609907947194769e+00;
	sol[7] = -2.82268059307864192e+01;
	sol[8] = 1.81841834628623644e+01;
	sol[9] = 1.16876988007583165e+02;
	sol[10] = 1.96143936526489815e-01;
	sol[11] = 3.61143595741667101e-01;
	sol[12] = -1.19613270510495917e-01;
	sol[13] = -2.45466501786646152e+00;
	sol[14] = 1.40187669688222609e+01;
	sol[15] = 2.25897033575195572e+01;
	sol[16] = 6.08714526721828975e+02;
	sol[17] = 8.47488020397782179e+02;
	sol[18] = -1.25071254195740835e+03;
	
	dp_tides_delta(std_lorenz_partials, std_lorenz_partials_PDData, 
				   nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, 
				   &std_lorenz_partials_DataMatrix, NULL);

	double aux = 0.;
	double error = -1.;
	int i;

	for (i=0; i<19; i++) {
		aux = fabs(sol[i]-std_lorenz_partials_DataMatrix.data[1][i]);
		if (aux > error) error = aux;
	}
	
	
	if (error < 1.e-10) {
		return 0;
	} else {
		printf("Error in test_std_lorenz_partials: %25.17e\n",error);
		return 1;
	}
}


