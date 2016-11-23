/****************************************************************************
	This file has been created by MathTIDES (2.00) October 28, 2015, 11:44

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
#include "vehicleODE.h"


long  vehicleODE(iteration_data *itd, double t, double v[], double p[], int ORDER, double *cvfd)
{

	int i;
	static int   VARIABLES        = 3;
	static int   PARAMETERS       = 2;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 19;
	static int   POS_FUNCTIONS[1] = {0};

	initialize_dp_case();

	double ct[] = {0.5, 1., 4.};


	for(i=0;  i<=ORDER; i++) {
		double_var_t(itd, link[18],var[1], i);
		double_var_t(itd, link[17],var[2], i);
		double_var_t(itd, link[6],var[3], i);
		double_sin_cos_t(itd, par[1],link[1],link[0],i);
		double_mul_t_cc(itd, ct[0],par[0],link[2],i);
		double_div_t(itd, link[1],link[0],link[3],i);
		double_mul_t_cc(itd, ct[0],link[3],link[4],i);
		double_mul_t(itd, link[3],link[3],link[5],i);
		double_mul_t(itd, link[3],par[0],link[6],i);
		double_mul_t(itd, link[4],link[4],link[7],i);
		double_add_t_cc(itd, ct[1],link[7],link[8],i);
		double_atan_t(itd, link[4],link[8],link[9],i);
		double_add_t_cc(itd, ct[2],link[5],link[10],i);
		double_add_t(itd, link[9],var[3],link[11],i);
		double_pow_t_cc(itd, link[10],ct[0],link[12],i);
		double_sin_cos_t(itd, link[11],link[14],link[13],i);
		double_mul_t(itd, link[2],link[12],link[15],i);
		double_mul_t(itd, link[2],link[13],link[16],i);
		double_mul_t(itd, link[14],link[15],link[17],i);
		double_mul_t(itd, link[12],link[16],link[18],i);
	}

	write_dp_solution();

	return NUM_COLUMNS;
}

