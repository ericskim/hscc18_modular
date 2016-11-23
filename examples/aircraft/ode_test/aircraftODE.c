/****************************************************************************
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


long  aircraftODE(iteration_data *itd, double t, double v[], double p[], int ORDER, double *cvfd)
{

	int i;
	static int   VARIABLES        = 3;
	static int   PARAMETERS       = 2;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 26;
	static int   POS_FUNCTIONS[1] = {0};

	initialize_dp_case();

	double ct[] = {-9.810000000000002, -9.81, -0.00090552, -0.000539, -0.00012520833333333332, 0.000016666666666666667, 0.00001666666666666667, 0.001429166666666667, 0.004802000000000001};

	for(i=0;  i<=ORDER; i++) {
		double_var_t(itd, link[24],var[1], i);
		double_var_t(itd, link[25],var[2], i);
		double_var_t(itd, link[17],var[3], i);
		double_sin_cos_t(itd, par[1],link[2],link[0],i);
		double_sin_cos_t(itd, var[2],link[3],link[1],i);
		double_mul_t_cc(itd, ct[3],par[1],link[4],i);
		double_mul_t_cc(itd, ct[5],par[0],link[5],i);
		double_mul_t_cc(itd, ct[6],par[0],link[6],i);
		double_mul_t_cc(itd, ct[8],par[1],link[7],i);
		double_mul_t(itd, par[1],par[1],link[8],i);
		double_mul_t(itd, var[1],var[1],link[9],i);
		double_add_t_cc(itd, ct[4],link[4],link[10],i);
		double_add_t_cc(itd, ct[7],link[7],link[11],i);
		double_mul_t_cc(itd, ct[0],link[1],link[12],i);
		double_mul_t_cc(itd, ct[0],link[3],link[13],i);
		double_mul_t_cc(itd, ct[2],link[8],link[14],i);
		double_mul_t(itd, link[0],link[5],link[15],i);
		double_mul_t(itd, link[2],link[6],link[16],i);
		double_mul_t(itd, link[3],var[1],link[17],i);
		double_add_t(itd, link[10],link[14],link[18],i);
		double_add_t(itd, link[13],link[15],link[19],i);
		double_mul_t(itd, link[9],link[11],link[20],i);
		double_add_t(itd, link[12],link[20],link[21],i);
		double_mul_t(itd, link[9],link[18],link[22],i);
		double_add_t(itd, link[16],link[21],link[23],i);
		double_add_t(itd, link[19],link[22],link[24],i);
		double_div_t(itd, link[23],var[1],link[25],i);
	}

	write_dp_solution();

	return NUM_COLUMNS;
}

