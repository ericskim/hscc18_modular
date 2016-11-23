/****************************************************************************
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


long  sincosDP(iteration_data *itd, double t, double v[], double p[], int ORDER, double *cvfd)
{

	int i;
	static int   VARIABLES        = 2;
	static int   PARAMETERS       = 0;
	static int   FUNCTIONS        = 0;
	static int   LINKS            = 1;
	static int   POS_FUNCTIONS[1] = {0};

	initialize_dp_case();

	double ct[] = {-1.};

	for(i=0;  i<=ORDER; i++) {
		double_var_t(itd, var[2],var[1], i);
		double_var_t(itd, link[0],var[2], i);
		double_mul_t_cc(itd, ct[0],var[1],link[0],i);
	}

	write_dp_solution();

	return NUM_COLUMNS;
}

