/************************************************************************
 *
 * Basic Interval Arithmetic Subroutines Standard Functions Test Program
 * ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 Christian Keil
 *
 * This file is part of PROFIL/BIAS.
 *
 * PROFIL/BIAS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
 * USA.
 *
 * $Id $
 *
 ************************************************************************/

/*
 * For the moment being this is no exhaustive test, but just a test for known problems at some
 * special x-values (e.g., exp(1) with 64bit and GNU math library).
 */

#include <BiasRnd.h>
#include <BIAS/BiasF.h>
#include <stdio.h>

char
Test(char *sz, VOID (*f)(BIASINTERVAL * const, const BIASINTERVAL * const), REAL x, REAL y)
{
	BIASINTERVAL X, R1, R2, R3;
	char st = 0;

	X.inf = X.sup = x;
	BiasRoundDown();
	f(&R1, &X);
	BiasRoundNear();
	f(&R2, &X);
	BiasRoundUp();
	f(&R3, &X);

	printf("%s(%e = %a) = %e = %a\n", sz, x, x, y, y);
	printf(" D: [%e, %e] = [%a, %a]", R1.inf, R1.sup, R1.inf, R1.sup);
	if (!BiasInR(&y, &R1)) {
		st = 2;
		printf(" !");
	}
	printf("\n N: [%e, %e] = [%a, %a]", R2.inf, R2.sup, R2.inf, R2.sup);
	if (!BiasInR(&y, &R2)) {
		st = 1;
		printf(" !");
	}
	printf("\n U: [%e, %e] = [%a, %a]", R3.inf, R3.sup, R3.inf, R3.sup);
	if (!BiasInR(&y, &R3)) {
		st = 2;
		printf(" !");
	}
	if (st == 2) {
#ifdef __BIASSETROUNDTONEAREST__
		printf("\n*** WARNING: Only valid when rounding to nearest!\n");
		printf("*** Please make sure that you don't change rounding mode prior calling standard "
					 "functions.");
#else
		st = 1;
#endif
	}

	printf("\n");

	return st;
}

int
main(int argc, char *argv[])
{
	char st = 0;

	BiasFuncInit();

	printf("-= Test for known standard functions problems\n\n");

	/* Testing for exp(1) bug (e.g., in 64Bit GNU math library). */
	st |= Test("BiasExp", BiasExp, 1.0, BiasE);

	if (st & 1) {
		printf("\nFAILED!\n");
#ifndef __BIASSETROUNDTONEAREST__
		printf("*** Please try to define __BIASSETROUNDTONEAREST__ in "
					 "config/(your architecture)/BiasRnd.h\n");
#endif
		return st;
	} else if (st & 2) {
		printf("\n-= Passed with WARNING. Please check.\n");
		return st;
	}

	printf("\n-= Passed\n");
	return 0;
}
