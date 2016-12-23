/************************************************************************
 * 
 * Test routines for the sparse matrix BIAS functions
 * --------------------------------------------------
 *
 * Copyright (C) 2005 Christian Keil
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
 * $Id: sparseBiasTest.c 572 2008-11-04 08:01:34Z christian $
 *
 ************************************************************************/
#include <BIAS/Bias2S.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define P1 (1 - (BiasEpsilon / 2))
#define S1 (1 + BiasEpsilon)

int
int_failure(int line, INT *result, INT *shouldbe, int size)
{
	if (memcmp(result, shouldbe, size * sizeof (INT)) == 0)
		return 0;

	int i;
	INT *walk;

	printf("FAILED!!!\n");
	printf("at line %d:\n", line);

	printf("the result is\n");
	for (i = 0, walk = result; i < size; ++i, ++walk)
		printf("%i ", *walk);
	printf("\n");

	printf("but should be\n");
	for (i = 0, walk = shouldbe; i < size; ++i, ++walk)
		printf("%i ", *walk);
	printf("\n");

	return 1;
}

int
real_failure(int line, REAL *result, REAL *shouldbe, int size)
{
	if (memcmp(result, shouldbe, size * sizeof (REAL)) == 0)
		return 0;

	int i;
	REAL *walk;

	printf("FAILED!!!\n");
	printf("at line %d:\n", line);

	printf("the result is\n");
	for (i = 0, walk = result; i < size; ++i, ++walk)
		printf("%a ", *walk);
	printf("\n");

	printf("but should contain\n");
	for (i = 0, walk = shouldbe; i < size; ++i, ++walk)
		printf("%a ", *walk);
	printf("\n");

	return 1;
}

int
interval_failure(int line, BIASINTERVAL *result, BIASINTERVAL *shouldbe, int size)
{
	if (memcmp(result, shouldbe, size * sizeof (BIASINTERVAL)) == 0)
		return 0;

	int i;
	BIASINTERVAL *walk;

	printf("FAILED!!!\n");
	printf("at line %d:\n", line);

	printf("the result is\n");
	for (i = 0, walk = result; i < size; ++i, ++walk)
		printf("[%a %a]", walk->inf, walk->sup);
	printf("\n");

	printf("but should contain\n");
	for (i = 0, walk = shouldbe; i < size; ++i, ++walk)
		printf("[%a %a]", walk->inf, walk->sup);
	printf("\n");

	return 1;
}

int
main(void)
{
	printf("-= Test of BIAS's sparse matrix routines\n\n");

	BiasInit();
	/* allocating these dynamically ensures that valgrind finds more errors -
	 * otherwise they seem to be located in a rather "friendly" memory region
	 * */
	INT *emptyCS = (INT *) calloc(1, sizeof (INT));
	INT *emptyCS2 = (INT *) calloc(1, sizeof (INT));

	printf("sparse matrix functions\n");
	printf("-----------------------\n");
	{
		{
			REAL pa[] = {1.0, 1.0};
			BIASINTERVAL pA[] = {{1.0, 2.0}, {1.0, 2.0}};
			INT paARI[] = {0, 1};
			INT paACS[] = {0, 1, 2};
			INT dummy;
			INT prR1RI[] = {0, 1};
			INT prR1CS[] = {0, 1, 2};
			REAL *pr2 = (REAL *) malloc(2 * sizeof (REAL));
			BIASINTERVAL *pR2 = (BIASINTERVAL *) malloc(2 * sizeof (BIASINTERVAL));
			INT *prR2RI = (INT *) malloc(2 * sizeof (INT));
			INT *prR2CS = (INT *) malloc(3 * sizeof (INT));

			printf("BiasPredS[RI]... ");
			{
				{
					REAL pr1[] = {P1, P1};

					BiasPredSR(pr2, prR2RI, prR2CS, &dummy, pa, paARI, paACS, 2, 2);
					if (real_failure(__LINE__, pr2, pr1, 2))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, prR2RI, prR1RI, 2))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, prR2CS, prR1CS, 3))
						exit(EXIT_FAILURE);

					BiasPredSR(NULL, NULL, emptyCS, &dummy, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BIASINTERVAL pR1[] = {{S1, 2.0 - BiasEpsilon},
						{S1, 2.0 - BiasEpsilon}};
					BiasPredSI(pR2, prR2RI, prR2CS, &dummy, pA, paARI, paACS, 2, 2);

					if (interval_failure(__LINE__, pR2, pR1, 2))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, prR2RI, prR1RI, 2))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, prR2CS, prR1CS, 3))
						exit(EXIT_FAILURE);

					BiasPredSI(NULL, NULL, emptyCS, &dummy, NULL, NULL, emptyCS, 0, 0);
				}
			}
			printf("OK\n");

			printf("BiasSuccS[RI]... ");
			{
				{
					REAL pr1[] = {S1, S1};

					BiasSuccSR(pr2, prR2RI, prR2CS, &dummy, pa, paARI, paACS, 2, 2);
					if (real_failure(__LINE__, pr2, pr1, 2))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, prR2RI, prR1RI, 2))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, prR2CS, prR1CS, 3))
						exit(EXIT_FAILURE);

					BiasSuccSR(NULL, NULL, emptyCS, &dummy, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BIASINTERVAL pR1[] = {{P1, 2.0 + (BiasEpsilon * 2)},
						{P1, 2.0 + (BiasEpsilon * 2)}};

					BiasSuccSI(pR2, prR2RI, prR2CS, &dummy, pA, paARI, paACS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 2))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, prR2RI, prR1RI, 2))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, prR2CS, prR1CS, 3))
						exit(EXIT_FAILURE);

					BiasSuccSI(NULL, NULL, emptyCS, &dummy, NULL, NULL, emptyCS, 0, 0);
				}
			}
			printf("OK\n");

			printf("BiasInfS... ");
			{
				REAL pr1[] = {1.0, 1.0};

				BiasInfS(pr2, prR2RI, prR2CS, &dummy, pA, paARI, paACS, 2, 2);
				if (real_failure(__LINE__, pr2, pr1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2CS, prR1CS, 3))
					exit(EXIT_FAILURE);

				BiasInfS(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, 0);
			}
			printf("OK\n");

			printf("BiasSupS... ");
			{
				REAL pr1[] = {2.0, 2.0};

				BiasSupS(pr2, prR2RI, prR2CS, &dummy, pA, paARI, paACS, 2, 2);
				if (real_failure(__LINE__, pr2, pr1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2CS, prR1CS, 3))
					exit(EXIT_FAILURE);

				BiasSupS(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, 0);
			}
			printf("OK\n");

			printf("BiasNegS... ");
			{
				BIASINTERVAL pR1[] = {{-2.0, -1.0}, {-2.0, -1.0}};

				BiasNegS(pR2, prR2RI, prR2CS, &dummy, pA, paARI, paACS, 2, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2CS, prR1CS, 3))
					exit(EXIT_FAILURE);

				BiasNegS(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, 0);
			}
			printf("OK\n");

			printf("BiasMidS... ");
			{
				REAL pr1[] = {1.5, 1.5};

				BiasMidS(pr2, prR2RI, prR2CS, &dummy, pA, paARI, paACS, 2, 2);
				if (real_failure(__LINE__, pr2, pr1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2CS, prR1CS, 3))
					exit(EXIT_FAILURE);

				BiasMidS(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, 0);
			}
			printf("OK\n");

			printf("BiasMidRadS... ");
			{
				REAL pm1[] = {1.5, 1.5};
				REAL *pm2 = (REAL *) malloc(2 * sizeof (REAL));
				INT *pm2RI = (INT *) malloc(2 * sizeof (INT));
				INT *pm2CS = (INT *) malloc(3 * sizeof (INT));
				REAL pr1[] = {0.5, 0.5};

				BiasMidRadS(pm2, pm2RI, pm2CS, &dummy, pr2, prR2RI, prR2CS, &dummy,
										pA, paARI, paACS, 2, 2);
				if (real_failure(__LINE__, pm2, pm1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pm2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pm2CS, prR1CS, 3))
					exit(EXIT_FAILURE);
				if (real_failure(__LINE__, pr2, pr1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2CS, prR1CS, 3))
					exit(EXIT_FAILURE);

				BiasMidRadS(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS2, &dummy,
										NULL, NULL, emptyCS, 0, 0);
			}
			printf("OK\n");

			printf("BiasDiamS... ");
			{
				REAL pr1[] = {1.0, 1.0};

				BiasDiamS(pr2, prR2RI, prR2CS, &dummy, pA, paARI, paACS, 2, 2);
				if (real_failure(__LINE__, pr2, pr1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2CS, prR1CS, 3))
					exit(EXIT_FAILURE);

				BiasDiamS(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, 0);
			}
			printf("OK\n");

			printf("BiasAbsS... ");
			{
				REAL pr1[] = {2.0, 2.0};

				BiasAbsS(pr2, prR2RI, prR2CS, &dummy, pA, paARI, paACS, 2, 2);
				if (real_failure(__LINE__, pr2, pr1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2CS, prR1CS, 3))
					exit(EXIT_FAILURE);

				BiasAbsS(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, 0);
			}
			printf("OK\n");

			printf("BiasHullSR... ");
			{
				BIASINTERVAL pR1[] = {{1.0, 1.0}, {1.0, 1.0}};

				BiasHullSR(pR2, prR2RI, prR2CS, &dummy, pa, paARI, paACS, 2, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2RI, prR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, prR2CS, prR1CS, 3))
					exit(EXIT_FAILURE);

				BiasHullSR(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, 0);
			}
			printf("OK\n");
		}

		printf("BiasSetToZeroS... ");
		{
			INT dummy = 0;
			INT pR1CS[] = {0, 0, 0};
			INT pR2CS[] = {0, 0, 2};

			BiasSetToZeroS(pR2CS, &dummy, 2);
			if (int_failure(__LINE__, pR2CS, pR1CS, 3))
				exit(EXIT_FAILURE);

			BiasSetToZeroS(emptyCS, &dummy, 0);
		}
		printf("OK\n");
	}
	printf("\n");

	printf("sparse matrix x scalar functions\n");
	printf("--------------------------------\n");
	{
		REAL a = 1.0;
		BIASINTERVAL A = {1.0, 2.0};
		REAL pb[] = {1.0, 1.0};
		BIASINTERVAL pB[] = {{1.0, 2.0}, {1.0, 2.0}};
		INT pbBRI[] = {0, 1};
		INT pbBCS[] = {0, 1, 2};
		INT dummy = 0;
		INT pR1RI[] = {0, 1};
		INT pR1CS[] = {0, 1, 2};
		INT *pR2RI = (INT *) malloc(2 * sizeof (INT));
		INT *pR2CS = (INT *) malloc(3 * sizeof (INT));
		BIASINTERVAL *pR2 = (BIASINTERVAL *) malloc(2 * sizeof (BIASINTERVAL));

		printf("BiasMul[RI]S[RI]... ");
		{
			{
				BIASINTERVAL pR1[] = {{1.0, 1.0}, {1.0, 1.0}};

				BiasMulRSR(pR2, pR2RI, pR2CS, &dummy, &a, pb, pbBRI, pbBCS, 2, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 2))
					exit(EXIT_FAILURE);

				BiasMulRSR(NULL, NULL, emptyCS2, &dummy, &a, NULL, NULL, emptyCS, 0, 0);
			}
			{
				BIASINTERVAL pR1[] = {{1.0, 2.0}, {1.0, 2.0}};

				BiasMulRSI(pR2, pR2RI, pR2CS, &dummy, &a, pB, pbBRI, pbBCS, 2, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 2))
					exit(EXIT_FAILURE);

				BiasMulRSI(NULL, NULL, emptyCS2, &dummy, &a, NULL, NULL, emptyCS, 0, 0);
			}
			{
				BIASINTERVAL pR1[] = {{1.0, 2.0}, {1.0, 2.0}};

				BiasMulISR(pR2, pR2RI, pR2CS, &dummy, &A, pb, pbBRI, pbBCS, 2, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 2))
					exit(EXIT_FAILURE);

				BiasMulISR(NULL, NULL, emptyCS2, &dummy, &A, NULL, NULL, emptyCS, 0, 0);
			}
			{
				BIASINTERVAL pR1[] = {{1.0, 4.0}, {1.0, 4.0}};

				BiasMulISI(pR2, pR2RI, pR2CS, &dummy, &A, pB, pbBRI, pbBCS, 2, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 2))
					exit(EXIT_FAILURE);

				BiasMulISI(NULL, NULL, emptyCS2, &dummy, &A, NULL, NULL, emptyCS, 0, 0);
			}
		}
		printf("OK\n");

		printf("BiasDivS[RI][RI]... ");
		{
			{
				BIASINTERVAL pR1[] = {{1.0, 1.0}, {1.0, 1.0}};

				BiasDivSRR(pR2, pR2RI, pR2CS, &dummy, pb, pbBRI, pbBCS, 2, &a, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 2))
					exit(EXIT_FAILURE);

				BiasDivSRR(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, &a, 0);
			}
			{
				BIASINTERVAL pR1[] = {{0.5, 1.0}, {0.5, 1.0}};

				BiasDivSRI(pR2, pR2RI, pR2CS, &dummy, pb, pbBRI, pbBCS, 2, &A, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 2))
					exit(EXIT_FAILURE);

				BiasDivSRI(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, &A, 0);
			}
			{
				BIASINTERVAL pR1[] = {{1.0, 2.0}, {1.0, 2.0}};

				BiasDivSIR(pR2, pR2RI, pR2CS, &dummy, pB, pbBRI, pbBCS, 2, &a, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 2))
					exit(EXIT_FAILURE);

				BiasDivSIR(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, &a, 0);
			}
			{
				BIASINTERVAL pR1[] = {{0.5, 2.0}, {0.5, 2.0}};

				BiasDivSII(pR2, pR2RI, pR2CS, &dummy, pB, pbBRI, pbBCS, 2, &A, 2);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 2))
					exit(EXIT_FAILURE);

				BiasDivSII(NULL, NULL, emptyCS2, &dummy, NULL, NULL, emptyCS, 0, &A, 0);
			}
		}
		printf("OK\n");
	}
	printf("\n");

	printf("sparse matrix x vector functions\n");
	printf("--------------------------------\n");
	{
		printf("BiasMulS[RI]V[RI]... ");
		{
			/*
			 * test cases contain:
			 *   empty matrix column
			 *
			 *        /       \
			 *        | X o X |
			 *    a = |       |
			 *        | o o X |
			 *        \       /
			 */
			REAL pa[] = {1.0, BiasEpsilon / 2, 1.0};
			BIASINTERVAL pA[] = {{1.0, 1.0}, {BiasEpsilon / 2, BiasEpsilon / 2},
				{1.0, 1.0}};
			INT paARI[] = {0, 0, 1};
			INT paACS[] = {0, 1, 1, 3};
			REAL pb[] = {1.0, 1.0, 1.0};
			BIASINTERVAL pB[] = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
			BIASINTERVAL *pR2 = (BIASINTERVAL *) malloc(2 * sizeof (BIASINTERVAL));

			BIASINTERVAL pR1[] = {{1.0, S1}, {1.0, 1.0}};
			{
				BiasMulSRVR(pR2, pa, paARI, paACS, pb, 2, 3);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);

				BiasMulSRVR(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
			}
			{
				BiasMulSRVI(pR2, pa, paARI, paACS, pB, 2, 3);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);

				BiasMulSRVI(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
			}
			{
				BiasMulSIVR(pR2, pA, paARI, paACS, pb, 2, 3);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);

				BiasMulSIVR(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
			}
			{
				BiasMulSIVI(pR2, pA, paARI, paACS, pB, 2, 3);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);

				BiasMulSIVI(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
			}
		}
		printf("OK\n");
	}
	printf("\n");

	printf("sparse matrix x matrix functions\n");
	printf("--------------------------------\n");
	{
		{
			/*
			 * test cases contain:
			 *   empty matrix column
			 *
			 *       /     \
			 *       | X o |
			 *   b = |     |
			 *       | X o |
			 *       \     /
			 * */
			REAL pa[] = {1.0, 1.0, 0.0, 1.0};
			BIASINTERVAL pA[] = {{1.0, 1.0}, {1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}};
			REAL pb[] = {BiasEpsilon / 4, 1.0};
			BIASINTERVAL pB[] = {{BiasEpsilon / 4, BiasEpsilon / 4}, {1.0, 1.0}};
			INT pbBRI[] = {0, 1};
			INT pbBCS[] = {0, 2, 2};
			BIASINTERVAL *pR2 = (BIASINTERVAL *) malloc(4 * sizeof (BIASINTERVAL));

			printf("BiasAddM[RI]S[RI]... ");
			{
				BIASINTERVAL pR1[] = {{1.0, S1}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
				{
					BiasAddMRSR(pR2, pa, pb, pbBRI, pbBCS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasAddMRSR(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasAddMRSI(pR2, pa, pB, pbBRI, pbBCS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasAddMRSI(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasAddMISR(pR2, pA, pb, pbBRI, pbBCS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasAddMISR(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasAddMISI(pR2, pA, pB, pbBRI, pbBCS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasAddMISI(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
			}
			printf("OK\n");

			printf("BiasAddS[RI]M[RI]... ");
			{
				BIASINTERVAL pR1[] = {{1.0, S1}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
				{
					BiasAddSRMR(pR2, pb, pbBRI, pbBCS, pa, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasAddSRMR(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasAddSRMI(pR2, pb, pbBRI, pbBCS, pA, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasAddSRMI(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasAddSIMR(pR2, pB, pbBRI, pbBCS, pa, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasAddSIMR(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasAddSIMI(pR2, pB, pbBRI, pbBCS, pA, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasAddSIMI(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
			}
			printf("OK\n");

			printf("BiasSubM[RI]S[RI]... ");
			{
				BIASINTERVAL pR1[] = {{P1, 1.0}, {1.0, 1.0}, {-1.0, -1.0}, {1.0, 1.0}};
				{
					BiasSubMRSR(pR2, pa, pb, pbBRI, pbBCS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasSubMRSR(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasSubMRSI(pR2, pa, pB, pbBRI, pbBCS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasSubMRSI(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasSubMISR(pR2, pA, pb, pbBRI, pbBCS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasSubMISR(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasSubMISI(pR2, pA, pB, pbBRI, pbBCS, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasSubMISI(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
			}
			printf("OK\n");

			printf("BiasSubS[RI]M[RI]... ");
			{
				BIASINTERVAL pR1[] = {{-1.0, -P1}, {-1.0, -1.0}, {1.0, 1.0}, {-1.0, -1.0}};
				{
					BiasSubSRMR(pR2, pb, pbBRI, pbBCS, pa, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasSubSRMR(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasSubSRMI(pR2, pb, pbBRI, pbBCS, pA, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasSubSRMI(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasSubSIMR(pR2, pB, pbBRI, pbBCS, pa, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasSubSIMR(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasSubSIMI(pR2, pB, pbBRI, pbBCS, pA, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasSubSIMI(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
			}
			printf("OK\n");

			printf("BiasMulM[RI]S[RI]... ");
			{
				BIASINTERVAL pR1[] = {{1.0, S1}, {0.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}};
				{
					BiasMulMRSR(pR2, pa, pb, pbBRI, pbBCS, 2, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasMulMRSR(NULL, NULL, NULL, NULL, emptyCS, 0, 0, 0);
				}
				{
					BiasMulMRSI(pR2, pa, pB, pbBRI, pbBCS, 2, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasMulMRSI(NULL, NULL, NULL, NULL, emptyCS, 0, 0, 0);
				}
				{
					BiasMulMISR(pR2, pA, pb, pbBRI, pbBCS, 2, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasMulMISR(NULL, NULL, NULL, NULL, emptyCS, 0, 0, 0);
				}
				{
					BiasMulMISI(pR2, pA, pB, pbBRI, pbBCS, 2, 2, 2);
					if (interval_failure(__LINE__, pR2, pR1, 4))
						exit(EXIT_FAILURE);

					BiasMulMISI(NULL, NULL, NULL, NULL, emptyCS, 0, 0, 0);
				}
			}
			printf("OK\n");
		}

		printf("BiasMulS[RI]M[RI]... ");
		{
			/*
			 * test cases contain:
			 *   empty matrix column
			 *
			 *       /       \
			 *       | X o X |
			 *   b = |       |
			 *       | X o o |
			 *       \       /
			 * */
			REAL pa[] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0};
			BIASINTERVAL pA[] = {{1.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0},
				{1.0, 1.0}, {0.0, 0.0}};
			REAL pb[] = {1.0, 1.0, BiasEpsilon / 2};
			BIASINTERVAL pB[] = {{1.0, 1.0}, {1.0, 1.0},
				{BiasEpsilon / 2, BiasEpsilon / 2}};
			INT pbBRI[] = {0, 1, 0};
			INT pbBCS[] = {0, 2, 2, 3};
			BIASINTERVAL *pR2 = (BIASINTERVAL *) malloc(4 * sizeof (BIASINTERVAL));

			BIASINTERVAL pR1[] = {{1.0, S1}, {0.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}};
			{
				BiasMulSRMR(pR2, pb, pbBRI, pbBCS, pa, 2, 3, 2);
				if (interval_failure(__LINE__, pR2, pR1, 4))
					exit(EXIT_FAILURE);

				BiasMulSRMR(NULL, NULL, NULL, emptyCS, NULL, 0, 0, 0);
			}
			{
				BiasMulSRMI(pR2, pb, pbBRI, pbBCS, pA, 2, 3, 2);
				if (interval_failure(__LINE__, pR2, pR1, 4))
					exit(EXIT_FAILURE);

				BiasMulSRMI(NULL, NULL, NULL, emptyCS, NULL, 0, 0, 0);
			}
			{
				BiasMulSIMR(pR2, pB, pbBRI, pbBCS, pa, 2, 3, 2);
				if (interval_failure(__LINE__, pR2, pR1, 4))
					exit(EXIT_FAILURE);

				BiasMulSIMR(NULL, NULL, NULL, emptyCS, NULL, 0, 0, 0);
			}
			{
				BiasMulSIMI(pR2, pB, pbBRI, pbBCS, pA, 2, 3, 2);
				if (interval_failure(__LINE__, pR2, pR1, 4))
					exit(EXIT_FAILURE);

				BiasMulSIMI(NULL, NULL, NULL, emptyCS, NULL, 0, 0, 0);
			}
		}
		printf("OK\n");

		{
			/*
			 * test cases contain:
			 *   b does not start in 1st column
			 *   Inf(a_ij) > 0, b_ij empty
			 *   Sup(a_ij) < 0, b_ij empty
			 *   Inf(a_ij) > Inf(b_ij), b_ij nonempty
			 *   Sup(a_ij) < Sup(b_ij), b_ij nonempty
			 *   Inf(a_ij) <= Inf(b_ij) <= Sup(b_ij) <= Sup(a_ij), b_ij nonempty
			 *   empty column in b after nonempty one
			 *   elements in b ending column - not ending column
			 *
			 *       /         \
			 *       | o X X o |
			 *   b = |         |
			 *       | o X o o |
			 *       \         /
			 */
			REAL pa[] = {1.0, 2.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0};
			BIASINTERVAL pA[] = {{1.0, 1.0}, {2.0, 2.0}, {1.0, 1.0}, {0.0, 0.0},
				{-1.0, -1.0}, {1.0, 1.0}, {0.0, 0.0}, {0.0, 0.0}};
			REAL pb[] = {1.0, 2.0, 1.0};
			BIASINTERVAL pB[] = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 1.0}};
			INT pbBRI[] = {0, 1, 0};
			INT pbBCS[] = {0, 0, 2, 3, 3};
			BIASINTERVAL pR1[] = {{0.0, 1.0}, {1.0, 2.0}, {1.0, 1.0}, {0.0, 0.0},
				{-1.0, 0.0}, {1.0, 2.0}, {0.0, 0.0}, {0.0, 0.0}};
			BIASINTERVAL *pR2 = (BIASINTERVAL *) malloc(8 * sizeof (BIASINTERVAL));

			printf("BiasHullM[RI]S[RI]... ");
			{
				{
					BiasHullMRSR(pR2, pa, pb, pbBRI, pbBCS, 2, 4);
					if (interval_failure(__LINE__, pR2, pR1, 8))
						exit(EXIT_FAILURE);

					BiasHullMRSR(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasHullMRSI(pR2, pa, pB, pbBRI, pbBCS, 2, 4);
					if (interval_failure(__LINE__, pR2, pR1, 8))
						exit(EXIT_FAILURE);

					BiasHullMRSI(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasHullMISR(pR2, pA, pb, pbBRI, pbBCS, 2, 4);
					if (interval_failure(__LINE__, pR2, pR1, 8))
						exit(EXIT_FAILURE);

					BiasHullMISR(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
				{
					BiasHullMISI(pR2, pA, pB, pbBRI, pbBCS, 2, 4);
					if (interval_failure(__LINE__, pR2, pR1, 8))
						exit(EXIT_FAILURE);

					BiasHullMISI(NULL, NULL, NULL, NULL, emptyCS, 0, 0);
				}
			}
			printf("OK\n");

			printf("BiasHullS[RI]M[RI]... ");
			{
				{
					BiasHullSRMR(pR2, pb, pbBRI, pbBCS, pa, 2, 4);
					if (interval_failure(__LINE__, pR2, pR1, 8))
						exit(EXIT_FAILURE);

					BiasHullSRMR(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasHullSRMI(pR2, pb, pbBRI, pbBCS, pA, 2, 4);
					if (interval_failure(__LINE__, pR2, pR1, 8))
						exit(EXIT_FAILURE);

					BiasHullSRMI(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasHullSIMR(pR2, pB, pbBRI, pbBCS, pa, 2, 4);
					if (interval_failure(__LINE__, pR2, pR1, 8))
						exit(EXIT_FAILURE);

					BiasHullSIMR(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					BiasHullSIMI(pR2, pB, pbBRI, pbBCS, pA, 2, 4);
					if (interval_failure(__LINE__, pR2, pR1, 8))
						exit(EXIT_FAILURE);

					BiasHullSIMI(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
				}
			}
			printf("OK\n");
		}

		{
			/*
			 * test cases contain:
			 *   sparse matrix starting with / containing empty column
			 *   empty sparse element succeding / failing test
			 *   nonempty sparse element succeding / failing test
			 *
			 *   s2 = ( o X X )   s3 = ( X o X )
			 */
			REAL pd1[] = {1.0};
			BIASINTERVAL pD1[] = {{1.0, 1.0}};
			REAL ps1[] = {2.0};
			BIASINTERVAL pS1[] = {{2.0, 2.0}};
			INT psS1RI[] = {0};
			INT psS1CS[] = {0, 1};

			REAL pd2[] = {1.0, 0.0, 1.0};
			REAL ps23[] = {1.0, 1.0};
			BIASINTERVAL pS23[] = {{1.0, 1.0}, {1.0, 1.0}};
			INT psS23RI[] = {0, 0};
			INT psS2CS[] = {0, 0, 1, 2};
			INT psS3CS[] = {0, 1, 1, 2};

			{
				BIASINTERVAL pD2[] = {{1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}};

				printf("BiasInM[RI]S... ");
				{
					{
						if (BiasInMRS(pd1, pS1, psS1RI, psS1CS, 1, 1) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasInMRS(pd2, pS23, psS23RI, psS2CS, 1, 3) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasInMRS(pd2, pS23, psS23RI, psS3CS, 1, 3) != 1) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}

						BiasInMRS(NULL, NULL, NULL, emptyCS, 0, 0);
					}
					{
						if (BiasInMIS(pD1, pS1, psS1RI, psS1CS, 1, 1) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasInMIS(pD2, pS23, psS23RI, psS2CS, 1, 3) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasInMIS(pD2, pS23, psS23RI, psS3CS, 1, 3) != 1) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}

						BiasInMIS(NULL, NULL, NULL, emptyCS, 0, 0);
					}
				}
				printf("OK\n");

				printf("BiasInS[RI]... ");
				{
					{
						if (BiasInSR(ps1, psS1RI, psS1CS, pD1, 1, 1) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasInSR(ps23, psS23RI, psS2CS, pD2, 1, 3) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasInSR(ps23, psS23RI, psS3CS, pD2, 1, 3) != 1) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}

						BiasInSR(NULL, NULL, emptyCS, NULL, 0, 0);
					}
					{
						if (BiasInSI(pS1, psS1RI, psS1CS, pD1, 1, 1) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasInSI(pS23, psS23RI, psS2CS, pD2, 1, 3) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasInSI(pS23, psS23RI, psS3CS, pD2, 1, 3) != 1) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}

						BiasInSI(NULL, NULL, emptyCS, NULL, 0, 0);
					}
				}
				printf("OK\n");

				printf("BiasIsEqual[MS][SM]... ");
				{
					{
						if (BiasIsEqualMS(pD1, pS1, psS1RI, psS1CS, 1, 1) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasIsEqualMS(pD2, pS23, psS23RI, psS2CS, 1, 3) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasIsEqualMS(pD2, pS23, psS23RI, psS3CS, 1, 3) != 1) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}

						BiasIsEqualMS(NULL, NULL, NULL, emptyCS, 0, 0);
					}
					{
						if (BiasIsEqualSM(pS1, psS1RI, psS1CS, pD1, 1, 1) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasIsEqualSM(pS23, psS23RI, psS2CS, pD2, 1, 3) != 0) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}
						if (BiasIsEqualSM(pS23, psS23RI, psS3CS, pD2, 1, 3) != 1) {
							printf("FAILED!!!\nat line %d\n", __LINE__);
							exit(EXIT_FAILURE);
						}

						BiasIsEqualSM(NULL, NULL, emptyCS, NULL, 0, 0);
					}
				}
				printf("OK\n");
			}

			printf("BiasInInteriorS[RI]... ");
			{
				BIASINTERVAL pD2[] = {{0.5, 1.5}, {-0.5, 0.5}, {0.5, 1.5}};

				{
					if (BiasInInteriorSR(ps1, psS1RI, psS1CS, pD1, 1, 1) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasInInteriorSR(ps23, psS23RI, psS2CS, pD2, 1, 3) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasInInteriorSR(ps23, psS23RI, psS3CS, pD2, 1, 3) != 1) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}

					BiasInInteriorSR(NULL, NULL, emptyCS, NULL, 0, 0);
				}
				{
					if (BiasInInteriorSI(pS1, psS1RI, psS1CS, pD1, 1, 1) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasInInteriorSI(pS23, psS23RI, psS2CS, pD2, 1, 3) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasInInteriorSI(pS23, psS23RI, psS3CS, pD2, 1, 3) != 1) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}

					BiasInInteriorSI(NULL, NULL, emptyCS, NULL, 0, 0);
				}
			}
			printf("OK\n");
		}
	}
	printf("\n");

	printf("sparse matrix x sparse matrix functions\n");
	printf("---------------------------------------\n");
	{
		{
			/*
			 * test cases contain:
			 *   a has an empty column
			 *   a starts in column 2, b in column 3
			 *   elements only in a / b 
			 *     not ending a / b column
			 *     ending a / b column
			 *       next column empty
			 *         not ending r column
			 *         ending r column with 0, 1 empty following column
			 *       next column nonempty
			 *         not ending r column
			 *         ending r column
			 *   elements in both matrices
			 *     evaluating to nonzero
			 *     evaluating to zero
			 *     ending no column
			 *     ending a, b column 
			 *       next a column nonempty, b column empty
			 *
			 *       /                                    \
			 *       | o  X  o  X  o  o  o  X  X  o  X  X |
			 *   a = |                                    |
			 *       | o  X  o  o  o  X  o  o  o  X  X  o |
			 *       \                                    /
			 *
			 *       /                                    \
			 *       | o  o  X  o  o  X  o  o  o  X  X  o |
			 *   b = |                                    |
			 *       | o  o  X  X  o  o  o  X  o  o -X  o |
			 *       \                                    /
			 */
			INT paARI[] = {0, 1, 0, 1, 0, 0, 1, 0, 1, 0};
			INT paACS[] = {0, 0, 2, 2, 3, 3, 4, 4, 5, 6, 7, 9, 10};
			INT pbBRI[] = {0, 1, 1, 0, 1, 0, 0, 1};
			INT pbBCS[] = {0, 0, 0, 2, 3, 3, 4, 4, 5, 5, 6, 8, 8};
			INT dummy;
			BIASINTERVAL *pR2 = (BIASINTERVAL *) malloc(16 * sizeof (BIASINTERVAL));
			INT *pR2RI = (INT *) malloc(16 * sizeof (INT));
			INT *pR2CS = (INT *) malloc(13 * sizeof (INT));

			printf("BiasAddS[RI]S[RI]... ");
			{
				REAL pa[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
				BIASINTERVAL pA[] = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0},
					{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}}; 
				REAL pb[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, BiasEpsilon / 2, -1.0};
				BIASINTERVAL pB[] = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0},
					{1.0, 2.0}, {1.0, 2.0}, {BiasEpsilon / 2, BiasEpsilon / 2},
					{-1.0, -1.0}}; 
				INT pR1RI[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0};
				INT pR1CS[] = {0, 0, 2, 4, 6, 6, 8, 8, 10, 11, 13, 14, 15};

				{
					BIASINTERVAL pR1[] = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0},
						{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0},
						{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, S1}, {1.0, 1.0}};

					BiasAddSRSR(pR2, pR2RI, pR2CS, &dummy, pa, paARI, paACS, pb, pbBRI, pbBCS,
											12);
					if (interval_failure(__LINE__, pR2, pR1, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasAddSRSR(NULL, NULL, emptyCS2, &dummy,
											NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BIASINTERVAL pR1[] = {{1.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}, {1.0, 2.0},
						{1.0, 1.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 2.0},
						{1.0, 1.0}, {1.0, 2.0}, {1.0, 1.0}, {1.0, S1}, {1.0, 1.0}};

					BiasAddSRSI(pR2, pR2RI, pR2CS, &dummy, pa, paARI, paACS, pB, pbBRI, pbBCS,
											12);
					if (interval_failure(__LINE__, pR2, pR1, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasAddSRSI(NULL, NULL, emptyCS2, &dummy,
											NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BIASINTERVAL pR1[] = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 1.0}, {1.0, 1.0},
						{1.0, 2.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 1.0},
						{1.0, 2.0}, {1.0, 1.0}, {1.0, 2.0}, {1.0, S1}, {1.0, 2.0}};

					BiasAddSISR(pR2, pR2RI, pR2CS, &dummy, pA, paARI, paACS, pb, pbBRI, pbBCS,
											12);
					if (interval_failure(__LINE__, pR2, pR1, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasAddSISR(NULL, NULL, emptyCS2, &dummy,
											NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BIASINTERVAL pR1[] = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0},
						{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0},
						{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, S1}, {1.0, 2.0}};

					BiasAddSISI(pR2, pR2RI, pR2CS, &dummy, pA, paARI, paACS, pB, pbBRI, pbBCS,
											12);
					if (interval_failure(__LINE__, pR2, pR1, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasAddSISI(NULL, NULL, emptyCS2, &dummy,
											NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
			}
			printf("OK\n");

			printf("BiasSubS[RI]S[RI]... ");
			{
				REAL pa[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
				BIASINTERVAL pA[] = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0},
					{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 2.0}}; 
				REAL pb[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -(BiasEpsilon / 2)};
				BIASINTERVAL pB[] = {{1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0}, {1.0, 2.0},
					{1.0, 2.0}, {1.0, 2.0}, {1.0, 1.0},
					{-(BiasEpsilon / 2), -(BiasEpsilon / 2)}}; 
				INT pR1RI[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0};
				INT pR1CS[] = {0, 0, 2, 4, 6, 6, 8, 8, 10, 11, 13, 14, 15};

				{
					BIASINTERVAL pR1[] = {{1.0, 1.0}, {1.0, 1.0}, {-1.0, -1.0}, {-1.0, -1.0},
						{1.0, 1.0}, {-1.0, -1.0}, {-1.0, -1.0}, {1.0, 1.0}, {1.0, 1.0}, {-1.0, -1.0},
						{1.0, 1.0}, {-1.0, -1.0}, {1.0, 1.0}, {1.0, S1}, {1.0, 1.0}};

					BiasSubSRSR(pR2, pR2RI, pR2CS, &dummy, pa, paARI, paACS, pb, pbBRI, pbBCS,
											12);
					if (interval_failure(__LINE__, pR2, pR1, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasSubSRSR(NULL, NULL, emptyCS2, &dummy,
											NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BIASINTERVAL pR1[] = {{1.0, 1.0}, {1.0, 1.0}, {-2.0, -1.0}, {-2.0, -1.0},
						{1.0, 1.0}, {-2.0, -1.0}, {-2.0, -1.0}, {1.0, 1.0}, {1.0, 1.0},
						{-2.0, -1.0}, {1.0, 1.0}, {-2.0, -1.0}, {1.0, 1.0}, {1.0, S1},
						{1.0, 1.0}};

					BiasSubSRSI(pR2, pR2RI, pR2CS, &dummy, pa, paARI, paACS, pB, pbBRI, pbBCS,
											12);
					if (interval_failure(__LINE__, pR2, pR1, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasSubSRSI(NULL, NULL, emptyCS2, &dummy,
											NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BIASINTERVAL pR1[] = {{1.0, 2.0}, {1.0, 2.0}, {-1.0, -1.0}, {-1.0, -1.0},
						{1.0, 2.0}, {-1.0, -1.0}, {-1.0, -1.0}, {1.0, 2.0}, {1.0, 2.0}, {-1.0, -1.0},
						{1.0, 2.0}, {-1.0, -1.0}, {1.0, 2.0}, {1.0, S1}, {1.0, 2.0}};

					BiasSubSISR(pR2, pR2RI, pR2CS, &dummy, pA, paARI, paACS, pb, pbBRI, pbBCS,
											12);
					if (interval_failure(__LINE__, pR2, pR1, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasSubSISR(NULL, NULL, emptyCS2, &dummy,
											NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BIASINTERVAL pR1[] = {{1.0, 2.0}, {1.0, 2.0}, {-2.0, -1.0}, {-2.0, -1.0},
						{1.0, 2.0}, {-2.0, -1.0}, {-2.0, -1.0}, {1.0, 2.0}, {1.0, 2.0},
						{-2.0, -1.0}, {1.0, 2.0}, {-2.0, -1.0}, {1.0, 2.0}, {1.0, S1},
						{1.0, 2.0}};

					BiasSubSISI(pR2, pR2RI, pR2CS, &dummy, pA, paARI, paACS, pB, pbBRI, pbBCS,
											12);
					if (interval_failure(__LINE__, pR2, pR1, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 15))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasSubSISI(NULL, NULL, emptyCS2, &dummy,
											NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
			}
			printf("OK\n");

			printf("BiasIntersectionS... ");
			{
				{
					/*
					 * test cases contain:
					 *   elements only in a, b, both matrices not intersecting
					 */
					BIASINTERVAL pA1[] = {{1.0, 1.0}};
					BIASINTERVAL pA2[] = {};
					INT pA1RI[] = {0};
					INT pA1CS[] = {0,1,1};
					INT pA2RI[] = {};
					INT pA2CS[] = {0, 0, 0};
					BIASINTERVAL pB1[] = {{2.0, 2.0}};
					BIASINTERVAL pB2[] = {};
					INT pB1RI[] = {0};
					INT pB1CS[] = {0,1,1};
					INT pB2RI[] = {};
					INT pB2CS[] = {0, 0, 0};

					if (BiasIntersectionS(pR2, pR2RI, pR2CS, &dummy,
																pA1, pA1RI, pA1CS, pB2, pB2RI, pB2CS, 2)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasIntersectionS(pR2, pR2RI, pR2CS, &dummy,
																pA2, pA2RI, pA2CS, pB1, pB1RI, pB1CS, 2)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasIntersectionS(pR2, pR2RI, pR2CS, &dummy,
																pA1, pA1RI, pA1CS, pB1, pB1RI, pB1CS, 2)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
				}
				{
					/*
					 * test cases see AddSRSR
					 */
					BIASINTERVAL pA[] = {{-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
						{-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 0.0}, {-1.0, 1.0}}; 
					BIASINTERVAL pB[] = {{-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
						{-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}, {0.0, 1.0}}; 
					BIASINTERVAL pR1[] = {{-1.0, 1.0}};
					INT pR1RI[] = {0};
					INT pR1CS[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1};

					if (BiasIntersectionS(pR2, pR2RI, pR2CS, &dummy,
																pA, paARI, paACS, pB, pbBRI, pbBCS, 12) != 1)
					{
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (interval_failure(__LINE__, pR2, pR1, 1))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 1))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);
				}

				BiasIntersectionS(NULL, NULL, emptyCS2, &dummy,
													NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
			}
			printf("OK\n");

			printf("BiasHullS[RI]S[RI]... ");
			{
				/* 
				 * test cases contain:
				 *   see AddSRSR
				 *   elements only in a, b with Inf > 0, Sup < 0
				 *
				 *       /                                      \
				 *       | o  X>0  o  X  o  o  o  X  X  o  X  X |
				 *   a = |                                      |
				 *       | o  X<0  o  o  o  X  o  o  o  X  X  o |
				 *       \                                      /
				 *
				 *       /                                      \
				 *       | o  o  X>0  o  o  X  o  o  o  X  X  o |
				 *   b = |                                      |
				 *       | o  o  X<0  X  o  o  o  X  o  o  X  o |
				 *       \                                      /
				 */
				REAL pa[] = {1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
				BIASINTERVAL pA[] = {{1.0, 1.0}, {-1.0, -1.0}, {1.0, 1.0}, {1.0, 1.0},
					{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0},
					{1.0, 1.0}};
				REAL pb[] = {1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
				BIASINTERVAL pB[] = {{1.0, 1.0}, {-1.0, -1.0}, {1.0, 1.0}, {1.0, 1.0},
					{1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}, {1.0, 1.0}};
				BIASINTERVAL pR1[] = {{0.0, 1.0}, {-1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0},
					{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0},
					{0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {0.0, 1.0}, {1.0, 1.0},
					{1.0, 1.0}, {0.0, 1.0}};
				INT pR1RI[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0};
				INT pR1CS[] = {0, 0, 2, 4, 6, 6, 8, 8, 10, 11, 13, 15, 16};

				{
					BiasHullSRSR(pR2, pR2RI, pR2CS, &dummy,
											 pa, paARI, paACS, pb, pbBRI, pbBCS, 12);
					if (interval_failure(__LINE__, pR2, pR1, 16))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 16))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasHullSRSR(NULL, NULL, emptyCS2, &dummy,
											 NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BiasHullSRSI(pR2, pR2RI, pR2CS, &dummy,
											 pa, paARI, paACS, pB, pbBRI, pbBCS, 12);
					if (interval_failure(__LINE__, pR2, pR1, 16))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 16))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasHullSRSI(NULL, NULL, emptyCS2, &dummy,
											 NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BiasHullSISR(pR2, pR2RI, pR2CS, &dummy,
											 pA, paARI, paACS, pb, pbBRI, pbBCS, 12);
					if (interval_failure(__LINE__, pR2, pR1, 16))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 16))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasHullSISR(NULL, NULL, emptyCS2, &dummy,
											 NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				{
					BiasHullSISI(pR2, pR2RI, pR2CS, &dummy,
											 pA, paARI, paACS, pB, pbBRI, pbBCS, 12);
					if (interval_failure(__LINE__, pR2, pR1, 16))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2RI, pR1RI, 16))
						exit(EXIT_FAILURE);
					if (int_failure(__LINE__, pR2CS, pR1CS, 13))
						exit(EXIT_FAILURE);

					BiasHullSISI(NULL, NULL, emptyCS2, &dummy,
											 NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
			}
			printf("OK\n");
		}

		printf("BiasInS[RI]S... ");
		{
			{
				/*
				 * test cases contain:
				 *   empty elements failing
				 *   nonempty elements failing
				 */
				REAL ps1[] = {};
				REAL ps2[] = {1.0};
				INT ps1RI[] = {};
				INT ps1CS[] = {0, 0};
				INT ps2RI[] = {0};
				INT ps2CS[] = {0, 1};

				BIASINTERVAL pS1[] = {};
				BIASINTERVAL pS2[] = {{2.0, 2.0}};
				BIASINTERVAL pS3[] = {{3.0, 3.0}};
				INT pS1RI[] = {};
				INT pS1CS[] = {0, 0};
				INT pS23RI[] = {0};
				INT pS23CS[] = {0, 1};

				{
					if (BiasInSRS(ps1, ps1RI, ps1CS, pS2, pS23RI, pS23CS, 1) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasInSRS(ps2, ps2RI, ps2CS, pS1, pS1RI, pS1CS, 1) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasInSRS(ps2, ps2RI, ps2CS, pS2, pS23RI, pS23CS, 1) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
				}
				{
					if (BiasInSIS(pS1, pS1RI, pS1CS, pS2, pS23RI, pS23CS, 1) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasInSIS(pS2, pS23RI, pS23CS, pS1, pS1RI, pS1CS, 1) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (BiasInSIS(pS2, pS23RI, pS23CS, pS3, pS23RI, pS23CS, 1) != 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
				}
			}
			{
				/*
				 * test cases contain:
				 *   a starts in column 4, b in column 2
				 *   elements only in b
				 *     not ending b column
				 *     ending b column
				 *       next column empty
				 *       next column nonempty
				 *   elements in both matrices
				 *     ending no column
				 *     ending a/b column
				 *       next a column empty, b column nonempty
				 *
				 *       /             \       /             \
				 *       | o o o o X o |       | o X o X X X |
				 *   a = |             |   b = |             |
				 *       | o o o o X o |       | o X o o X o |
				 *       \             /       \             /
				 */
				REAL pa[] = {1.0, 1.0};
				BIASINTERVAL pA[] = {{0.0, 1.0}, {0.0, 1.0}};
				INT paARI[] = {0, 1};
				INT paACS[] = {0, 0, 0, 0, 0, 2, 2};

				BIASINTERVAL pB[] = {{-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0},
					{-1.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0}};
				INT pBRI[] = {0, 1, 0, 0, 1, 0};
				INT pBCS[] = {0, 0, 2, 2, 3, 5, 6};

				{
					if (BiasInSRS(pa, paARI, paACS, pB, pBRI, pBCS, 6) != 1) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
				}
				{
					if (BiasInSIS(pA, paARI, paACS, pB, pBRI, pBCS, 6) != 1) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
				}
			}

			BiasInSRS(NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
			BiasInSIS(NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
		}
		printf("OK\n");

		printf("BiasMulS[RI]S[RI]... ");
		{
			/*
			 * test cases contain:
			 *   empty columns in a, b
			 *   even and odd number of elements in b column 
			 *     -> short circuit because of rounding mode switch optimization
			 *   0 elements in r
			 *
			 *       /       \
			 *       | X X o |
			 *       |       |
			 *   a = | o o o |
			 *       |       |
			 *       | o o o |
			 *       \       /
			 *       
			 *       /       \
			 *       | o X X |
			 *       |       |
			 *   b = | o o X |
			 *       |       |
			 *       | o o o |
			 *       \       /
			 */
			REAL pa[] = {0.5, 0.5};
			BIASINTERVAL pA[] = {{0.5, 0.5}, {0.5, 0.5}};
			INT paARI[] = {0, 0};
			INT paACS[] = {0, 1, 2, 2};
			REAL pb[] = {BiasEta, 2.0, BiasEpsilon};
			BIASINTERVAL pB[] = {{BiasEta, BiasEta}, {2.0, 2.0},
				{BiasEpsilon, BiasEpsilon}};
			INT pbBRI[] = {0, 0, 1};
			INT pbBCS[] = {0, 0, 1, 3};
			INT pR1RI[] = {0, 0};
			INT pR1CS[] = {0, 0, 1, 2};
			INT dummy;
			BIASINTERVAL *pR2 = (BIASINTERVAL *) malloc(2 * sizeof (BIASINTERVAL));
			INT *pR2RI = (INT *) malloc(2 * sizeof (INT));
			INT *pR2CS = (INT *) malloc(4 * sizeof (INT));

			BIASINTERVAL pR1[] = {{0.0, BiasEta}, {1.0, S1}};
			{
				BiasMulSRSR(pR2, pR2RI, pR2CS, &dummy,
										pa, paARI, paACS, pb, pbBRI, pbBCS,
										3, 3);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 4))
					exit(EXIT_FAILURE);

				BiasMulSRSR(NULL, NULL, emptyCS2, &dummy,
										NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0, 0);
			}
			{
				BiasMulSRSI(pR2, pR2RI, pR2CS, &dummy,
										pa, paARI, paACS, pB, pbBRI, pbBCS,
										3, 3);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 4))
					exit(EXIT_FAILURE);

				BiasMulSRSI(NULL, NULL, emptyCS2, &dummy,
										NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0, 0);
			}
			{
				BiasMulSISR(pR2, pR2RI, pR2CS, &dummy,
										pA, paARI, paACS, pb, pbBRI, pbBCS,
										3, 3);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 4))
					exit(EXIT_FAILURE);

				BiasMulSISR(NULL, NULL, emptyCS2, &dummy,
										NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0, 0);
			}
			{
				BiasMulSISI(pR2, pR2RI, pR2CS, &dummy,
										pA, paARI, paACS, pB, pbBRI, pbBCS,
										3, 3);
				if (interval_failure(__LINE__, pR2, pR1, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2RI, pR1RI, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pR2CS, pR1CS, 4))
					exit(EXIT_FAILURE);

				BiasMulSISI(NULL, NULL, emptyCS2, &dummy,
										NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0, 0);
			}
		}
		printf("OK\n");

		printf("BiasIsEqualS... ");
		{
			{
				/*
				 * test cases contain:
				 *   empty elements failing
				 *   nonemtpy elements failing
				 */
				BIASINTERVAL pS1[] = {};
				INT pS1RI[] = {};
				INT pS1CS[] = {0, 0};
				BIASINTERVAL pS2[] = {{1.0, 1.0}};
				INT pS2RI[] = {0};
				INT pS2CS[] = {0, 1};

				if (BiasIsEqualS(pS1, pS1RI, pS1CS, pS2, pS2RI, pS2CS, 1) != 0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (BiasIsEqualS(pS2, pS2RI, pS2CS, pS1, pS1RI, pS1CS, 1) != 0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
			{
				/*
				 * test cases contain:
				 *   elements in both matrices
				 *     ending no column
				 *     ending both columns
				 */
				BIASINTERVAL pS1[] = {{1.0, 1.0}, {1.0, 1.0}};
				BIASINTERVAL pS2[] = {{1.0, 1.0}, {1.0, 1.0}};
				INT pS12RI[] = {0, 1};
				INT pS12CS[] = {0, 2};

				if (BiasIsEqualS(pS1, pS12RI, pS12CS, pS2, pS12RI, pS12CS, 1) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}

			BiasIsEqualS(NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
		}
		printf("OK\n");
	}

	free(emptyCS);
	free(emptyCS2);
	printf("-= Passed\n");

	exit(EXIT_SUCCESS);
}
