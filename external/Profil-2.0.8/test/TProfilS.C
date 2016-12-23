/************************************************************************
 * 
 * Test routines for sparse PROFIL operaions
 * -----------------------------------------
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
 * $Id: sparseProfilTest.C 572 2008-11-04 08:01:34Z christian $
 * 
 ************************************************************************/
#include <SLSS.h>
#include <SparseUtilities.h>
#include <SparseIntervalMatrix.h>
#include <SparseMatrix.h>
#include <SparseRealOp.h>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
using namespace std;

/* If the sparse constructors clear the matrix, i.e. set the 'colStarts' array
 * to zero we don't need to do this ourselves. Otherwise the constructor saves
 * the memset making it necessary to call 'Clear(...)' prior to almost all other
 * calls if we don't assign to the matrix, which initializes the 'colStart'
 * array.
 *
 * To additionally clear the matrix define 'SAFETYCLEAR(x)' to be 'Clear(x)'
 * */
#define SAFETYCLEAR(x)

bool
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

	printf("but should be\n");
	for (i = 0, walk = shouldbe; i < size; ++i, ++walk)
		printf("%a ", *walk);
	printf("\n");

	return 1;
}

bool
interval_failure(int line, BIASINTERVAL *result, BIASINTERVAL *shouldbe, int
								 size)
{
	if (memcmp(result, shouldbe, size * sizeof (BIASINTERVAL)) == 0)
		return false;

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

	return true;
}

int
main(int argc, char *argv[])
{
	printf("-= Test of PROFIL's sparse matrix routines\n\n");

	printf("SparseRealOp\n");
	printf("------------\n");
	{
		INT emptyCS[] = {0};
		{
			/*
			 * test cases contain:
			 *   empty matrices
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
			REAL pa[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
			INT paRI[] = {0, 1, 0, 1, 0, 0, 1, 0, 1, 0};
			INT paCS[] = {0, 0, 2, 2, 3, 3, 4, 4, 5, 6, 7, 9, 10};
			REAL pb[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0};
			INT pbRI[] = {0, 1, 1, 0, 1, 0, 0, 1};
			INT pbCS[] = {0, 0, 0, 2, 3, 3, 4, 4, 5, 5, 6, 8, 8};
			INT pr1CS[] = {0, 0, 2, 4, 6, 6, 8, 8, 10, 11, 13, 14, 15};
			REAL *pr2  = (REAL *) malloc(15 * sizeof (REAL));
			INT *pr2RI =  (INT *) malloc(15 * sizeof (INT));
			INT *pr2CS =  (INT *) malloc(13 * sizeof (INT));

			printf("SparseRealOpAdd... ");
			{
				REAL pr1[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
					1.0, 2.0, 1.0};
				INT pr1RI[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0};

				SparseRealOpAdd(pr2, pr2RI, pr2CS, pa, paRI, paCS, pb, pbRI, pbCS, 12);
				if (real_failure(__LINE__, pr2, pr1, 15))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pr2RI, pr1RI, 15))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pr2CS, pr1CS, 13))
					exit(EXIT_FAILURE);

				SparseRealOpAdd(NULL, NULL, emptyCS,
												NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
			}
			printf("OK\n");

			printf("SparseRealOpSub... ");
			{
				REAL pr1[] = {1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0,
					-1.0, 1.0, 2.0, 1.0};
				INT pr1RI[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0};

				SparseRealOpSub(pr2, pr2RI, pr2CS, pa, paRI, paCS, pb, pbRI, pbCS, 12);
				if (real_failure(__LINE__, pr2, pr1, 15))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pr2RI, pr1RI, 15))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, pr2CS, pr1CS, 13))
					exit(EXIT_FAILURE);

				SparseRealOpSub(NULL, NULL, emptyCS,
												NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
			}
			printf("OK\n");

			{
				/*
				 * test cases contain:
				 *   elements only in a/b failing
				 *   elements in both matrices failing
				 *   test cases from surrounding block
				 */
				REAL ps1[] = {-1.0};
				REAL ps2[] = {};
				REAL ps3[] = {1.0};
				INT ps13RI[] = {0};
				INT ps2RI[] = {1};
				INT ps13CS[] = {0, 1};
				INT ps2CS[] = {0, 0};

				REAL pa[] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
				REAL pb[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0};

				printf("SparseRealOpLessThan... ");
				{
					if (SparseRealOpLessThan(ps2, ps2RI, ps2CS, ps1, ps13RI, ps13CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpLessThan(ps3, ps13RI, ps13CS, ps2, ps2RI, ps2CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpLessThan(ps3, ps13RI, ps13CS, ps1, ps13RI, ps13CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpLessThan(pa, paRI, paCS, pb, pbRI, pbCS, 12) != 1) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}

					SparseRealOpLessThan(NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				printf("OK\n");

				printf("SparseRealOpLessEqual... ");
				{
					if (SparseRealOpLessEqual(ps2, ps2RI, ps2CS, ps1, ps13RI, ps13CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpLessEqual(ps3, ps13RI, ps13CS, ps2, ps2RI, ps2CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpLessEqual(ps3, ps13RI, ps13CS, ps1, ps13RI, ps13CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpLessEqual(pa, paRI, paCS, pb, pbRI, pbCS, 12) != 1) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}

					SparseRealOpLessEqual(NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				printf("OK\n");

				printf("SparseRealOpGreaterThan... ");
				{
					if (SparseRealOpGreaterThan(ps1, ps13RI, ps13CS, ps2, ps2RI, ps2CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpGreaterThan(ps2, ps2RI, ps2CS, ps3, ps13RI, ps13CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpGreaterThan(ps1, ps13RI, ps13CS, ps3, ps13RI, ps13CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpGreaterThan(pb, pbRI, pbCS, pa, paRI, paCS, 12) != 1) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}

					SparseRealOpGreaterThan(NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				printf("OK\n");

				printf("SparseRealOpGreaterEqual... ");
				{
					if (SparseRealOpGreaterEqual(ps1, ps13RI, ps13CS, ps2, ps2RI, ps2CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpGreaterEqual(ps2, ps2RI, ps2CS, ps3, ps13RI, ps13CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpGreaterEqual(ps1, ps13RI, ps13CS, ps3, ps13RI, ps13CS, 1)
							!= 0) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}
					if (SparseRealOpGreaterEqual(pb, pbRI, pbCS, pa, paRI, paCS, 12) != 1) {
						printf("FAILED!!!\nat line %d\n", __LINE__);
						exit(EXIT_FAILURE);
					}

					SparseRealOpGreaterEqual(NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0);
				}
				printf("OK\n");
			}
		}

		printf("SparseRealOpVecMul... ");
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
			REAL pa[] = {1.0, 1.0, 1.0};
			INT paRI[] = {0, 0, 1};
			INT paCS[] = {0, 1, 1, 3};
			REAL pb[] = {1.0, 1.0, 1.0};
			REAL pr1[] = {2.0, 1.0};
			REAL *pr2 = (REAL *) malloc(2 * sizeof (REAL));

			SparseRealOpVecMul(pr2, pa, paRI, paCS, pb, 2, 3);
			if (real_failure(__LINE__, pr2, pr1, 2))
				exit(EXIT_FAILURE);

			SparseRealOpVecMul(NULL, NULL, NULL, emptyCS, NULL, 0, 0);
		}
		printf("OK\n");

		printf("SparseRealOpMatMul... ");
		{
			/*
			 * test cases contain:
			 *   empty columns in a, b
			 *   0 elements in r
			 *
			 *       /     \       /     \
			 *       | X o |       | o X |
			 *   a = |     |   b = |     |
			 *       | o o |       | o o |
			 *       \     /       \     /
			 */
			REAL pa[] = {1.0};
			INT paRI[] = {0};
			INT paCS[] = {0, 1, 1};
			REAL pb[] = {1.0};
			INT pbRI[] = {0};
			INT pbCS[] = {0, 0, 1};
			REAL pr[] = {1.0};
			REAL pr1[] = {1.0};
			INT pr1RI[] = {0};
			INT pr1CS[] = {0, 0, 1};
			REAL *pr2  = (REAL *) malloc(1 * sizeof (REAL));
			INT *pr2RI =  (INT *) malloc(1 * sizeof (INT));
			INT *pr2CS =  (INT *) malloc(3 * sizeof (INT));

			SparseRealOpMatMul(pr2, pr2RI, pr2CS, pa, paRI, paCS, pb, pbRI, pbCS, 2, 2);
			if (real_failure(__LINE__, pr2, pr1, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, pr2RI, pr1RI, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, pr2CS, pr1CS, 3))
				exit(EXIT_FAILURE);

			SparseRealOpMatMul(NULL, NULL, emptyCS,
												 NULL, NULL, emptyCS, NULL, NULL, emptyCS, 0, 0);
		}
		printf("OK\n");
	}
	printf("\n");

	printf("SPARSE_MATRIX\n");
	printf("-------------\n");
	{
		printf("SPARSE_MATRIX()... ");
		{
			SPARSE_MATRIX t1;
			SAFETYCLEAR(t1);
			if ((RowDimension(t1) != 0) || (ColDimension(t1) != 0)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (t1.colStarts[ColDimension(t1)] != 0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}

			SPARSE_MATRIX t2(2, 2);
			SAFETYCLEAR(t2);
			if ((RowDimension(t2) != 2) || (ColDimension(t2) != 2)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (t2.colStarts[ColDimension(t2)] != 0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}

			SPARSE_MATRIX t3(2, 2, 2);
			SAFETYCLEAR(t3);
			if ((RowDimension(t3) != 2) || (ColDimension(t3) != 2)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (t3.colStarts[ColDimension(t3)] != 0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("operator ()... ");
		{
			/*
			 * test cases contain:
			 *   return of empty elements
			 *     with later nonempty elements in column
			 *     without later nonempty elements in column
			 *   return of nonempty elements
			 */
			SPARSE_MATRIX t(3, 1, 1);
			t.theElements[0] = 1.0;
			t.rowIndices[0] = 1;
			t.colStarts[0] = 0;
			t.colStarts[1] = 1;

			if (t(1, 1) != 0.0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (t(2, 1) != 1.0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (t(3, 1) != 0.0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("IncAlloc... ");
		{
			/*
			 * test cases contain:
			 *   chopping arrays
			 *   extending arrays
			 */
			SPARSE_MATRIX t(2, 2, 2);
			SAFETYCLEAR(t);
			SetElement(t, 1, 1, 1.0);
			SetElement(t, 2, 2, 1.0);
			REAL theElements[] = {1.0, 1.0};
			INT  rowIndices[]  = {0, 1};

			IncAlloc(t, 1);
			if (Allocated(t) != 3) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (real_failure(__LINE__, t.theElements, theElements, 2))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t.rowIndices, rowIndices, 2))
				exit(EXIT_FAILURE);

			ErrorReport(0);
			ErrorHandler::LastErrorCode = 0;
			IncAlloc(t, -2);
			if (ErrorHandler::LastErrorCode != 1) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (Allocated(t) != 1) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (real_failure(__LINE__, t.theElements, theElements, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t.rowIndices, rowIndices, 1))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("SetElement... ");
		{
			/*
			 * test cases contain:
			 *   setting nonzero element to zero / nonzero value
			 *   setting zero element to nonzero value
			 *     with / without later nonzero element in column
			 *     with / without enough memory allocated
			 *   previous bugs
			 */
			SPARSE_MATRIX t(3, 1, 1);
			t.theElements[0] = 1.0;
			t.rowIndices[0] = 0;
			t.colStarts[0] = 0;
			t.colStarts[1] = 1;

			SetElement(t, 1, 1, 2.0);
			if (t(1, 1) != 2.0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			SetElement(t, 1, 1, 0.0);
			if (t(1, 1) != 0.0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			SetElement(t, 2, 1, 1.0);
			if (t(2, 1) != 1.0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			SetElement(t, 1, 1, 1.0);
			SetElement(t, 3, 1, 3.0);
			if (Allocated(t) != 3) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if ((t(1, 1) != 1.0) || (t(3, 1) != 3.0)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}

			printf("testing previous bugs... ");
			/* bug 1 in SparseMatrix.C:SetElement:
			 *   Setting a zero element to a nonzero value with enough memory
			 *   allocated the number of elements to move had wrong parenthesis.
			 *   was:
			 *     memmove([...], *nonZeros - offset * sizeof (REAL))
			 *   must be:
			 *     memmove([...], (*nonZeros - offset) * sizeof (REAL))
			 *   Applies analogously to SparseIntervalMatrix.C:SetElement.
			 * */
			SetElement(t, 3, 1, 0.0);
			SetElement(t, 3, 1, 3.0);
		}
		printf("OK\n");

		printf("operator =... ");
		{
#ifdef __DONTCOPY__
			{
				printf("(with __DONTCOPY__) ");
				/*
				 * test cases contain:
				 *   assignment temporary matrix to
				 *     nonempty / empty matrix
				 *
				 * note: MakeTemporary and MakePermanent not really testable
				 */
				SPARSE_MATRIX t1(1, 1, 1);
				SAFETYCLEAR(t1);
				SPARSE_MATRIX t2(2, 2, 2);

				MakeTemporary(t1);
				MakeTemporary(t2);
				SetElement (t1, 1, 1, 1.0);
				t2 = t1;
				if ((RowDimension(t2) != 1) || (ColDimension(t2) != 1)
						|| (RowDimension(t1) != 0) || (ColDimension(t1) != 0)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if ((Allocated(t2) != 1) || (Allocated(t1) != 0)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t1.colStarts[ColDimension(t1)] != 0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t2(1, 1) != 1.0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				t1 = t2;
				if ((RowDimension(t1) != 1) || (ColDimension(t1) != 1)
						|| (RowDimension(t2) != 0) || (ColDimension(t2) != 0)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if ((Allocated(t1) != 1) || (Allocated(t2) != 0)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t2.colStarts[ColDimension(t2)] != 0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t1(1, 1) != 1.0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
#endif
			{
				/*
				 * test cases contain:
				 *   assign non-temporary matrix to
				 *     empty matrix
				 *     nonempty matrix
				 *       with / without sufficient memory allocated
				 */
				SPARSE_MATRIX t1(2, 1, 1);
				SAFETYCLEAR(t1);
				SPARSE_MATRIX t2;
				SetElement(t1, 1, 1, 1.0);

				t2 = t1;
				if ((RowDimension(t2) != 2) || (ColDimension(t2) != 1)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (Allocated(t2) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t2(1, 1) != 1.0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				SetElement(t1, 1, 1, 2.0);
				t2 = t1;
				if ((RowDimension(t2) != 2) || (ColDimension(t2) != 1)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (Allocated(t2) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t2(1, 1) != 2.0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				SetElement(t1, 2, 1, 1.0);
				t2 = t1;
				if (Allocated(t2) != 2) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t2(2, 1) != 1.0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
		}
		printf("OK\n");

		{
			/*
			 * test cases contain:
			 *   check of
			 *     number of nonzero elements
			 *     size of allocated memory
			 *   in result matrix
			 */
			SPARSE_MATRIX t1(1, 1, 1);
			SAFETYCLEAR(t1);
			SPARSE_MATRIX t2(1, 1, 1);
			SAFETYCLEAR(t2);
			SetElement(t1, 1, 1, 1.0);
			SetElement(t2, 1, 1, 2.0);

			SPARSE_MATRIX t3(2, 2, 1);
			SAFETYCLEAR(t3);
			SPARSE_MATRIX t4(2, 2, 1);
			SAFETYCLEAR(t4);
			SetElement(t3, 1, 1, 1.0);
			SetElement(t4, 1, 1, 2.0);

			printf("operator +=... ");
			{
				t1 += t2;
				if (Allocated(t1) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t1.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				t3 += t4;
#ifdef __DONTCOPY__
				printf("(with __DONTCOPY__) ");
				if (Allocated(t3) != 2) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#else
				printf("(without __DONTCOPY__) ");
				if (Allocated(t3) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#endif
				if (t3.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
			printf("OK\n");

			printf("operator -=... ");
			{
				t1 -= t2;
				if (Allocated(t1) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (t1.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				t3 -= t4;
#ifdef __DONTCOPY__
				printf("(with __DONTCOPY__) ");
				if (Allocated(t3) != 2) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#else
				printf("(without __DONTCOPY__) ");
				if (Allocated(t3) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#endif
				if (t3.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
			printf("OK\n");
		}

		printf("Row... ");
		{
			/*
			 * test cases contain:
			 *   empty columns
			 *   nonempty columns
			 *     with / without an entry in the row
			 */
			SPARSE_MATRIX t(2, 3, 2);
			SAFETYCLEAR(t);
			SetElement(t, 1, 2, 1.0);
			SetElement(t, 2, 3, 1.0);
			VECTOR v1(3);
			Clear(v1);
			v1(3) = 1.0;

			VECTOR v2 = Row(t, 2);
			if (real_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("Col... ");
		{
			/*
			 * test cases contain:
			 *   empty / nonempty elements in column
			 */
			SPARSE_MATRIX t(2, 2, 1);
			SAFETYCLEAR(t);
			SetElement(t, 1, 2, 1.0);
			VECTOR v1(2);
			Clear(v1);
			v1(1) = 1.0;

			VECTOR v2 = Col(t, 2);
			if (real_failure(__LINE__, v2.theElements, v1.theElements, 2))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("SetRow... ");
		{
			/*
			 * test cases contain:
			 *   setting a row
			 *     with / without sufficient memory allocated
			 *     setting zero element to nonzero value
			 *       in the middle / at the end of a column
			 *     setting nonzero element
			 *       to zero
			 *       to nonzero value
			 *
			 *       /       \
			 *       | o o o |
			 *   t = |       |
			 *       | X X o |
			 *       \       /
			 */
			SPARSE_MATRIX t(2, 3, 5);
			SAFETYCLEAR(t);
			SetElement(t, 2, 1, 1.0);
			SetElement(t, 2, 2, 1.0);
			VECTOR v1(3);
			Clear(v1);
			v1(1) = v1(3) = 1.0;
			VECTOR v2(3);

			SetRow(t, 1, v1);
			v2 = Row(t, 1);
			if (real_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
			SetRow(t, 2, v1);
			v2 = Row(t, 2);
			if (real_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("SetCol... ");
		{
			/*
			 * test cases contain:
			 *   setting a column
			 *     with zero elements
			 *     with / without sufficient memory allocated
			 *     increasing / decreasing element count in column
			 *
			 *       /      \
			 *       | X  o |
			 *       |      |
			 *   t = | X  o |
			 *       |      |
			 *       | o  o |
			 *       \      /
			 */
			SPARSE_MATRIX t(3, 2, 2);
			SAFETYCLEAR(t);
			SetElement(t, 1, 1, 1.0);
			SetElement(t, 2, 1, 1.0);
			VECTOR v1(3);
			Initialize(v1, 1.0);
			v1(3) = 0.0;
			VECTOR v2(3);

			SetCol(t, 2, v1);
			v2 = Col(t, 2);
			if (real_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
			v1(2) = 0.0;
			SetCol(t, 1, v1);
			v2 = Col(t, 1);
			if (real_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("Abs... ");
		{
			SPARSE_MATRIX t1(2, 3, 2);
			SAFETYCLEAR(t1);
			SetElement(t1, 1, 2, -1.0);
			SetElement(t1, 2, 3, 1.0);
			SPARSE_MATRIX t2(t1);
			MakePermanent(t2);
			SetElement(t2, 1, 2, 1.0);

			if (!((Abs(t1) >= t2) && (Abs(t1) <= t2))) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("operator <<... ");
		{
			/*
			 * test cases contain:
			 *   empty column
			 */
			SPARSE_MATRIX t(2, 2, 2);
			SAFETYCLEAR(t);
			SetElement(t, 1, 1, 1.0);
			SetElement(t, 2, 1, 2.0);
			ostringstream os;

			os << t;
			if (os.str() != "((1, 1, 1) ; (2, 1, 2))") {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				printf("output stream reads:\n");
				printf("%s\n", os.str().c_str());
				printf("but should read:\n");
				printf("((1, 1, 1) ; (2, 1, 2))\n");
				exit(EXIT_FAILURE);
			}

			Clear(t);
			os.str("");
			os << t;
			if (os.str() != "()") {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				printf("output stream reads:\n");
				printf("%s\n", os.str().c_str());
				printf("but should read:\n");
				printf("()\n");
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("operator >>... ");
		{
			/*
			 * test cases contain:
			 *   empty columns at the beginning and at the end of the matrix
			 */
			SPARSE_MATRIX t1(2, 3, 1);
			SAFETYCLEAR(t1);
			SetElement(t1, 1, 2, 1.0);
			SPARSE_MATRIX t2(2, 3, 1);
			SAFETYCLEAR(t2);
			istringstream is("1 2 1.0");

			is >> t2;
			if (real_failure(__LINE__, t2.theElements, t1.theElements, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t2.rowIndices, t1.rowIndices, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t2.colStarts, t1.colStarts, 4))
				exit(EXIT_FAILURE);

			Clear(t1);
			Resize(t2, 2, 3, 0);
			is.str("");

			is >> t2;
		}
		printf("OK\n");

		printf("calling trivial operations:\n");
		printf("  MakeTemporary\n");
		printf("  MakePermanent\n");
		printf("  operator *=\n");
		printf("  operator /=\n");
		printf("  RowDimension\n");
		printf("  ColDimension\n");
		printf("  Allocated\n");
		printf("  Resize\n");
		printf("  Clear\n");
		printf("  operator +(CONST SPARSE_MATRIX &)\n");
		printf("  operator -(CONST SPARSE_MATRIX &)\n");
		printf("  operator +(CONST SPARSE_MATRIX &, CONST SPARSE_MATRIX &)\n");
		printf("  operator -(CONST SPARSE_MATRIX &, CONST SPARSE_MATRIX &)\n");
		printf("  operator *(REAL, CONST SPARSE_MATRIX &)\n");
		printf("  operator /(CONST SPARSE_MATRIX &, CONST REAL)\n");
		printf("  operator *(CONST SPARSE_MATRIX &, CONST VECTOR &)\n");
		printf("  operator *(CONST SPARSE_MATRIX &, CONST SPARSE_MATRIX &)\n");
		printf("  operator  <(CONST SPARSE_MATRIX &, CONST SPARSE_MATRIX &)\n");
		printf("  operator <=(CONST SPARSE_MATRIX &, CONST SPARSE_MATRIX &)\n");
		printf("  operator  >(CONST SPARSE_MATRIX &, CONST SPARSE_MATRIX &)\n");
		printf("  operator >=(CONST SPARSE_MATRIX &, CONST SPARSE_MATRIX &)\n");
		printf("... ");
		{
			SPARSE_MATRIX t(2, 2, 2);
			SAFETYCLEAR(t);
			VECTOR v(2);
			Clear(v);

			MakeTemporary(t);
			MakePermanent(t);
			t *= 2.0;
			t /= 2.0;
			RowDimension(t);
			ColDimension(t);
			Allocated(t);
			Resize(t, 2, 2);
			Clear(t);
			+t;
			-t;
			t + t;
			t - t;
			2 * t;
			t / 2;
			t * v;
			t * t;
			t < t;
			t <= t;
			t > t;
			t >= t;

			printf("testing previous bugs... ");
			{
				/* bug 2 in SparseMatrix.C:Resize(SPARSE_MATRIX &, INT, INT, INT):
				 *   When Resize is called with the same row and column dimensions no
				 *   check was made for the number of allocated elements.
				 *   was:
				 *     if ((r * c) && (x.nRows == r) && (x.nCols == c)) return;
				 *   must be:
				 *     if ((r * c) 
				 *         && (x.nRows == r) && (x.nCols == c) && (x.nAlloc == a))
				 *       return;
				 *   Applies analogously to SparseIntervalMatrix.C:Resize([...]).
				 * */
				SPARSE_MATRIX t(2, 2, 2);
				SAFETYCLEAR(t);

				Resize(t, 2, 2, 1);
				if (Allocated(t) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				/* bug 3 in 
				 * SparseMatrix.C:operator *(CONST REAL, CONST SPARSE_MATRIX &)
				 * and
				 * SparseMaitrx.C:operator /(CONST SPARSE_MATRIX &, CONST REAL):
				 *   The result matrix t is initialized from the argument matrix. If
				 *   this is a temporary one, it is emptied resulting in RealOpMul or
				 *   RealOpDiv doing nothing and t just being the argument matrix times
				 *   or divided by 1. We could simply declare the result and copy the
				 *   rowIndices and colStarts arrays, which the dense routines don't
				 *   take care of. But using RealOpMulWith and RealOpDivBy instead saves
				 *   us these copy operations.
				 *   was:
				 *     operator *: RealOpMul([...]);
				 *     operator /: RealOpDiv([...]);
				 *   must be:
				 *     operator *: RealOpMulWith([...]);
				 *     operator /: RealOpDivBy([...]);
				 * */
				SPARSE_MATRIX i;
				i = SpId(2);
				SPARSE_MATRIX t1;
				t1 = 0.5 * i;
				SPARSE_MATRIX t2 = i;
				t2 = 0.5 * t2;

				if (real_failure(__LINE__, t2.theElements, t1.theElements, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, t2.rowIndices, t1.rowIndices, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, t2.colStarts, t1.colStarts, 3))
					exit(EXIT_FAILURE);

				t1 = i / 2;
				SPARSE_MATRIX t3 = i;
				t3 = t3 / 2;

				if (real_failure(__LINE__, t3.theElements, t1.theElements, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, t3.rowIndices, t1.rowIndices, 2))
					exit(EXIT_FAILURE);
				if (int_failure(__LINE__, t3.colStarts, t1.colStarts, 3))
					exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");
	}
	printf("\n");

	printf("SPARSE_INTERVAL_MATRIX\n");
	printf("----------------------\n");
	{
		printf("SPARSE_INTERVAL_MATRIX()... ");
		{
			SPARSE_INTERVAL_MATRIX T1;
			SAFETYCLEAR(T1);
			if ((RowDimension(T1) != 0) || (ColDimension(T1) != 0)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (T1.colStarts[ColDimension(T1)] != 0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}

			SPARSE_INTERVAL_MATRIX T2(2, 2);
			SAFETYCLEAR(T2);
			if ((RowDimension(T2) != 2) || (ColDimension(T2) != 2)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (T2.colStarts[ColDimension(T2)] != 0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}

			SPARSE_INTERVAL_MATRIX T3(2, 2, 2);
			SAFETYCLEAR(T3);
			if ((RowDimension(T3) != 2) || (ColDimension(T3) != 2)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (T3.colStarts[ColDimension(T3)] != 0) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("~SPARSE_INTERVAL_MATRIX()... ");
		{
			SPARSE_INTERVAL_MATRIX *T = new SPARSE_INTERVAL_MATRIX;
			delete T;
		}
		printf("OK\n");

		printf("operator ()... ");
		{
			/*
			 * test cases contain:
			 *   return of empty elements
			 *     with later nonempty elements in column
			 *     without later nonempty elements in column
			 *   return of nonempty elements
			 */
			SPARSE_INTERVAL_MATRIX T(3, 1, 1);
			INTERVAL zero(0.0, 0.0), one(1.0, 1.0);
			T.theElements[0] = *Bias(one);
			T.rowIndices[0] = 1;
			T.colStarts[0] = 0;
			T.colStarts[1] = 1;

			if (T(1, 1) != zero) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (T(2, 1) != one) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (T(3, 1) != zero) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("IncAlloc... ");
		{
			/*
			 * test cases contain:
			 *   chopping arrays
			 *   extending arrays
			 */
			INTERVAL one(1.0, 1.0);
			SPARSE_INTERVAL_MATRIX T(2, 2, 2);
			SAFETYCLEAR(T);
			SetElement(T, 1, 1, one);
			SetElement(T, 2, 2, one);
			BIASINTERVAL theElements[] = {{1.0, 1.0}, {1.0, 1.0}};
			INT rowIndices[]  = {0, 1};

			IncAlloc(T, 1);
			if (Allocated(T) != 3) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (interval_failure(__LINE__, T.theElements, theElements, 2))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, T.rowIndices, rowIndices, 2))
				exit(EXIT_FAILURE);

			ErrorReport(0);
			ErrorHandler::LastErrorCode = 0;
			IncAlloc(T, -2);
			if (ErrorHandler::LastErrorCode != 1) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (Allocated(T) != 1) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if (interval_failure(__LINE__, T.theElements, theElements, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, T.rowIndices, rowIndices, 1))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("SetElement... ");
		{
			/*
			 * test cases contain:
			 *   setting nonzero element to zero / nonzero value
			 *   setting zero element to nonzero value
			 *     with / without later nonzero element in column
			 *     with / without enough memory allocated
			 *   previous bugs
			 */
			SPARSE_INTERVAL_MATRIX T(3, 1, 1);
			INTERVAL zero(0.0, 0.0), one(1.0, 1.0), two(2.0, 2.0), three(3.0, 3.0);
			T.theElements[0] = *Bias(one);
			T.rowIndices[0] = 0;
			T.colStarts[0] = 0;
			T.colStarts[1] = 1;

			SetElement(T, 1, 1, two);
			if (T(1, 1) != two) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			SetElement(T, 1, 1, zero);
			if (T(1, 1) != zero) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			SetElement(T, 2, 1, one);
			if (T(2, 1) != one) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			SetElement(T, 1, 1, one);
			SetElement(T, 3, 1, three);
			if (Allocated(T) != 3) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
			if ((T(1, 1) != one) || (T(3, 1) != three)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}

			printf("testing previous bugs... ");
			/* bug 1 in SparseMatrix.C:SetElement:
			 *   Setting a zero element to a nonzero value with enough memory
			 *   allocated the number of elements to move had wrong parenthesis.
			 *   was:
			 *     memmove([...], *nonZeros - offset * sizeof (REAL))
			 *   must be:
			 *     memmove([...], (*nonZeros - offset) * sizeof (REAL))
			 *   Applies analogously to SparseIntervalMatrix.C:SetElement.
			 */
			SetElement(T, 3, 1, zero);
			SetElement(T, 3, 1, three);
		}
		printf("OK\n");

		printf("operator =(CONST SPARSE_INTERVAL_MATRIX &)... ");
		{
			INTERVAL one(1.0, 1.0);
#ifdef __DONTCOPY__
			{
				printf("(with __DONTCOPY__) ");
				/*
				 * test cases contain:
				 *   assignment temporary matrix to
				 *     nonempty / empty matrix
				 *
				 * note: MakeTemporary and MakePermanent not really testable
				 */
				SPARSE_INTERVAL_MATRIX T1(1, 1, 1);
				SAFETYCLEAR(T1);
				SPARSE_INTERVAL_MATRIX T2(2, 2, 2);

				MakeTemporary(T1);
				MakeTemporary(T2);
				SetElement (T1, 1, 1, one);
				T2 = T1;
				if ((RowDimension(T2) != 1) || (ColDimension(T2) != 1)
						|| (RowDimension(T1) != 0) || (ColDimension(T1) != 0)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if ((Allocated(T2) != 1) || (Allocated(T1) != 0)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T1.colStarts[ColDimension(T1)] != 0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T2(1, 1) != one) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				T1 = T2;
				if ((RowDimension(T1) != 1) || (ColDimension(T1) != 1)
						|| (RowDimension(T2) != 0) || (ColDimension(T2) != 0)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if ((Allocated(T1) != 1) || (Allocated(T2) != 0)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T2.colStarts[ColDimension(T2)] != 0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T1(1, 1) != one) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
#endif
			{
				/*
				 * test cases contain:
				 *   assign non-temporary matrix to
				 *     empty matrix
				 *     nonempty matrix
				 *       with / without sufficient memory allocated
				 */
				SPARSE_INTERVAL_MATRIX T1(2, 1, 1);
				SAFETYCLEAR(T1);
				SPARSE_INTERVAL_MATRIX T2;
				SetElement(T1, 1, 1, one);

				T2 = T1;
				if ((RowDimension(T2) != 2) || (ColDimension(T2) != 1)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (Allocated(T2) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T2(1, 1) != one) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				SetElement(T1, 1, 1, 2.0);
				T2 = T1;
				if ((RowDimension(T2) != 2) || (ColDimension(T2) != 1)) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (Allocated(T2) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T2(1, 1) != 2.0) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				SetElement(T1, 2, 1, one);
				T2 = T1;
				if (Allocated(T2) != 2) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T2(2, 1) != one) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
		}
		printf("OK\n");

		printf("operator =(CONST SPARSE_MATRIX &)... ");
		{
			/*
			 * test cases contain:
			 *   assign to empty matrix
			 *   assign to nonempty matrix
			 *     with / without sufficient memory allocated
			 */
			SPARSE_MATRIX t(2, 2, 2);
			SAFETYCLEAR(t);
			SetElement(t, 1, 1, 1.0);
			SPARSE_INTERVAL_MATRIX T;

			T = t;
			T = t;
			SetElement(t, 2, 2, 1.0);
			T = t;
		}
		printf("OK\n");

		{
			/*
			 * test cases contain:
			 *   check of
			 *     number of nonzero elements
			 *     size of allocated memory
			 *       with / without trivial upper bound active
			 *   in result matrix
			 */
			INTERVAL one(1.0, 1.0), two(2.0, 2.0);

			SPARSE_INTERVAL_MATRIX T1(1, 1, 1);
			SAFETYCLEAR(T1);
			SPARSE_MATRIX t2(1, 1, 1);
			SAFETYCLEAR(t2);
			SPARSE_INTERVAL_MATRIX T2(1, 1, 1);
			SAFETYCLEAR(T2);
			SetElement(T1, 1, 1, one);
			SetElement(t2, 1, 1, 2.0);
			SetElement(T2, 1, 1, two);

			SPARSE_INTERVAL_MATRIX T3(2, 2, 1);
			SAFETYCLEAR(T3);
			SPARSE_MATRIX t4(2, 2, 1);
			SAFETYCLEAR(t4);
			SPARSE_INTERVAL_MATRIX T4(2, 2, 1);
			SAFETYCLEAR(T4);
			SetElement(T3, 1, 1, one);
			SetElement(t4, 1, 1, 2.0);
			SetElement(T4, 1, 1, two);

			printf("operator +=... ");
			{
				T1 += t2;
				if (Allocated(T1) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T1.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				T1 += T2;
				if (Allocated(T1) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T1.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				T3 += t4;
#ifdef __DONTCOPY__
				printf("(with __DONTCOPY__) ");
				if (Allocated(T3) != 2) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#else
				printf("(without __DONTCOPY__) ");
				if (Allocated(T3) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#endif
				if (T3.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				T3 += T4;
#ifdef __DONTCOPY__
				printf("(with __DONTCOPY__) ");
				if (Allocated(T3) != 2) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#else
				printf("(without __DONTCOPY__) ");
				if (Allocated(T3) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#endif
				if (T3.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
			printf("OK\n");

			printf("operator -=... ");
			{
				T1 -= t2;
				if (Allocated(T1) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T1.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				T1 -= T2;
				if (Allocated(T1) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
				if (T1.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				T3 -= t4;
#ifdef __DONTCOPY__
				printf("(with __DONTCOPY__) ");
				if (Allocated(T3) != 2) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#else
				printf("(without __DONTCOPY__) ");
				if (Allocated(T3) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#endif
				if (T3.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}

				T3 -= T4;
#ifdef __DONTCOPY__
				printf("(with __DONTCOPY__) ");
				if (Allocated(T3) != 2) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#else
				printf("(without __DONTCOPY__) ");
				if (Allocated(T3) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
#endif
				if (T3.colStarts[1] != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
			printf("OK\n");
		}

		printf("Row... ");
		{
			/*
			 * test cases contain:
			 *   empty columns
			 *   nonempty columns
			 *     with / without an entry in the row
			 */
			INTERVAL one(1.0, 1.0);
			SPARSE_INTERVAL_MATRIX T(2, 3, 2);
			SAFETYCLEAR(T);
			SetElement(T, 1, 2, one);
			SetElement(T, 2, 3, one);
			INTERVAL_VECTOR v1(3);
			Clear(v1);
			v1(3) = one;

			INTERVAL_VECTOR v2 = Row(T, 2);
			if (interval_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("Col... ");
		{
			/*
			 * test cases contain:
			 *   empty / nonempty elements in column
			 */
			INTERVAL one(1.0, 1.0);
			SPARSE_INTERVAL_MATRIX T(2, 2, 1);
			SAFETYCLEAR(T);
			SetElement(T, 1, 2, one);
			INTERVAL_VECTOR v1(2);
			Clear(v1);
			v1(1) = one;

			INTERVAL_VECTOR v2 = Col(T, 2);
			if (interval_failure(__LINE__, v2.theElements, v1.theElements, 2))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("SetRow... ");
		{
			/*
			 * test cases contain:
			 *   setting a row
			 *     with / without sufficient memory allocated
			 *     setting zero element to nonzero value
			 *       in the middle / at the end of a column
			 *     setting nonzero element
			 *       to zero
			 *       to nonzero value
			 *
			 *       /       \
			 *       | o o o |
			 *   T = |       |
			 *       | X X o |
			 *       \       /
			 */
			INTERVAL one(1.0, 1.0);
			SPARSE_INTERVAL_MATRIX T(2, 3, 5);
			SAFETYCLEAR(T);
			SetElement(T, 2, 1, one);
			SetElement(T, 2, 2, one);
			INTERVAL_VECTOR v1(3);
			Clear(v1);
			v1(1) = v1(3) = one;
			INTERVAL_VECTOR v2(3);

			SetRow(T, 1, v1);
			v2 = Row(T, 1);
			if (interval_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
			SetRow(T, 2, v1);
			v2 = Row(T, 2);
			if (interval_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("SetCol... ");
		{
			/*
			 * test cases contain:
			 *   setting a column
			 *     with zero elements
			 *     with / without sufficient memory allocated
			 *     increasing / decreasing element count in column
			 *
			 *       /      \
			 *       | X  o |
			 *       |      |
			 *   T = | X  o |
			 *       |      |
			 *       | o  o |
			 *       \      /
			 */
			INTERVAL one(1.0, 1.0);
			SPARSE_INTERVAL_MATRIX T(3, 2, 2);
			SAFETYCLEAR(T);
			SetElement(T, 1, 1, one);
			SetElement(T, 2, 1, one);
			INTERVAL_VECTOR v1(3);
			Initialize(v1, one);
			v1(3) = 0.0;
			INTERVAL_VECTOR v2(3);

			SetCol(T, 2, v1);
			v2 = Col(T, 2);
			if (interval_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
			v1(2) = 0.0;
			SetCol(T, 1, v1);
			v2 = Col(T, 1);
			if (interval_failure(__LINE__, v2.theElements, v1.theElements, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("Intersection... ");
		{
			/*
			 * test cases contain:
			 *   assign intersection to
			 *     empty matrix
			 *     nonempty matrix w/o sufficient mempory allocated
			 */
			INTERVAL one(1.0, 1.0);
			SPARSE_INTERVAL_MATRIX T;
			SAFETYCLEAR(T);
			SPARSE_INTERVAL_MATRIX a(2, 2, 2);
			SAFETYCLEAR(a);
			SetElement(a, 1, 1, one);

			Intersection(T, a, a);
			SetElement(a, 2, 2, one);
			Intersection(T, a, a);
		}
		printf("OK\n");

		printf("operator <<... ");
		{
			/*
			 * test cases contain:
			 *   empty column
			 */
			INTERVAL one(1.0, 1.0), two(2.0, 2.0);
			SPARSE_INTERVAL_MATRIX t(2, 2, 2);
			SAFETYCLEAR(t);
			SetElement(t, 1, 1, one);
			SetElement(t, 2, 1, two);
			ostringstream os;

			os << t;
			if (os.str() != "((1, 1, [1,1]) ; (2, 1, [2,2]))") {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				printf("output stream reads:\n");
				printf("%s\n", os.str().c_str());
				printf("but should read:\n");
				printf("((1, 1, [1,1]) ; (2, 1, [2,2]))\n");
				exit(EXIT_FAILURE);
			}

			Clear(t);
			os.str("");
			os << t;
			if (os.str() != "()") {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				printf("output stream reads:\n");
				printf("%s\n", os.str().c_str());
				printf("but should read:\n");
				printf("()\n");
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("operator >>... ");
		{
			/*
			 * test cases contain:
			 *   empty columns at the beginning and at the end of the matrix
			 */
			INTERVAL one(1.0, 1.0);
			SPARSE_INTERVAL_MATRIX t1(2, 3, 1);
			SAFETYCLEAR(t1);
			SetElement(t1, 1, 2, one);
			SPARSE_INTERVAL_MATRIX t2(2, 3, 1);
			SAFETYCLEAR(t2);
			istringstream is("1 2 1.0 1.0");

			is >> t2;
			if (interval_failure(__LINE__, t2.theElements, t1.theElements, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t2.rowIndices, t1.rowIndices, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t2.colStarts, t1.colStarts, 4))
				exit(EXIT_FAILURE);

			Clear(t1);
			Resize(t2, 2, 3, 0);
			is.str("");

			is >> t2;
		}
		printf("OK\n");

		printf("\n");
		printf("calling trivial operations:\n");
		printf("  MakeTemporary\n");
		printf("  MakePermanent\n");
		printf("  operator *=\n");
		printf("  operator /=\n");
		printf("  RowDimension\n");
		printf("  ColDimension\n");
		printf("  Allocated\n");
		printf("  Resize\n");
		printf("  Clear\n");
		printf("  Inf\n");
		printf("  Sup\n");
		printf("  Mid\n");
		printf("  Diam\n");
		printf("  Abs\n");
		printf("  operator <=\n");
		printf("  operator <\n");
		printf("  operator ==\n");
		printf("  operator !=\n");
		printf("  Hull\n");
		printf("  SymHull\n");
		printf("  operator +(CONST SPARSE_INTERVAL_MATRIX &)\n");
		printf("  operator -(CONST SPARSE_INTERVAL_MATRIX &)\n");
		printf("  operator +\n");
		printf("  AddBounds\n");
		printf("  operator -\n");
		printf("  SubBounds\n");
		printf("  operator * [scalar x matrix]\n");
		printf("  MulBounds  [scalar x matrix]\n");
		printf("  operator /\n");
		printf("  DivBounds\n");
		printf("  operator * [matrix x vector]\n");
		printf("  MulBounds  [matrix x vector]\n");
		printf("  operator * [matrix x matrix]\n");
		printf("  MulBounds  [matrix x matrix]\n");
		printf("... ");
		{
			SPARSE_INTERVAL_MATRIX T(2, 2, 2);
			SAFETYCLEAR(T);
			SPARSE_MATRIX t(2, 2, 2);
			SAFETYCLEAR(t);
			INTERVAL_MATRIX D(2, 2);
			Clear(D);
			INTERVAL two(2.0, 2.0);
			VECTOR v(2);
			Clear(v);
			INTERVAL_VECTOR V(2);
			Clear(V);

			MakeTemporary(T);
			MakePermanent(T);
			T *= 2;
			T /= 2;
			RowDimension(T);
			ColDimension(T);
			Allocated(T);
			Resize(T, 2, 2);
			Clear(T);
			Inf(T);
			Sup(T);
			Mid(T);
			Diam(T);
			Abs(T);
			t <= T; T <= T;
			t <  D; T <  D;
			T == T;
			T != T;
			Hull(t); Hull(t, t); Hull(t, T); Hull(T, t); Hull(T, T);
			SymHull(t);
			+T;
			-T;
			t + T; T + t; T + T;
			AddBounds(t, t);
			t - T; T - t; T - T;
			SubBounds(t, t);
			2 * T; two * t; two * T;
			MulBounds(t, 2);
			T / 2; t / two; T / two;
			DivBounds(t, 2);
			t * V; T * v; T * V;
			MulBounds(t, v);

			printf("testing previous bugs... ");
			{
				/* bug 2 in SparseMatrix.C:Resize(SPARSE_MATRIX &, INT, INT, INT):
				 *   When Resize is called with the same row and column dimensions no
				 *   check was made for the number of allocated elements.
				 *   was:
				 *     if ((r * c) && (x.nRows == r) && (x.nCols == c)) return;
				 *   must be:
				 *     if ((r * c) 
				 *         && (x.nRows == r) && (x.nCols == c) && (x.nAlloc == a))
				 *       return;
				 *   Applies analogously to SparseIntervalMatrix.C:Resize([...]).
				 * */
				SPARSE_MATRIX t(2, 2, 2);

				Resize(t, 2, 2, 1);
				if (Allocated(t) != 1) {
					printf("FAILED!!!\nat line %d\n", __LINE__);
					exit(EXIT_FAILURE);
				}
			}
		}
		printf("OK\n");
	}
	printf("\n");

	printf("Sparse utility routines (SparseUtilities.h)\n");
	printf("-------------------------------------------\n");
	{
		printf("SpId... ");
		{
			SPARSE_MATRIX t2(3, 3, 3);
			SAFETYCLEAR(t2);
			SetElement(t2, 1, 1, 1.0);
			SetElement(t2, 2, 2, 1.0);
			SetElement(t2, 3, 3, 1.0);
			SPARSE_MATRIX t1 = SpId(3);

			if (real_failure(__LINE__, t2.theElements, t1.theElements, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t2.rowIndices, t1.rowIndices, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t2.colStarts, t1.colStarts, 4))
				exit(EXIT_FAILURE);

		}
		printf("OK\n");

		printf("MaxColSum... ");
		{
			SPARSE_MATRIX t = SpId(2);
			SAFETYCLEAR(t);
			SetElement(t, 1, 2, -1.0);

			if (MaxColSum(t) != 2) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("AAT... ");
		{
			SPARSE_MATRIX t1(2, 1, 2);
			SAFETYCLEAR(t1);
			SetElement(t1, 1, 1, 1.0);
			SetElement(t1, 2, 1, 2.0);
			SPARSE_MATRIX t2(2, 2, 4);
			SAFETYCLEAR(t2);
			SetElement(t2, 1, 1, 1.0);
			SetElement(t2, 1, 2, 2.0);
			SetElement(t2, 2, 1, 2.0);
			SetElement(t2, 2, 2, 4.0);

			SPARSE_MATRIX t3 = AAT(t1, 2);

			if (real_failure(__LINE__, t3.theElements, t2.theElements, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t3.rowIndices, t2.rowIndices, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t3.colStarts, t2.colStarts, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("ATA... ");
		{
			SPARSE_MATRIX t1(1, 2, 2);
			SAFETYCLEAR(t1);
			SetElement(t1, 1, 1, 1.0);
			SetElement(t1, 1, 2, 2.0);
			SPARSE_MATRIX t2(2, 2, 4);
			SAFETYCLEAR(t2);
			SetElement(t2, 1, 1, 1.0);
			SetElement(t2, 1, 2, 2.0);
			SetElement(t2, 2, 1, 2.0);
			SetElement(t2, 2, 2, 4.0);

			SPARSE_MATRIX t3 = ATA(t1, 2);

			if (real_failure(__LINE__, t3.theElements, t2.theElements, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t3.rowIndices, t2.rowIndices, 1))
				exit(EXIT_FAILURE);
			if (int_failure(__LINE__, t3.colStarts, t2.colStarts, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");
	}
	printf("\n");

	printf("Sparse verified solver routines(SLSS.h)\n");
	printf("---------------------------------------\n");
	{
		printf("Cholesky/SolveCholesky... ");
		{
			/* These routine is based on the LDL package by Timothy Davis, so let's
			 * hope he knows what he's doing and just check rudimentary =).
			 * */
			SPARSE_MATRIX t = SpId(3);
			SAFETYCLEAR(t);
			SetElement(t, 1, 2, 1.0);
			SetElement(t, 2, 1, 1.0);
			SPARSE_MATRIX c(3, 3);

			if (Cholesky(t, NULL, c) != 1) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}

			SetElement(t, 2, 2, 2.0);
			if (Cholesky(t, NULL, c) != 3) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}

			VECTOR e(3);
			Initialize(e, 1.0);
			VECTOR b = t * e;

			if (real_failure(__LINE__, 
											 (SolveLLT(c, b)).theElements, e.theElements, 3))
				exit(EXIT_FAILURE);
		}
		printf("OK\n");

		printf("MinSing... ");
		{
#warning MinSing test missing
			/*
			 * Test cases include:
			 *   singular matrix
			 *   indefinite matrix
			 *   s.p.d. matrix
			 */
			SPARSE_MATRIX t = SpId(3);
			SetElement(t, 1, 2, 1.0);
			SetElement(t, 2, 1, 1.0);
			INT p[3];
			SPARSE_MATRIX c(3, 3);

			REAL sv = MinSing(t, p, c);

			if ((converged != 1) || (fabs(ev - 1.0) > BiasEpsilon)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");

		printf("SLSS... ");
		{
			SPARSE_MATRIX t = SpId(3);
			SAFETYCLEAR(t);
			SetElement(t, 1, 2, 1.0);
			SetElement(t, 2, 1, 1.0);
			SetElement(t, 2, 2, 2.0);
			VECTOR e(3);
			Initialize(e, 1.0);

			INT info;
			INTERVAL_VECTOR X = SLSS(Hull(t), Hull(t * e), info);

			if ((info != 1) || !(e <= X)) {
				printf("FAILED!!!\nat line %d\n", __LINE__);
				exit(EXIT_FAILURE);
			}
		}
		printf("OK\n");
	}
	printf("\n");
	printf("-= Passed\n");

	exit(EXIT_SUCCESS);
}

