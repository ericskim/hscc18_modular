/************************************************************************
 *
 * Basic Interval Arithmetic Subroutines x86-64-Rounding (Linux ELF)
 * -----------------------------------------------------------------
 *
 * Copyright (C) 2009 Christian Keil
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
 * $Id$
 *
 ************************************************************************/

#ifndef __BIASRND__
#define __BIASRND__

#include <BIAS/BiasType.h>

#undef __BIASINLINEROUNDING__

/*
 * If the following macro is defined,
 * division by zero aborts with an error message;
 * otherwise, bounds containing +/-oo
 * or NaN are computed.
 */
#define __BIASRAISEDIVIDEBYZERO__

/*
 * If the following macro is defined,
 * the rounding mode is set to nearest after
 * each interval operation.
 *
 * 2008-10-27
 * The libm math library from glibc is in error when standard functions are computed in non-standard
 * rounding mode (Sources Bugzilla bug #3976) on x86-64.
 * For now Keep __BIASSETROUNDTONEAREST__ defined.
 */
#define  __BIASSETROUNDTONEAREST__

/*
 * The following macro defines the number of invalid bits
 * obtained by the standard C library functions (e.g. sin(), cos(), ...)
 */
#define __BIASSTDFUNCINVALIDBITS__      0

#ifdef __cplusplus
extern "C" {
#endif

void BiasRoundInit (void);
void BiasRoundUp   (void);
void BiasRoundDown (void);
void BiasRoundNear (void);

#ifdef __cplusplus
}
#endif

#endif /* __BIASRND__ */
