/************************************************************************
 *
 * Implementation of standard functions (REAL and INTERVAL)
 * --------------------------------------------------------
 *
 * Copyright (C) 1993, 1997 Olaf Knueppel
 *               2005 Christian Keil
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
 * $Id: Functions.C 478 2006-08-09 13:13:30Z keil $
 *
 ************************************************************************/

static char rcs_id[] = "$Id: Functions.C 478 2006-08-09 13:13:30Z keil $";

#include <Functions.h>
#include <BIAS/BiasF.h>

REAL Constant::Pi         = 0.0;
REAL Constant::TwoPi      = 0.0;
REAL Constant::PiHalf     = 0.0;
REAL Constant::PiQuarter  = 0.0;
REAL Constant::e          = 0.0;
REAL Constant::Sqrt2      = 0.0;
REAL Constant::InvSqrt2   = 0.0;
REAL Constant::Ln10       = 0.0;

Constant UsefulConstants; // constructor is called automatically

Constant::Constant ()
{
  BiasFuncInit ();  

  Constant::Pi         = BiasPi;
  Constant::TwoPi      = BiasTwoPi;
  Constant::PiHalf     = BiasPiHalf;
  Constant::PiQuarter  = BiasPiQuarter;
  Constant::e          = BiasE;
  Constant::Sqrt2      = BiasSqrt2;
  Constant::InvSqrt2   = BiasInvSqrt2;
  Constant::Ln10       = BiasLn10;
}

REAL ArSinh (REAL x)
{
  return log (x + sqrt (x * x + 1.0));
}

REAL ArCosh (REAL x)
{
  return log (x + sqrt (x * x - 1.0));
}

REAL ArTanh (REAL x)
{
  return 0.5 * log ((1.0 + x) / (1.0 - x));
}

REAL ArCoth (REAL x)
{
  return 0.5 * log ((x + 1.0) / (x - 1.0));
}

REAL Power (REAL x, INT n)
// Calculates x^n, all cases are considered
{
  INT i, absn;
  REAL y = 1.0;

  absn = (n < 0) ? (-n) : n;
  for (i = 1; i <= absn; i++) y *= x;
  if (n < 0) return 1.0 / y;
  else return y;
}

INTERVAL Sin (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasSin (Bias(r), Bias(x));
  return r;
}

INTERVAL Cos (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasCos (Bias(r), Bias(x));
  return r;
}

INTERVAL Tan (CONST INTERVAL & x)
{
  INTERVAL r;
  BiasTan (Bias(r), Bias(x));
  return r;
}

INTERVAL Cot (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasCot (Bias(r), Bias(x));
  return r;
}

INTERVAL ArcSin (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasArcSin (Bias(r), Bias(x));
  return r;
}

INTERVAL ArcCos (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasArcCos (Bias(r), Bias(x));
  return r;
}

INTERVAL ArcTan (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasArcTan (Bias(r), Bias(x));
  return r;
}

INTERVAL ArcCot (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasArcCot (Bias(r), Bias(x));
  return r;
}

INTERVAL Sinh (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasSinh (Bias(r), Bias(x));
  return r;
}

INTERVAL Cosh (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasCosh (Bias(r), Bias(x));
  return r;
}

INTERVAL Tanh (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasTanh (Bias(r), Bias(x));
  return r;
}

INTERVAL Coth (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasCoth (Bias(r), Bias(x));
  return r;
}

INTERVAL ArSinh (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasArSinh (Bias(r), Bias(x));
  return r;
}

INTERVAL ArCosh (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasArCosh (Bias(r), Bias(x));
  return r;
}

INTERVAL ArTanh (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasArTanh (Bias(r), Bias(x));
  return r;
}

INTERVAL ArCoth (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasArCoth (Bias(r), Bias(x));
  return r;
}

INTERVAL Exp (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasExp (Bias(r), Bias(x));
  return r;
}

INTERVAL Log (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasLog (Bias(r), Bias(x));
  return r;
}

INTERVAL Log10 (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasLog10 (Bias(r), Bias(x));
  return r;
}

INTERVAL IAbs (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasIAbs (Bias(r), Bias(x));
  return r;
}

INTERVAL Sqr (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasSqr (Bias(r), Bias(x));
  return r;
}

INTERVAL Sqrt (CONST INTERVAL & x)
{
  INTERVAL r;

  BiasSqrt (Bias(r), Bias(x));
  return r;
}

INTERVAL Root (CONST INTERVAL & x, INT n)
{
  INTERVAL r;

  BiasRoot (Bias(r), Bias(x), n);
  return r;
}

INTERVAL Power (CONST INTERVAL & x, INT n)
{
  INTERVAL r;

  BiasPowerN (Bias(r), Bias(x), n);
  return r;
}

INTERVAL Power (CONST INTERVAL & x, CONST INTERVAL & y)
{
  INTERVAL r;

  BiasPowerI (Bias(r), Bias(x), Bias(y));
  return r;
}
