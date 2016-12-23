/************************************************************************
 *
 * Definition of type INTERVAL_VECTOR
 * ----------------------------------
 *
 * Copyright (C) 1993, 1995, 1997 Olaf Knueppel
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
 * $Id: IntervalVector.h 478 2006-08-09 13:13:30Z keil $
 *
 ************************************************************************/

#ifndef __INTERVAL_VECTOR__
#define __INTERVAL_VECTOR__

#include <BIAS/Bias1.h>
#include <Configuration.h>
#include <Error.h>
#include <Vector.h>
#include <Interval.h>

class INTERVAL_VECTOR {
  INT nElements;
#ifdef __DONTCOPY__
  INT IsTemporary;
#endif
public:
  BIASINTERVAL *theElements;
  INTERVAL_VECTOR ()
    {
#ifdef __DONTCOPY__
      IsTemporary = 0;
#endif
      nElements = 0; theElements = NULL; }
  INTERVAL_VECTOR (INT n)
    {
#ifdef __DONTCOPY__
      IsTemporary = 0;
#endif
      CheckVectorRange (n);
      nElements = n;
      theElements = new BIASINTERVAL[n];
      CheckMemory (theElements);
    }
  ~INTERVAL_VECTOR ()
    { if (theElements != NULL) delete [] theElements; }
  INTERVAL_VECTOR (CONST INTERVAL_VECTOR &);
  INTERVAL & operator () (INT i) CONST
    { CheckVectorIndex (i, *this); return * ((INTERVAL *)(theElements + i-1)); }
  INTERVAL_VECTOR & operator =  (CONST INTERVAL_VECTOR &);
  INTERVAL_VECTOR & operator =  (CONST VECTOR &);
  INTERVAL_VECTOR & operator += (CONST INTERVAL_VECTOR &);
  INTERVAL_VECTOR & operator -= (CONST INTERVAL_VECTOR &);
  INTERVAL_VECTOR & operator *= (CONST INTERVAL &);
  INTERVAL_VECTOR & operator /= (CONST INTERVAL &);
  INTERVAL_VECTOR & operator += (CONST VECTOR &);
  INTERVAL_VECTOR & operator -= (CONST VECTOR &);
  INTERVAL_VECTOR & operator *= (REAL);
  INTERVAL_VECTOR & operator /= (REAL);
  friend INT Dimension (CONST INTERVAL_VECTOR &);
  friend VOID Resize (INTERVAL_VECTOR &, INT);
  friend VOID MakeTemporary (INTERVAL_VECTOR &);
  friend VOID MakePermanent (INTERVAL_VECTOR &);
  friend VOID Clear (INTERVAL_VECTOR &);
  friend VOID Initialize (INTERVAL_VECTOR &, CONST INTERVAL &);
  friend VECTOR Inf (CONST INTERVAL_VECTOR &);
  friend VECTOR Sup (CONST INTERVAL_VECTOR &);
  friend VECTOR Mid (CONST INTERVAL_VECTOR &);
  friend VECTOR Diam (CONST INTERVAL_VECTOR &);
  friend VECTOR Abs (CONST INTERVAL_VECTOR &);
  friend INT Intersection (INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INT operator <= (CONST VECTOR &         , CONST INTERVAL_VECTOR &);
  friend INT operator <= (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INT operator <  (CONST VECTOR &         , CONST INTERVAL_VECTOR &);
  friend INT operator <  (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INT operator == (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INT operator != (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR Hull (CONST VECTOR &);
  friend INTERVAL_VECTOR Hull (CONST VECTOR &, CONST VECTOR &);
  friend INTERVAL_VECTOR Hull (CONST VECTOR &, CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR Hull (CONST INTERVAL_VECTOR &, CONST VECTOR &);
  friend INTERVAL_VECTOR Hull (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR SymHull (CONST VECTOR &);
  friend INTERVAL_VECTOR operator + (CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR operator - (CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR operator + (CONST VECTOR &         , CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR operator + (CONST INTERVAL_VECTOR &, CONST VECTOR &         );
  friend INTERVAL_VECTOR operator + (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR AddBounds  (CONST VECTOR &, CONST VECTOR &);
  friend INTERVAL_VECTOR operator - (CONST VECTOR &         , CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR operator - (CONST INTERVAL_VECTOR &, CONST VECTOR &         );
  friend INTERVAL_VECTOR operator - (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR SubBounds  (CONST VECTOR &, CONST VECTOR &);
  friend INTERVAL_VECTOR operator * (REAL                   , CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR operator * (CONST INTERVAL &       , CONST VECTOR &         );
  friend INTERVAL_VECTOR operator * (CONST INTERVAL &       , CONST INTERVAL_VECTOR &);
  friend INTERVAL_VECTOR MulBounds  (REAL, CONST VECTOR &);
  friend INTERVAL_VECTOR operator / (CONST INTERVAL_VECTOR &, REAL                   );
  friend INTERVAL_VECTOR operator / (CONST VECTOR &         , CONST INTERVAL &       );
  friend INTERVAL_VECTOR operator / (CONST INTERVAL_VECTOR &, CONST INTERVAL &       );
  friend INTERVAL_VECTOR DivBounds  (CONST VECTOR &, REAL);
  friend INTERVAL operator * (CONST VECTOR &         , CONST INTERVAL_VECTOR &);
  friend INTERVAL operator * (CONST INTERVAL_VECTOR &, CONST VECTOR &         );
  friend INTERVAL operator * (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
  friend INTERVAL MulBounds  (CONST VECTOR &, CONST VECTOR &);
  friend INTERVAL Sqr  (CONST INTERVAL_VECTOR &);
  friend INTERVAL Norm (CONST INTERVAL_VECTOR &);
  friend ostream & operator << (ostream &, CONST INTERVAL_VECTOR &);
  friend istream & operator >> (istream &, INTERVAL_VECTOR &);
};

inline INT Dimension (CONST INTERVAL_VECTOR & x) {
  return x.nElements; 
}

VOID Resize (INTERVAL_VECTOR &, INT);

#ifdef __DONTCOPY__
inline VOID MakeTemporary (INTERVAL_VECTOR & x) { 
  x.IsTemporary = 1; 
}

inline VOID MakePermanent (INTERVAL_VECTOR & x) { 
  x.IsTemporary = 0; 
}
#else
inline VOID MakeTemporary (INTERVAL_VECTOR & x) { }

inline VOID MakePermanent (INTERVAL_VECTOR & x) { }
#endif

VOID Clear (INTERVAL_VECTOR &);
VOID Initialize (INTERVAL_VECTOR &, CONST INTERVAL &);
VECTOR Inf (CONST INTERVAL_VECTOR &);
VECTOR Sup (CONST INTERVAL_VECTOR &);
VECTOR Mid (CONST INTERVAL_VECTOR &);
VECTOR Diam (CONST INTERVAL_VECTOR &);
VECTOR Abs (CONST INTERVAL_VECTOR &);
INT Intersection (INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
INT operator <= (CONST VECTOR &, CONST INTERVAL_VECTOR &);
INT operator <= (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
INT operator <  (CONST VECTOR &, CONST INTERVAL_VECTOR &);
INT operator <  (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
INT operator == (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
INT operator != (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR Hull (CONST VECTOR &);
INTERVAL_VECTOR Hull (CONST VECTOR &, CONST VECTOR &);
INTERVAL_VECTOR Hull (CONST VECTOR &, CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR Hull (CONST INTERVAL_VECTOR &, CONST VECTOR &);
INTERVAL_VECTOR Hull (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);

inline INTERVAL_VECTOR SymHull (CONST VECTOR & x) {
  return Hull(-x, x);
}

inline INTERVAL_VECTOR operator + (CONST INTERVAL_VECTOR & a) { 
  return (INTERVAL_VECTOR &) a; 
}

INTERVAL_VECTOR operator - (CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR operator + (CONST VECTOR &, CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR operator + (CONST INTERVAL_VECTOR &, CONST VECTOR &);
INTERVAL_VECTOR operator + (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR AddBounds  (CONST VECTOR &, CONST VECTOR &);
INTERVAL_VECTOR operator - (CONST VECTOR &, CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR operator - (CONST INTERVAL_VECTOR &, CONST VECTOR &);
INTERVAL_VECTOR operator - (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR SubBounds  (CONST VECTOR &, CONST VECTOR &);
INTERVAL_VECTOR operator * (REAL, CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR operator * (CONST INTERVAL &, CONST VECTOR &);
INTERVAL_VECTOR operator * (CONST INTERVAL &, CONST INTERVAL_VECTOR &);
INTERVAL_VECTOR MulBounds  (REAL, CONST VECTOR &);
INTERVAL_VECTOR operator / (CONST INTERVAL_VECTOR &, REAL);
INTERVAL_VECTOR operator / (CONST VECTOR &, CONST INTERVAL &);
INTERVAL_VECTOR operator / (CONST INTERVAL_VECTOR &, CONST INTERVAL &);
INTERVAL_VECTOR DivBounds  (CONST VECTOR &, REAL);
INTERVAL operator * (CONST VECTOR &, CONST INTERVAL_VECTOR &);
INTERVAL operator * (CONST INTERVAL_VECTOR &, CONST VECTOR &);
INTERVAL operator * (CONST INTERVAL_VECTOR &, CONST INTERVAL_VECTOR &);
INTERVAL MulBounds  (CONST VECTOR &, CONST VECTOR &);
INTERVAL Sqr  (CONST INTERVAL_VECTOR &);
INTERVAL Norm (CONST INTERVAL_VECTOR &);
ostream & operator << (ostream &, CONST INTERVAL_VECTOR &);
istream & operator >> (istream &, INTERVAL_VECTOR &);
  
#endif /* __INTERVAL_VECTOR__ */
