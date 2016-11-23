C****************************************************************************
C     This file has been created by MathTIDES (2.00) April 16, 2011, 20:06
C       
C     Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
C     Grupo de Mecanica Espacial
C     University of Zaragoza
C     SPAIN
C       
C     http://gme.unizar.es/software/tides
C     Contact: <tides@unizar.es>
C        
C     This file is part of TIDES.
C        
C     TIDES is free software: you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published by
C     the Free Software Foundation, either version 3 of the License, or
C     (at your option) any later version.
C        
C     TIDES is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C        
C     You should have received a copy of the GNU General Public License
C     along with TIDES.  If not, see <http://www.gnu.org/licenses/>.
C        
C*****************************************************************************


      SUBROUTINE minfseries(t,v,NVAR,p,NPAR,XVAR,ORDER,MO)
      IMPLICIT   NONE
      INTEGER    i,inext,j,ORDER,NVAR,NPAR,TT,MO
      REAL*8     XVAR(0:MO,0:NVAR)
      REAL*8     XX(0:MO,0:14)
      REAL*8     v(NVAR)
      REAL*8     p(NPAR)
      REAL*8     t

      REAL*8     pr(4)
C-------------------------------------------------------------------------------
      REAL*8     mul_mf
      REAL*8     div_mf
      REAL*8     inv_mf
      REAL*8     exp_mf
      REAL*8     pow_mf_c
      REAL*8     sin_mf
      REAL*8     cos_mf
      REAL*8     log_mf
C-------------------------------------------------------------------------------
      DO j = 1, NPAR
          pr(j) = p(j)
      END DO
      pr(4) = -1.d0 * pr(3)
      TT = 14
      DO j=0, TT
          DO i=0, ORDER
              XX(i,j) = 0.d0
          END DO
      END DO
      XX(0,0) = t
      XX(1,0) = 1.d0
      DO j = 1, NVAR
          XX(0,j) = v(j)
      END DO
C-------------------------------------------------------------------------------
      DO i=0, ORDER-1
          XX(i,4) = -1.d0* XX(i,1)
          XX(i,5) = -1.d0* XX(i,2)
          XX(i,6) = pr(2)* XX(i,1)
          XX(i,7) = mul_mf(1,2,i,XX,TT,MO)
          XX(i,8) = XX(i,4)+XX(i,2)
          XX(i,9) = XX(i,5)+XX(i,6)
          XX(i,10) = pr(4)* XX(i,3)
          XX(i,11) = mul_mf(4,3,i,XX,TT,MO)
          XX(i,12) = XX(i,7)+XX(i,10)
          XX(i,13) = XX(i,9)+XX(i,11)
          XX(i,14) = pr(1)* XX(i,8)
          inext = i + 1
          XX(inext,1) = XX(i,14)/inext
          XX(inext,2) = XX(i,13)/inext
          XX(inext,3) = XX(i,12)/inext
      END DO
C-------------------------------------------------------------------------------
      DO j=0, NVAR
          DO i=0, ORDER
              XVAR(i,j) = XX(i,j)
          END DO
      END DO
      RETURN
      END SUBROUTINE
C-------------------------------------------------------------------------------