C****************************************************************************
C     Driver file of the minf_tides program
C     This file has been created by MathTIDES (2.00) January 26, 2011, 18:16
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


      Program  dr_minf_lorenz
      IMPLICIT NONE
      INTEGER  i,j
C --- NUMBER OF VARIABLES AND PARAMETERS
      INTEGER  NVAR,NPAR
      PARAMETER  (NVAR = 3)
      PARAMETER  (NPAR = 3)
C --- TOLERANCES
      REAL*8 tolabs,tolrel
C --- TIMES: INITIAL, FINAL, INCREMENT
      REAL*8 tini, tend, dt
C --- VARIABLES AND PARAMETERS
      REAL*8 v(NVAR)
      REAL*8 x(NVAR)
      REAL*8 p(NPAR)
C --- FILE NAME OF DENSE AND COEFFICIENTS OUTPUT 
      CHARACTER ofname*20, cfname*20
C --- OPTIONS
      LOGICAL defect_error_control
C --- COUNTERS
      INTEGER accepted_steps, rejected_steps
C --- CONSTANTS OF THE METHOD (safety factors, maximum order, ...)
      REAL*8 fac1,fac2,fac3,rminstep,rmaxstep, aux, error
      INTEGER nitermax,nordinc,minord,maxord
C --- GLOBALS
      COMMON /OPT/ defect_error_control
      COMMON /ARS/ accepted_steps, rejected_steps
      COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
      COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
      COMMON /NFILES/ ofname, cfname




C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C     INITIAL CONDITIONS,  INTEGRATION TIMES, TOLERANCES
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------

C --- PARAMETERS VALUE
      p(1) = 10.d0
      p(2) = 28.d0
      p(3) = 2.666666666666667d0

C --- INITIAL VALUES
      v(1) = -13.7636106821342d0
      v(2) = -19.5787519424518d0
      v(3) = 27.d0

C --- INITIAL INTEGRATION POINT
      tini = 0.d0

C --- ENDPOINT OF INTEGRATION
      tend = 1.558652210716175d0

C --- DELTA t FOR DENSE OUTPUT
      dt   = 1.558652210716175d0


C --- REQUIRED TOLERANCES
      tolrel = 1.d-16
      tolabs = 1.d-16


C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
C       CALL THE INTEGRATOR
C-------------------------------------------------------------------------------
C-------------------------------------------------------------------------------
      DO 10 i = 1,NVAR
        x(i) = v(i)
 10   CONTINUE

      CALL minf_tides(v,NVAR,p,NPAR,tini,tend,dt,
     &   tolrel,tolabs)



      aux = 0.d0
      error = -1.d0

      DO 20 i = 1,NVAR
        aux = ABS(v(i)-x(i))
        IF (aux .gt. error) error = aux
 20   CONTINUE

      IF (error .lt. 1.d-10) THEN
       CALL exit(0) 
       STOP
      ELSE
       write (*,*) "Error in test_minf_lorenz = ", error
       CALL exit(1) 
       STOP
      ENDIF


      STOP
      END




