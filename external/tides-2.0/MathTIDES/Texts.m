(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`Texts`*)


(* ::Text:: *)
(*Copyright (C) 2010  Alberto Abad, Roberto Barrio, Fernando Blesa, Marcos Rodriguez*)
(*Grupo de Mec\[AAcute]nica Espacial.  IUMA.*)
(*University of Zaragoza*)
(*50009 Zaragoza. Spain.*)
(**)
(*http://gme.unizar.es/software/tides*)


(* ::Text:: *)
(*This file is part of TIDES.*)
(*  	*)
(*TIDES is free software : you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.*)
(*  	*)
(*TIDES is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.*)
(*  	*)
(*You should have received a copy of the GNU General Public License along with TIDES.  If not, see < http://www.gnu.org/licenses/ > .*)


(* ::Text:: *)
(*v - 20*)


(* ::Title::Closed:: *)
(*Texts*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`Texts`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
mincgen,
minhgen,
minfgen,
headstdDP,
headstdMP
}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`Texts`*"]
Clear @@ Names["MathTIDES`Texts`*"]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsection::Closed:: *)
(*C*)


mincgen ="
/****************************************************************************
    
    minc_tides: kernel of the C Minimal version of TIDES

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

#include \"minc_tides.h\" 


double    *v, *p, **XVAR, **XVAR2, TOL_ABS, TOL_REL, tzero, deltat; 
int       VARS, PARS, MAX_ORDER, tflag;
double    fac1=0.95e0,fac2=10.e0,fac3=0.8e0 ; 
double    rmaxstep=1.e2,rminstep=1.e-2; 
int       nitermax=5, nordinc=5, minord=6, maxord=26; 
int       defect_error_control = 0;
int       dense_output = -1, coef_output = 0;
int       accepted_steps = 0, rejected_steps = 0;
int       ipos = 1;
char      ofname[20]=\"\", cfname[20]=\"\";
FILE      *fd, *fc;



/************************************************************************/ 

void minc_tides(double *var, int nvar, double *par, int npar,  double tini, double tend, double dt,
                double tolrel, double tolabs) 
{ 
	double    tol, tolo;
	volatile  double  temp, extra;
	double    t0,t1,t2, step, nstep;
	int       i, j, ORDER; 

	
	VARS = nvar;
	PARS = npar; 
	MAX_ORDER = maxord;
	TOL_REL = tolrel;
	TOL_ABS = tolabs;
	tzero = tini;
	deltat = dt;
	v = var; 
	p = par;
	
	declare_matrix_coefs_mc();

	t0 = tini; 
	t1= 0.e0;
	t2 = 0.e0; 

	extra = 0.e0; 
	
	if(strcmp(ofname,\"\") == 0) 
		dense_output = 0;
	else if (strcmp(ofname,\"Screen\") == 0 || strcmp(ofname,\"screen\") == 0)
		fd = stdout;
	else 
		fd = fopen(ofname, \"w\");

	if(strcmp(cfname,\"\") != 0) {
		coef_output = -1;
		fc = fopen(cfname, \"w\");
	}
	
	if(t0 < tend) {
		tflag = 0;
		if(deltat < 0) deltat = -deltat;
	} else {
		tflag = 1; 
		if(deltat > 0) deltat = -deltat;
	}

	if(dense_output) {
		fprintf(fd, \"%25.15le \", tini ); 
		for(i = 0; i < VARS; i++) fprintf(fd, \"%25.15le \" , v[i]); 
		fprintf(fd, \"\\n\"); 	
	}

	while (((t0 < tend) && (tflag == 0)) ||
		   ((t0 > tend) && (tflag == 1))  ) { 
		tolerances_mc(&tol, &tolo, &ORDER);
		mincseries(t0, v, p, XVAR, ORDER, MAX_ORDER); 
		step = steps_mc(ORDER, tol);
		if(defect_error_control) steps_DEC_mc(t0, tol, ORDER, &step);
		if(coef_output)  fprintf(fc, \"%d %25.15le \", ORDER, t0); 
		temp = t0; 
		nstep = step + extra; 
		t0 = temp + nstep; 
		extra = (temp-t0) +nstep; 
		if(((t0 > tend) && (tflag == 0)) ||
		   ((t0 < tend) && (tflag == 1)))  nstep = (tend-temp); 
		accepted_steps +=1; 
		if(dense_output) dense_output_mc(temp, nstep, ORDER);
		horner_mc(v, nstep, ORDER) ; 
		if(coef_output)  {
		 	fprintf(fc, \"%25.15le \\n\", nstep); 
		 	for(i= 1; i <=VARS; i++) {
		 		fprintf(fc, \"%d \", i); 
		 		for(j = 0; j <= ORDER; j++) fprintf(fc, \"%25.15le \", XVAR[j][i]); 
		 		fprintf(fc, \" \\n\"); 
		 	}
			
		}

	} 
    ipos = 1;
    if(dense_output && strcmp(ofname,\"Screen\") != 0 && strcmp(ofname,\"screen\") != 0) {
		fclose(fd);	
	}
	if(coef_output)  fclose(fc);	


} 

/************************************************************************/ 
void dense_output_mc(double t0, double step, int ORDER)
{
	int i;
	double tend, ti, tit, vh[VARS]; 
	tend = t0 + step;
	ti = tzero + ipos*deltat;
	tit = ti - t0;
	while(((tit <= step) && (tflag == 0)) ||
		   ((tit >= step) && (tflag == 1))  ) {
		horner_mc(vh, tit, ORDER); 
		fprintf(fd, \"%25.15le \", ti ); 
		for(i = 0; i < VARS; i++) fprintf(fd, \"%25.15le \" , vh[i]); 
		fprintf(fd, \"\\n\");
		ipos++;
		ti = tzero + ipos*deltat;
		tit = ti -t0;
	}
}

/************************************************************************/ 

void horner_mc(double *v, double t, int ORDER) 
{ 
	int i,j; 
	double temp; 
	for(j = 1; j <= VARS; j++) { 
		temp = XVAR[ORDER][j] * t; 
		for(i = ORDER-1; i >= 1; i--) { 
			temp = t * (temp + XVAR[i][j]);
		} 
		v[j-1] = temp + XVAR[0][j]; 
	} 
} 

void hornerd_mc(double *v, double t, int ORDER) 
{ 
	int i,j; 
	double temp; 
	for(j = 1; j <= VARS; j++) { 
		temp = ORDER * XVAR[ORDER][j] * t; 
		for(i = ORDER-1; i > 1; i--)	
			temp = t * (temp + i * XVAR[i][j]); 
		v[j-1]= temp + XVAR[1][j]; 
	} 
} 
/************************************************************************/ 
double norm_inf_vec_mv(void)
{
	int i; 
	double norm;
	norm = fabs(v[0]);
	for(i = 1; i < VARS; i++) 
		norm = max_d(norm, fabs(v[i]));
	return norm;
}

double norm_inf_mat_mc(int ord)
{
	int i; 
	double norm;
	norm = 0.e0;
	for(i = 1; i <= VARS; i++) {
		norm = max_d(norm, fabs(XVAR[ord][i]));
	}
	return norm;
}
/************************************************************************/ 

void tolerances_mc(double *tol, double *tolo, int *ORDER)
{
	static double ynb =0.e0;
	double yna, miny;
	yna = norm_inf_vec_mv();
	*tol = TOL_ABS + max_d(yna,ynb) * TOL_REL;
	miny = min_d(yna,ynb);
	if(miny > 0.) 
		*tolo = min_d(TOL_ABS/miny, TOL_REL); 
	else
		*tolo = min_d(TOL_ABS, TOL_REL);	
	*ORDER = min_i(MAX_ORDER, floor(-log(*tolo)/2)+nordinc); 
	*ORDER = max_i(minord, *ORDER); 
	ynb = yna;
}

double steps_mc(int ORDER, double tol)
{
	static double stepant =0.e0;
	double ynu, ynp, tu, tp, step;
	int ord, orda, ordp;
	double  dord, dorda, dordp, rstep;
	
	ord = ORDER+1;
	do {
		ord--;
		ynu = norm_inf_mat_mc(ord);
	} while (ynu == 0.e0 && ord > 0) ;
	
	if(ord !=0 ) {
		orda = ord-1;
		ordp = ord+1;
		dord = 1.e0/ord;
		dorda = 1.e0/orda;
		dordp = 1.e0/ordp;
		ynp = norm_inf_mat_mc(orda);		
		if(ynp == 0.e0 ) 
			step = pow(tol, dordp) * pow(1.e0/ynu, dord); 
		else {	
			tp = pow(tol, dord) * pow(1.e0/ynp, dorda);
			tu = pow(tol, dordp) * pow(1.e0/ynu, dord); 
			step = min_d(tp,tu); 
		} 
		if(stepant != 0.e0) {
			rstep = step/stepant;
			if( rstep > rmaxstep )
				step = rmaxstep * stepant;
			else if( rstep <rminstep ) 
				step = rmaxstep * stepant;
	    }
		step *= fac1;
	} else {
		printf(\"*********Error*********\\n\");
		abort();
	}
	if(tflag ==1) step = -step;
	return step;	
}

void steps_DEC_mc(double t0, double tol, int ORDER, double *step)
{
	int iter, i;
	double norma, vh[VARS], vdh[VARS], t; 
	iter = 1; 
	norma = 1.e99;  
	t = *step;
	while((norma > fac2*tol) && (iter < nitermax)) { 
		horner_mc(vh, t, ORDER);
		hornerd_mc(vdh, t, ORDER);
		mincseries(t0+t, vh, p, XVAR2, 1, 1); 
		norma = fabs(XVAR2[1][1]-vdh[0]); 
		for(i=2; i <=VARS; i++) 
			norma +=max_d(fabs(XVAR2[1][i]-vdh[i-1]), norma); 
		if(iter >1) {
			t = fac3 *t; 
			rejected_steps +=1 ;
		}
		iter++; 
	} 
	if (iter >= nitermax) {
		printf(\"The maximum number of iterations in defect error control has been reached.\\n\");
		printf(\"The program continues with a stepsize that does not achieve the defect error control condition.\\n\");
	}
	*step = t;
}
/************************************************************************/ 
void	declare_matrix_coefs_mc(void)
{
	int i,j;
	XVAR  = (double **) calloc(MAX_ORDER+1,sizeof(double *));
	XVAR2 = (double **) calloc(2,sizeof(double *));
	
	for( i=0; i<=MAX_ORDER; i++ )
		XVAR[i] = (double *) calloc(VARS+1, sizeof(double));
	for( i=0; i< 2; i++ )
		XVAR2[i] = (double *) calloc(VARS+1, sizeof(double));
	
	for( i=0; i<=MAX_ORDER; i++ )
		for( j=0; j<=VARS; j++ )
			XVAR[i][j] = 0; 
	for( i=0; i< 2; i++ )
		for( j=0; j<=VARS; j++ )
			XVAR2[i][j] = 0; 
}

/************************************************************************/ 

double mul_mc(double* u, double* v, int k) 

{ 
	int j; 
	double w = 0.e0; 
	for(j = 0; j <= k; j++) 
		w += (u[j] * v[k-j]); 
	return w; 
} 

double div_mc(double* u, double* v, double* w, int k) 
{ 

	int j; 
	double ww; 
	if(v[0] == 0.e0) { 
		printf(\"*** Function div_mc found division by zero ***\\n\");
		exit(EXIT_FAILURE); 
	} 
	ww = u[k]; 
	for(j = 1; j <= k; j++) ww -= (v[j] *w[k-j]);
	ww /= v[0]; 
	return ww; 
} 

double inv_mc(double p, double* u, double* w, int k) 
{ 
	
	int j; 
	double ww = 0.e0; 
	if(u[0] == 0.e0) { 
		printf(\"*** Function inv_mc found division by zero ***\\n\");
		exit(EXIT_FAILURE); 
	} 
	if(k == 0) 
		ww = 1.e0; 
	else  
		for(j = 0; j < k; j++) ww -= (u[k-j] *w[j]);
	ww /= u[0]; 
	return ww*p; 
} 

double exp_mc(double* u, double* v, int k) 

{ 
	int j; 
	double w; 
	if(k == 0) 
		w = exp(u[0]); 
	else { 
		w = k *v[0]*u[k]; 
		for(j = 1; j < k; j++) 
			w += ( (k-j) * v[j] * u[k-j] );
		w /= k; 
	} 
	return w; 
} 

double pow_mc_c(double* u, double c, double* w, int k) 
{ 
	int j; 
	double ww; 
	if(k == 0) { 
		if(u[0] == 0.e0) { 
			printf(\"*** Function pow_mc _c found division by zero ***\\n\");
			exit(EXIT_FAILURE); 
		} 
		ww = pow(u[0], c);  
	} else { 
		if(u[0] != 0.e0) {
			ww = c * k * w[0] * u[k]; 
			for(j = 1; j < k; j++) 
				ww += ( (c * (k - j) -j ) * w[j] * u[k-j]); 
			ww /= (k * u[0]); 
		} else ww = 0.e0;   
	} 
	return ww; 
} 

double log_mc(double* u, double* w, int k) 
{ 
	int j; 
	double ww; 
	if (k == 0) { 
		if(u[0] == 0.e0) { 
			printf(\"*** Function log_mc found division by zero ***\\n\");
			exit(EXIT_FAILURE); 
		} 
		ww = log(u[0]); 
	} else { 
		ww=k*u[k]; 
		for(j = 1; j < k; j++) 
			ww -= ((k - j) * u[j] * w[k-j]); 
		ww /= (k * u[0]); 
	} 
	return ww; 
} 

double sin_mc(double* u, double* v, int k) 
{
	int j; 
	double ww; 
	if (k == 0) { 
		ww = sin(u[0]); 
	} else { 
		ww = u[1]*v[k-1];  
		for(j = 2; j <= k; j++) 
			ww += ( j * u[j] * v[k-j] ); 
		ww /= k; 
	} 
	return ww; 
} 

double cos_mc(double* u, double* v, int k) 

{ 

	int j; 
	double ww; 
	if(k == 0) { 
		ww = cos(u[0]); 
	} else { 
		ww = -u[1]*v[k-1]; 
		for(j = 2; j <= k; j++) 
			ww -= ( j * u [j] * v [k-j] ); 
		ww /=   k;   
	}     
	return ww;     
} 

"



(* ::Subsection::Closed:: *)
(*H*)


minhgen = 
"
/****************************************************************************
    
    minc_tides: kernel of the C Minimal version of TIDES

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

#ifndef minc_tides_HeadFile 
#define minc_tides_HeadFile 

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <string.h> 

double mul_mc(double* u,   double* v,   int k);     
double div_mc(double* u,   double* v,   double* w,int k);   
double inv_mc(double p, double* u, double* w, int k); 
double exp_mc(double* u,   double* v,   int k);     
double pow_mc_c(double* u,   double e,   double* w,   int k);   
double log_mc(double* u,   double* w,   int k);     
double sin_mc(double* u,   double* v,   int k);   
double cos_mc(double* u,   double* v,   int k);

void	horner_mc( double *v, double t, int ORDER) ;
void	hornerd_mc(double *v, double t, int ORDER) ; 
void	dense_output_mc(double t0, double step, int ORDER);

double	norm_inf_vec_mv(void);
double	norm_inf_mat_mc(int ord);
void	declare_matrix_coefs_mc(void);
void	tolerances_mc(double *tol, double *tolo, int *ORDER);
double	steps_mc(int ORDER, double tol);
void	steps_DEC_mc(double t0, double tol, int ORDER, double *step);


void	mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO); 
void	minc_tides(double *var, int nvar, double *par, int npar,  double tini, double tend, double dt,
                double tol_rel, double tol_abs); 


static inline int min_i(int a, int b) { 
	return a < b ? a: b; 
} 

static  inline int max_i(int a, int b) { 
	return a > b ? a: b; 
} 

static  inline double min_d(double a, double b) { 
	return a < b ? a: b; 
} 

static  inline double max_d(double a, double b) { 
	return a > b ? a: b; 
} 


#endif 



";


(* ::Subsection::Closed:: *)
(*Fortran*)


minfgen = "
C****************************************************************************
C       
C     minf_tides: kernel of the Fortran Minimal version of TIDES
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

      BLOCKDATA CONSTMETHOD
      REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
      INTEGER nitermax,nordinc,minord,maxord
      INTEGER accepted_steps, rejected_steps
      INTEGER IPOS
      LOGICAL defect_error_control
      LOGICAL dense_output, coef_output
      CHARACTER ofname*20, cfname*20
      COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
      COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
      COMMON /OPT/ defect_error_control
      COMMON /OPT2/ dense_output, coef_output
      COMMON /ARS/ accepted_steps, rejected_steps
      COMMON /DOPOS/IPOS
      COMMON /NFILES/ ofname, cfname
	  DATA fac1,fac2,fac3/0.95d0,10.d0,0.8d0/
	  DATA rminstep,rmaxstep/1.0d2,1.0d-2/
	  DATA nitermax,nordinc/5,5/
	  DATA minord,maxord/6,26/
      DATA defect_error_control/.FALSE./
      DATA dense_output,coef_output/.TRUE.,.FALSE./
      DATA accepted_steps, rejected_steps/0,0/
      DATA IPOS/1/
      DATA ofname, cfname/'',''/
      END


      SUBROUTINE minf_tides(v,numvar,p,numpar,tini,tend,dt,
     &   tolrel,tolabs)

        IMPLICIT NONE
        LOGICAL defect_error_control
        LOGICAL dense_output, coef_output
        CHARACTER ofname*20, cfname*20
        INTEGER accepted_steps, rejected_steps
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        INTEGER FL,FC, IPOS
        INTEGER numvar,numpar
        INTEGER I,J
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /VP/ NVAR, NPAR
        COMMON /OPT/ defect_error_control
        COMMON /OPT2/ dense_output, coef_output
        COMMON /ARS/ accepted_steps, rejected_steps
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /DOPOS/IPOS
        COMMON /FILES/ FL,FC
        COMMON /NFILES/ ofname, cfname
        REAL*8 v(numvar),p(numpar),XVAR(0:maxord,0:numvar)
        REAL*8 t0,tini,tend,dt,tolrel,tolabs
        REAL*8 tol,tolo, step, temp, nstep, extra
        INTEGER ORDER
        INTEGER tflag

        TZERO = tini
        t0 = tini
        DELTAT = dt
        TOL_REL = tolrel
        TOL_ABS = tolabs
        NVAR = numvar
        NPAR = numpar
        extra = 0.d0

        IF (t0 .LT. tend) THEN
            tflag = 0
            IF (DELTAT .LT. 0.d0) THEN
                DELTAT = -DELTAT
            END IF
        ELSE
            tflag = 1
            IF (DELTAT .GT. 0.d0) THEN
                DELTAT = -DELTAT
            END IF
        END IF
        
        
        IF(ofname == '') THEN
       		dense_output = .FALSE.
       	ELSE IF (ofname == 'Screen' .OR. ofname == 'screen') THEN
            FL = 6
       	ELSE 
            FL = 72
            OPEN (UNIT = FL, FILE = ofname, STATUS = 'UNKNOWN')
       	END IF
        
        IF(cfname .NE. '') THEN
       		coef_output = .TRUE.
            FC = 72
            OPEN (UNIT = FC, FILE = cfname, STATUS = 'UNKNOWN')
       	END IF
        

        IF(dense_output) THEN
       		WRITE(FL,'(90E25.16)') TZERO, v
        END IF
           
        DO WHILE(((t0 .LT. tend).AND.(tflag.EQ.0)) .OR.
     &           ((t0 .GT. tend).AND.(tflag.EQ.1)))
            CALL tolerances_mf(v, tol,tolo, ORDER)
            CALL minfseries(t0,v,NVAR,p,NPAR,XVAR,ORDER,maxord)
            CALL steps_mf(tol, XVAR, ORDER, step, tflag)

            IF (defect_error_control) THEN
                CALL steps_DEC_mf(t0,step,tolo,p,XVAR,ORDER)
            END IF
 
            temp = t0
            nstep = step + extra
            t0 = temp + nstep
            extra = (temp-t0)+nstep
            IF(((t0 .GT. tend).AND.(tflag.EQ.0)) .OR.
     &           ((t0 .LT. tend).AND.(tflag.EQ.1))) THEN
                nstep  = (tend-temp)
            END IF
            accepted_steps =  accepted_steps + 1
            IF(dense_output) THEN
                CALL dense_output_mf(temp,nstep,XVAR,ORDER,tflag)
            END IF
            CALL horner_mf(v,XVAR,ORDER,nstep)
         
            IF(coef_output ) THEN
       		  WRITE(FC,'(I2, 2E25.16)') ORDER, t0, nstep
       		  DO I = 1, NVAR
        		  WRITE(FC,'(I2, 50E25.16)') I, (XVAR(J,I), J=0,ORDER)
      		  END DO
            END IF
                
        END DO
           
        IPOS = 1
        IF(dense_output .AND. FL .NE. 6) THEN
            CLOSE(FL)
        END IF
        IF(coef_output ) THEN
            CLOSE(FC)
        END IF


        RETURN
      END SUBROUTINE
           
C--------------------------------------------------------------
C--------------------------------------------------------------
C     tolerances_mf
C--------------------------------------------------------------
C--------------------------------------------------------------

      SUBROUTINE tolerances_mf(v, tol, tolo, ORDER)       
        IMPLICIT NONE
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        INTEGER ORDER
        REAL*8 yna,ynb,tol,tolo,v(NVAR),miny
           DATA ynb /0.0d0/
           SAVE ynb
          
        CALL norm_inf_vec_mf(v, yna)
        tol = TOL_ABS + MAX(yna,ynb)*TOL_REL
        miny = MIN(yna,ynb)
        IF(miny .gt. 0.d0) THEN
             tolo = MIN(TOL_ABS/miny, TOL_REL)
        ELSE
             tolo = MIN(TOL_ABS, TOL_REL)
        END IF


        ORDER = MIN(maxord, int(-log(tolo)/2)+nordinc)
        ORDER = MAX(minord,ORDER)
        
        ynb = yna
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     steps_mf
C--------------------------------------------------------------
C--------------------------------------------------------------


      SUBROUTINE steps_mf(tol, XVAR, ORDER, step, tflag)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 XVAR(0:maxord,0:NVAR)
        REAL*8 tol, ynu,ynp,dord,dorda,dordp,sp,su
        INTEGER ORDER, ord, orda, ordp
        REAL*8 step, stepant, rstep
        DATA stepant /0.0d0/
        SAVE stepant
        INTEGER tflag
          
        ord = ORDER+1
        ynu = 0.d0
        DO WHILE((ynu .EQ. 0.d0) .AND. (ord .GT. 0))
           ord = ord-1
           CALL norm_inf_mat_mf(XVAR,ord, ynu)
        END DO
         
        IF(ord .NE. 0)  THEN
             orda = ORDER-1
             ordp = ORDER+1
             dord  = 1.d0/DFLOAT(ord)
             dorda = 1.d0/DFLOAT(orda)
             dordp = 1.d0/DFLOAT(ordp)
             CALL norm_inf_mat_mf(XVAR,orda, ynp)
             IF(ynp .eq. 0.d0) THEN
               step = (tol**dordp) * ((1.d00/ynu)**dord)
             ELSE
              sp = (tol**dord) * ((1.d00/ynp)**dorda)
              su = (tol**dordp) * ((1.d00/ynu)**dord)
              step = MIN(sp,su)
            END IF 
            IF(stepant. NE. 0.0d0) THEN
              rstep = step/stepant
              IF(rstep .GT. rmaxstep) THEN
                 step = rmaxstep*stepant
              ELSE IF(rstep .LT. rminstep) THEN
                 step = rminstep*stepant
              END IF
            END IF
            step = fac1*step
        ELSE
            WRITE(*,*) '*********Error*********'
            STOP
        END IF
        IF (tflag .EQ. 1) THEN
            step = -step
        END IF
        RETURN
      END SUBROUTINE

 

C--------------------------------------------------------------
C--------------------------------------------------------------
C     steps_DEC_mf
C--------------------------------------------------------------
C--------------------------------------------------------------


      SUBROUTINE steps_DEC_mf(t0,step,tolo,p,XVAR,ORDER)      
        IMPLICIT NONE
        INTEGER accepted_steps, rejected_steps
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER NVAR,NPAR
        REAL*8 fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /VP/ NVAR, NPAR
        COMMON /ARS/ accepted_steps, rejected_steps
        COMMON /CONSTMET1/ fac1,fac2,fac3,rminstep,rmaxstep
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 vh(NVAR),vdh(NVAR), NORM
        REAL*8 t0,step,tolo, p(NPAR)
        REAL*8 XVAR(0:maxord,0:NVAR), XVAR2(0:maxord,0:NVAR)
        INTEGER ITER, ORDER,I
c
        ITER = 1
        NORM = 1.d99
           
        DO WHILE ((NORM .GT. fac2*tolo) .AND. (ITER .LT. nitermax))
             
            CALL horner_mf(vh,XVAR,ORDER,step)
            CALL hornerd_mf(vdh,XVAR,ORDER,step)
            CALL minfseries(t0+step,vh,NVAR,p,NPAR,XVAR2,1,1)
            NORM = ABS(XVAR2(1,1)-vdh(1))
            DO I=2, NVAR
                    NORM = NORM + ABS(XVAR2(1,i)-vdh(i))
            END DO
            IF(ITER .GT.1) THEN
                rejected_steps = rejected_steps+1
                step = fac3*step
            END IF
            ITER = ITER+1
               
        END DO
        
        IF( ITER .GE. nitermax) THEN
        	WRITE(*,*) 'The maximum number of iterations in '
         	WRITE(*,*) 'defect_error_control has been reached.'
       		WRITE(*,*) 'The program continues with a stepsize that '
       		WRITE(*,*) 'does not achieve the defect_error_control '
       		WRITE(*,*) 'condition.'
		END IF
c
        RETURN
      END SUBROUTINE

         
C--------------------------------------------------------------
C--------------------------------------------------------------
C     norm_inf_vec_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE norm_inf_vec_mf(v, norm)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        COMMON /VP/ NVAR, NPAR
        REAL*8 norm
        REAL*8 v(NVAR)
        INTEGER I
        norm = 0.d0
        DO I =1,NVAR
             norm =MAX(ABS(v(I)), norm)
        END DO
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     norm_inf_mat_mf
C--------------------------------------------------------------
C--------------------------------------------------------------

      SUBROUTINE norm_inf_mat_mf(XVAR,ord, norm)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        COMMON /VP/ NVAR, NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        REAL*8 norm
        REAL*8 XVAR(0:maxord,0:NVAR)
        INTEGER I,ord
        norm = 0.e0
        DO I =1,NVAR
             norm =MAX(abs(XVAR(ord,I)), norm)
        END DO
        RETURN
      END SUBROUTINE


C--------------------------------------------------------------
C--------------------------------------------------------------
C     dense_output_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE dense_output_mf(t0,step,XVAR,ORDER,tflag)
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        INTEGER FL,FC,IPOS
        REAL*8 TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /TDO/ TZERO,DELTAT,TOL_REL,TOL_ABS
        COMMON /VP/ NVAR, NPAR
        COMMON /FILES/ FL,FC
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        COMMON /DOPOS/IPOS
        REAL*8 v(NVAR), XVAR(0:maxord,0:NVAR)
        REAL*8 t0,ti,tit,tend,step
        INTEGER  ORDER
        INTEGER tflag

        tend = t0 + step
        ti = TZERO + IPOS*DELTAT
        tit= ti    - t0
        DO WHILE(((tit .LE. step).AND.(tflag.EQ.0)) .OR.
     &           ((tit .GE. step).AND.(tflag.EQ.1)))
            CALL horner_mf(v,XVAR,ORDER,tit)
            WRITE(FL,'(90E25.16)') ti, v
            IPOS =IPOS + 1
            ti = TZERO + IPOS*DELTAT
            tit= ti    - t0
        END DO
        RETURN
      END SUBROUTINE
C--------------------------------------------------------------
C--------------------------------------------------------------
C     horner_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE horner_mf(v,XVAR,ORDER,t)       
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        INTEGER ORDER, i, nd
        REAL*8 XVAR(0:maxord,0:NVAR), t, au
        REAL*8 v(NVAR)
        DO nd=1, NVAR
          au = XVAR(ORDER,nd)*t
          DO i = ORDER-1, 1, -1
            au = (au+XVAR(i,nd))*t
          END DO
          v(nd) = au+XVAR(0,nd)
        END DO
        RETURN
      END SUBROUTINE

C--------------------------------------------------------------
C--------------------------------------------------------------
C     hornerd_mf 
C--------------------------------------------------------------
C--------------------------------------------------------------
      SUBROUTINE hornerd_mf(v,XVAR,ORDER,t)     
        IMPLICIT NONE
        INTEGER NVAR,NPAR
        INTEGER nitermax,nordinc,minord,maxord
        COMMON /VP/ NVAR, NPAR
        COMMON /CONSTMET2/ nitermax,nordinc,minord,maxord
        INTEGER ORDER, i, nd
        REAL*8 XVAR(0:maxord,0:NVAR), t, au
        REAL*8 v(NVAR)
        DO nd=1, NVAR
          au = ORDER*XVAR(ORDER,nd)*t
          DO i = ORDER-1, 2, -1
            au = (au+i*XVAR(i,nd))*t
          END DO
          v(nd) = au+XVAR(1,nd)
        END DO
        RETURN
       END SUBROUTINE

C--------------------------------------------------------------
C--------------------------------------------------------------
C     mul_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION mul_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 mul_mf
        au = XX(0,n1)*XX(i,n2)
        DO m = 1, i
          au = au + XX(m,n1)*XX(i-m,n2)
        END DO
        mul_mf = au
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     div_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION div_mf(n1,n2,n3,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,n3,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 div_mf
        IF(XX(0,n2) .EQ. 0.d0)THEN
          WRITE(*,*) ' Function div_mf found division by zero'
          STOP
        ELSE
          au = XX(i,n1)
          DO m = 1, i
            au = au - XX(m,n2)*XX(i-m,n3)
          END DO
          div_mf = au/XX(0,n2)
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     inv_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION inv_mf(p,n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au,p
        REAL*8 inv_mf
        inv_mf = 0.d0
        IF(i.EQ.0)THEN
          IF(XX(0,n1) .EQ. 0.d0) THEN
            WRITE(*,*) ' Function inv_mf found division by zero'
            STOP
          ELSE
            inv_mf = p/XX(0,n1)
          END IF
        ELSE
          au = 0.d0
          DO m = 0, i -1
            au = au - (XX(i-m,n1)*XX(m,n2))
          END DO
          inv_mf = au*p/XX(0,n1)
        END IF
        RETURN
      END


C--------------------------------------------------------------
C--------------------------------------------------------------
C     pow_mf_c
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION pow_mf_c(n1,ex,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 ex,au
        REAL*8 pow_mf_c
        IF(i.EQ.0)THEN
          IF(XX(0,n1) .EQ. 0.d0) THEN
            WRITE(*,*) ' Function pow_mf_c found division by zero'
            STOP
          ELSE
            pow_mf_c = XX(0,n1)**ex
          END IF
        ELSE
          IF(XX(0,n1) .NE. 0.d0) THEN
            au = ex*i*XX(0,n2)*XX(i,n1)
            DO m = 1, i - 1
              au = au+(ex*(i-m)-m)*XX(m,n2)*XX(i-m,n1)
            END DO
            pow_mf_c = au/(i*XX(0,n1))
          ELSE
            pow_mf_c = 0.D0
          END IF
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     exp_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION exp_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 exp_mf
        IF(i.EQ.0)THEN
          exp_mf = EXP(XX(0,n1))
        ELSE
          au = i*XX(0,n2)*XX(i,n1)
          DO m = 1, i - 1
            au = au+(i-m)*XX(m,n2)*XX(i-m,n1)
          END DO
          exp_mf = au/i
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     log_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION log_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL, maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 log_mf
        IF(i.EQ.0)THEN
          IF(XX(0,n1) .LE. 0.d0)THEN
            WRITE(*,*) 'Function log_mf found log of a negative value'
            STOP
          ELSE
            log_mf = LOG(XX(0,n1))
          END IF
        ELSE
          au = i*XX(i,n1)
          DO m = 1, i - 1
            au = au-(i-m)*XX(m,n1)*XX(i-m,n2)
          END DO
          log_mf = au/i/XX(0,n1)
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     sin_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION sin_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 sin_mf
        IF(i.EQ.0)THEN
          sin_mf = SIN(XX(0,n1))
        ELSE
          au = XX(1,n1)*XX(i-1,n2)
          DO m = 2, i
            au = au+m*XX(m,n1)*XX(i-m,n2)
          END DO
          sin_mf = au/i
        END IF
        RETURN
      END

C--------------------------------------------------------------
C--------------------------------------------------------------
C     cos_mf
C--------------------------------------------------------------
C--------------------------------------------------------------
      FUNCTION cos_mf(n1,n2,i,XX,TOTAL,maxord)
        IMPLICIT NONE
        INTEGER TOTAL,maxord     
        INTEGER n1,n2,i,m
        REAL*8 XX(0:maxord,0:TOTAL)
        REAL*8 au
        REAL*8 cos_mf
        IF(i.EQ.0)THEN
          cos_mf = COS(XX(0,n1))
        ELSE
          au = XX(1,n1)*XX(i-1,n2)
          DO m = 2, i
            au = au+m*XX(m,n1)*XX(i-m,n2)
          END DO
          cos_mf = -au/i
        END IF
        RETURN
      END


"


(* ::Subsection::Closed:: *)
(*Header Standard DP*)


headstdDP = 
"
/****************************************************************************
 libTIDES. 
 This file is part of TIDES.
 
 Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>


#ifndef Header_DP_TIDES_h
#define Header_DP_TIDES_h


//Extern

extern int		_info_steps_taylor_;
extern int		num_etapas;
extern int		order_series;
extern int		order_estimator;
extern double   fac1;
extern double   fac2;
extern double   fac3;
extern double   rmaxstep; 
extern double   rminstep; 
extern int      nitermax; 
extern int      nordinc; 
extern int      minord; 
extern int      defect_error_control;
extern int      kahan_summation;
extern int      compensated_horner;

extern int		test_relminstep ;
extern double   dp_relminstep;


extern double	etapa_dp_minima;
extern double	etapa_dp_maxima;
extern double	etapa_dp_total;

extern int		printError ;

//Defines
#define initialize_dp_case() \\
long  NUM_COLUMNS;\\
set_links(itd, LINKS, POS_FUNCTIONS);\\
NUM_COLUMNS = (itd->NVARS+itd->NFUNS)*itd->NDER;\\
if(ORDER < 0) return NUM_COLUMNS;\\
check_iteration_data_parameters(0, itd->NVARS, VARIABLES);\\
check_iteration_data_parameters(1, itd->NPARS, PARAMETERS);\\
check_iteration_data_parameters(2, itd->NFUNS, FUNCTIONS);\\
itd->LINK_ELEMENTS = itd->NDER*itd->MAX_ORDER;\\
double var[itd->NVARS+1][itd->LINK_ELEMENTS];\\
double par[itd->NPARS][itd->LINK_ELEMENTS];\\
double link[LINKS][itd->LINK_ELEMENTS];\\
varDB_init(itd,(double*)var,v,t);\\
parDB_init(itd,(double*)par,p);\\
linkDB_init(itd,(double*)link);\\
derDB_init(itd,(double*)var, (double*)par, v);

#define write_dp_solution()	write_sol_DB(itd,cvfd,(double*)var,(double*)link);

//typedef
typedef struct it_dt{
	int		MAX_ORDER;
	long	NDER;
	int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
	int		*PARTIAL_LIST, *FUNCTION_LIST;
	long	*PREV_ACCUM, *PREV_VI, *PREV_IV, *PREV_COEF; 
	long	*PREVSTAR_ACCUM, *PREVSTAR_VI, *PREVSTAR_IV, *PREVSTAR_COEF;
	long	LINK_ELEMENTS;
	int		clearPartials;
} iteration_data;

typedef long (*DBLinkedFunction)(iteration_data *itd, double t, double *v, double *p, int orden, double *cvfd);

typedef struct dp_DM {
	int rows;
	int columns;
	double **data;
} dp_data_matrix;


typedef struct dp_node_event {
	double *data;
	struct dp_node_event *next;
} dp_event_node;

typedef struct dp_linked_events_nodes {
	int total;
	int dim;
	dp_event_node *first;
} dp_event_list;

typedef  double*	Array1DB;
typedef  double**	Array2DB;

//CommonITER
int		same_der(int dim, int *der1, int *der2);
void	string_to_der(int dim, char *sder, int *der);
long	position_derivative(char* sder, int *pdd);
long	position_variable(int v, char* der, int NVARS, int NFUNS, int *pdd);
long	position_function(int f, char* der, int NVARS, int NFUNS, int *pdd);
int		is_variable(iteration_data *itd, int num);
void	set_iteration_parameters(iteration_data *itd, int v, int p, int f, int o, int *pdd);
void	delete_iteration_parameters(iteration_data *itd);
void	set_links(iteration_data *itd, int l, int *flst);
void	check_iteration_data_parameters(int t, int pa, int pb);
void	set_info_error(void);
void	unset_info_error(void);
void    set_maxnumsteps(unsigned long val);

//doubITER
void	varDB_init(iteration_data *itd, double *var,  double *v, double t);
void	parDB_init(iteration_data *itd, double *par, double *p);
void	linkDB_init(iteration_data *itd, double *lk);
void	derDB_init(iteration_data *itd, double *var, double *par, double *v);
void	write_sol_DB(iteration_data *itd, double *cvfd, double *var, double *link);
void	double_htilde(iteration_data *itd, double *h, long j, long v, long i, double *ht, int ORDER_INDEX);
void	double_var_t(iteration_data *itd, double *f, double *u, int ORDER_INDEX);
void	double_var_t_c(iteration_data *itd, char* cs, double *w, int ORDER_INDEX);
void	double_var_t_cc(iteration_data *itd, double c, double *w, int ORDER_INDEX);
void	double_add_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_sub_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_add_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_sub_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_add_t_cc(iteration_data *itd, double c, double *u,  double *w, int ORDER_INDEX);
void	double_sub_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);
void	double_mul_t(iteration_data *itd, double *u, double *v, double *w, int ORDER_INDEX);
void	double_mul_t_c(iteration_data *itd, char* cs, double *u, double *w, int ORDER_INDEX);
void	double_mul_t_cc(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);
void	double_div_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_div_t_vc(iteration_data *itd, double *u, double c, double *w, int ORDER_INDEX);
void	double_div_t_cv(iteration_data *itd, double c, double *u, double *w, int ORDER_INDEX);
void	double_inv_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_exp_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_pow_t_c(iteration_data *itd, double *u, char* cs, double *w, int ORDER_INDEX);
void	double_pow_t_cc(iteration_data *itd, double *u, double c, double *w, int ORDER_INDEX);
void	double_sct_0(iteration_data *itd, double *f, double *g, double *h, long i);
void	double_sct_i(iteration_data *itd, double *f, double *g, double *h, long i, int ORDER_INDEX);
void	double_sin_t(iteration_data *itd, double *s, double *c, double *w, int ORDER_INDEX);
void	double_cos_t(iteration_data *itd, double *c, double *s, double *w, int ORDER_INDEX);
void    double_sin_cos_t (iteration_data *itd, double *f, double *s, double *c,  int ORDER_INDEX);
void	double_sinh_t(iteration_data *itd, double *s, double *c, double *w, int ORDER_INDEX);
void	double_cosh_t(iteration_data *itd, double *c, double *s, double *w, int ORDER_INDEX);
void    double_sinh_cosh_t (iteration_data *itd, double *f, double *s, double *c, int ORDER_INDEX);
void	double_fgt_0(iteration_data *itd, double *f, double *g, double *h, long i);
void	double_fgt_i(iteration_data *itd, double *f, double *g, double *h, long i, int ORDER_INDEX);
void	double_log_t(iteration_data *itd, double *u, double *w, int ORDER_INDEX);
void	double_asin_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_acos_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_atan_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_asinh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_acosh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);
void	double_atanh_t(iteration_data *itd, double *f, double *g, double *h, int ORDER_INDEX);

//doubNUM
void	double_init(double *rop); 
void	double_set_d(double *rop, double op); 
void	double_set_str(double *rop, char *op); 
void	double_set(double *rop, double op); 
double	double_get_d(double op); 
void	double_clear (double op);
void	double_add(double *rop, double op1, double op2); 
void	double_sub(double *rop, double op1, double op2); 
void	double_mul(double *rop, double op1, double op2); 
void	double_div(double *rop, double op1, double op2); 
void	double_pow(double *rop, double op1, double op2); 
void	double_abs(double *rop, double op);  
void	double_add_i (double *rop, double op1, long   op2); 
void	double_sub_i (double *rop, double op1, long   op2); 
void	double_i_sub (double *rop, long   op1, double op2);
void	double_mul_i (double *rop, double op1, long   op2); 
void	double_div_i (double *rop, double op1, long  op2); 
void 	double_i_div (double *rop, long   op1, double op2); 
void 	double_pow_i (double *rop, double op1, long   op2);
void	double_i_pow (double *rop, unsigned long  op1, double op2);
int		double_greater(double op1, double op2); 
int		double_greaterequal(double op1, double op2);
int		double_less(double op1, double op2); 
int		double_lessequal(double op1, double op2);
int		double_equal(double op1, double op2); 
void	double_log(double *rop, double op); 
void	double_log10(double *rop, double op); 
void	double_exp(double *rop, double op); 
void	double_exp2(double *rop, double op); 
void	double_exp10(double *rop, double op); 
void	double_cos(double *rop, double op); 
void	double_sin(double *rop, double op);
void	double_sin_cos(double *rsin, double *rcos, double op); 
void	double_tan(double *rop, double op);
void	double_sec(double *rop, double op);
void	double_csc(double *rop, double op);
void	double_cot(double *rop, double op);
void	double_acos(double *rop, double op);
void	double_asin(double *rop, double op);
void	double_atan(double *rop, double op);
void	double_atan2(double *rop, double op1, double op2); 
void	double_cosh(double *rop, double op);
void	double_sinh(double *rop, double op);
void	double_tanh(double *rop, double op);
void	double_sech(double *rop, double op);
void	double_csch(double *rop, double op);
void	double_coth(double *rop, double op);
void	double_acosh(double *rop, double op);
void	double_asinh(double *rop, double op);
void	double_atanh(double *rop, double op);
void	Array1DB_init(Array1DB *vec, long dim); 
void	Array1DB_clear(Array1DB *vec, long dim);
void	Array2DB_init(Array2DB *vec, long rows, long columns); 
void	Array2DB_clear(Array2DB *vec, long rows, long columns);
void	Array3DB_init(Array2DB *vec, long dim, long rows, long columns); 
void	Array1DB_set(Array1DB rop, Array1DB op, long dim);
void	Array2DB_set(Array2DB rop, Array2DB op, long rows, long columns);
void	Array2DB_column_set(Array2DB rop, Array1DB op, long c, long dim);
void	Array2DB_row_set(Array2DB rop, Array1DB op, long r, long dim);
void	double_write (char *c, double op);

//doubPOL
void dp_pol_derivative(double *pol, int grado, double *dpol);
void dp_horner(double *pol, int grado, double t, double *eval);
void dp_horner_der(double *pol, int grado, double t, double *eval);
void dp_classic_horner(double *pol, int grado, double t, double *eval);
void dp_twosum(double a, double b, double *x, double *y);
void dp_split_real(double a, double *ah, double *al);
void dp_twoproduct(double a, double b, double *x, double *y);
void dp_compensated_horner(double *pol, int grado, double t, double *eval);
void dp_ch_evaluate_poly_der(double *pol, int grado, double t, double *vpol, double *vder);
int  dp_one_zero(double *pol, int grado, double ta, double tb, double tol, double *raiz);
int  dp_one_extremum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_one_maximum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_one_minimum(double *pol, int grado, double ta, double tb, double tol, double *extr);
int  dp_brent(double *pol, int grado, double ta, double tb, double tol, double *raiz);
int  dp_rtsafe(double *pol, int grado, double x1, double x2, double tol, double *raiz);

//doubTODE
void	check_relminstep(void);
void	uncheck_relminstep(void);
void	dp_set_relminstep(double val);

void use_default_step_estimator(void);
void dp_set_info_taylor(void);
void dp_unset_info_taylor(void);
void dp_str_info_taylor(void);
void dp_add_info_step(double tstep);
int  dp_taylor_order(double eps);
void dp_norm_inf_var(double *rop, int nvar, double *y) ;
void dp_norm_inf(double *rop, int n, int k, double *coef, int MAX_ORDER);
void dp_compute_step(double *rop, double *hant, double tol, int n, int ord, double *cvfd, int MAX_ORDER);
void dp_compute_tol (double *tol, double tolrel, double tolabs, int nvar, double *y) ;
void dp_taylor_horner(int n, int ord, double *cvfd, double t, double *x, int MAX_ORDER);
void dp_taylor_horner_der (int n, int ord, double *cvfd, double t, double *x, int MAX_ORDER);
void dp_write_taylor_solution(int n,  int j, double tini, double *x, double ***mat, FILE* fileout);
int  dp_valid_step (DBLinkedFunction fcn,  iteration_data *itd, double *step, double tip,double tol, int nvar, int ncol, int order, double *cvfd,double *p, int MAX_ORDER);
int  dp_tides(DBLinkedFunction fcn,int *pdd, int nvar,  int npar, int nfun, double *x, double *p, double *lt, int ntes, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout);
int  dp_tides_point(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double t0, double tf, double dt, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout);
int  dp_tides_list(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double *lt, int ntes, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout) ;
int  dp_tides_delta(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double t0, double dt, int ntot, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout);
int  dp_tides_delta_ft(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double t0, double tf, double dt, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout);
int  dp_tides_kernel(DBLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, double *x, double *p, double t0, double dt, double *dtl, int ntot, double tolrel, double tolabs, dp_data_matrix *mat, FILE* fileout, dp_data_matrix *der) ;
long dp_number_of_columns(DBLinkedFunction fcn);
void init_dp_data_matrix(dp_data_matrix *dm, int r, int c);
void delete_dp_data_matrix(dp_data_matrix *dm);

//doubEVENTS
void dp_init_event_list(dp_event_list *lista, int ncol);
void dp_add_to_event_list(dp_event_list *lista, double *val);
void dp_event_list_to_array(dp_event_list lista, dp_data_matrix *array);
int  dp_tides_find_zeros(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout) ;
int  dp_tides_find_extrema(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout) ;
int  dp_tides_find_minimum(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout) ;
int  dp_tides_find_maximum(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout);
int  dp_tides_events(DBLinkedFunction fcn, int nvar, int npar,double *x, double *p, double tini, double tend, double tol, int *numevents, dp_data_matrix *events, FILE* fileout, int evcase) ;


#endif





";


(* ::Subsection::Closed:: *)
(*Header Standard MP*)


headstdMP = 
"
/****************************************************************************
 libTIDES. 
 This file is part of TIDES.
 
 Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez
 Grupo de Mecanica Espacial
 University of Zaragoza
 SPAIN
 
 http://gme.unizar.es/software/tides
 Contact: <tides@unizar.es>
 
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include \"mpfr.h\"

#ifndef Header_DP_TIDES_h
#define Header_DP_TIDES_h

//Extern

extern int		_info_steps_taylor_;
extern int		num_etapas;
extern int		order_series;
extern int		order_estimator;
extern double   fac1;
extern double   fac2;
extern double   fac3;
extern double   rmaxstep; 
extern double   rminstep; 
extern int      nitermax; 
extern int      nordinc; 
extern int      minord; 
extern int      defect_error_control;
extern int      kahan_summation;
extern int      compensated_horner;

extern int		test_relminstep ;
extern int		printError ;

//typedef
typedef struct it_dt{
	int		MAX_ORDER;
	long	NDER;
	int		NVARS, NPARS, NFUNS, NLINKS, NPARTIALS; 
	int		*PARTIAL_LIST, *FUNCTION_LIST;
	long	*PREV_ACCUM, *PREV_VI, *PREV_IV, *PREV_COEF; 
	long	*PREVSTAR_ACCUM, *PREVSTAR_VI, *PREVSTAR_IV, *PREVSTAR_COEF;
	long	LINK_ELEMENTS;
	int		clearPartials;
} iteration_data;


//CommonITER
int		same_der(int dim, int *der1, int *der2);
void	string_to_der(int dim, char *sder, int *der);
long	position_derivative(char* sder, int *pdd);
long	position_variable(int v, char* der, int NVARS, int NFUNS, int *pdd);
long	position_function(int f, char* der, int NVARS, int NFUNS, int *pdd);
int		is_variable(iteration_data *itd, int num);
void	set_iteration_parameters(iteration_data *itd, int v, int p, int f, int o, int *pdd);
void	delete_iteration_parameters(iteration_data *itd);
void	set_links(iteration_data *itd, int l, int *flst);
void	check_iteration_data_parameters(int t, int pa, int pb);
void	set_info_error(void);
void	unset_info_error(void);

void	check_relminstep(void);
void	uncheck_relminstep(void);
void	use_default_step_estimator(void);
void    set_maxnumsteps(unsigned long val);

#endif



#ifndef Header_MP_TIDES_h
#define Header_MP_TIDES_h

//Extern

extern mpfr_t   mp_relminstep;

extern mpfr_t	etapa_mp_minima;
extern mpfr_t	etapa_mp_maxima;
extern mpfr_t	etapa_mp_total;

extern mpfr_rnd_t	TIDES_RND;
extern int		TIDES_PREC;
extern int		BINARY_PRECISION;
extern int		DECIMAL_PRECISION;


//Defines
#define initialize_mp_case() \\
long  NUM_COLUMNS;\\
set_links(itd, LINKS, POS_FUNCTIONS);\\
NUM_COLUMNS = (itd->NVARS+itd->NFUNS)*itd->NDER;\\
if(ORDER < 0) return NUM_COLUMNS;\\
check_iteration_data_parameters(0, itd->NVARS, VARIABLES);\\
check_iteration_data_parameters(1, itd->NPARS, PARAMETERS);\\
check_iteration_data_parameters(2, itd->NFUNS, FUNCTIONS);\\
itd->LINK_ELEMENTS = itd->NDER*itd->MAX_ORDER;\\
mpfr_t var[itd->NVARS+1][itd->LINK_ELEMENTS];\\
mpfr_t par[itd->NPARS][itd->LINK_ELEMENTS];\\
mpfr_t link[LINKS][itd->LINK_ELEMENTS];\\
varMP_init(itd,(mpfr_t*)var,v,t);\\
parMP_init(itd,(mpfr_t*)par,p);\\
linkMP_init(itd,(mpfr_t*)link);\\
derMP_init(itd,(mpfr_t*)var, (mpfr_t*)par, v);

#define write_mp_solution()	write_sol_MP(itd,cvfd,(mpfr_t*)var,(mpfr_t*)link);

#define clear_vpl()	\\
clear_mpfr_vpl(itd, (mpfr_t*)var, (mpfr_t*)par, (mpfr_t*)link);

#define clear_cts()	\\
for(i = 0; i < NCONST ; i++ ) mpfr_clear(ct[i]);

//typedef

typedef long (*MPLinkedFunction)(iteration_data *itd, mpfr_t t, mpfr_t *v, mpfr_t *p, int orden, mpfr_t *coef);

typedef struct mp_DM {
	int rows;
	int columns;
	mpfr_t **data;
} mp_data_matrix;


typedef struct mp_node_event {
	mpfr_t *data;
	struct mp_node_event *next;
} mp_event_node;

typedef struct mp_linked_events_nodes {
	int total;
	int dim;
	mp_event_node *first;
} mp_event_list;

typedef  mpfr_t*	Array1MP;
typedef  mpfr_t**	Array2MP;

//doubNUMdef
#define TIDES_RND GMP_RNDN

#define  set_precision_digits	mpfrts_set_prec
#define  set_rounding_mode  	mpfrts_set_rnd
	


//mpfrITER		
void	varMP_init(iteration_data *itd, mpfr_t *var, mpfr_t *v, mpfr_t t);
void	parMP_init(iteration_data *itd, mpfr_t *par, mpfr_t *p);
void	linkMP_init(iteration_data *itd, mpfr_t *lk);
void	derMP_init(iteration_data *itd, mpfr_t *var, mpfr_t *par, mpfr_t *v);
void 	clear_mpfr_vpl (iteration_data *itd, mpfr_t *var, mpfr_t *par, mpfr_t *link);
void	write_sol_MP(iteration_data *itd, mpfr_t *cvfd, mpfr_t *var, mpfr_t *link);
void	mpfrts_htilde(iteration_data *itd, mpfr_t *h, long j, long v, long i, mpfr_t *ht, int ORDER_INDEX);
void	mpfrts_var_t(iteration_data *itd, mpfr_t *f, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_var_t_c(iteration_data *itd, char* cs, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_var_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_add_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_add_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_add_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sub_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_mul_t(iteration_data *itd, mpfr_t *u, mpfr_t *v, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_mul_t_c(iteration_data *itd, char* cs, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_mul_t_cc(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_div_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_div_t_vc(iteration_data *itd, mpfr_t *u, mpfr_t c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_div_t_cv(iteration_data *itd, mpfr_t c, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_inv_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_exp_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_pow_t_c(iteration_data *itd, mpfr_t *u, char* cs, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_pow_t_cc(iteration_data *itd, mpfr_t *u, mpfr_t c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_sct_0(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i);
void	mpfrts_sct_i(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i, int ORDER_INDEX);
void	mpfrts_sin_t(iteration_data *itd, mpfr_t *s, mpfr_t *c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_cos_t(iteration_data *itd, mpfr_t *c, mpfr_t *s, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_sin_cos_t (iteration_data *itd, mpfr_t *f, mpfr_t *s, mpfr_t *c, int ORDER_INDEX);
void	mpfrts_sinh_t(iteration_data *itd, mpfr_t *s, mpfr_t *c, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_cosh_t(iteration_data *itd, mpfr_t *c, mpfr_t *s, mpfr_t *w, int ORDER_INDEX);
void    mpfrts_sinh_cosh_t (iteration_data *itd, mpfr_t *f,  mpfr_t *s, mpfr_t *c, int ORDER_INDEX);
void	mpfrts_fgt_0(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i);
void	mpfrts_fgt_i(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, long i, int ORDER_INDEX);
void	mpfrts_log_t(iteration_data *itd, mpfr_t *u, mpfr_t *w, int ORDER_INDEX);
void	mpfrts_asin_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_acos_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_atan_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_asinh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_acosh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);
void	mpfrts_atanh_t(iteration_data *itd, mpfr_t *f, mpfr_t *g, mpfr_t *h, int ORDER_INDEX);



//mpfrNUM
int		binary_precision(int prec) ;
void	mpfrts_init (mpfr_t *rop); 
void	mpfrts_set_i(mpfr_t  *rop, long op); 
void	mpfrts_set_d (mpfr_t *rop, double op); 
void	mpfrts_set_str (mpfr_t *rop, char *op); 
void	mpfrts_set (mpfr_t *rop, mpfr_t op); 
double	mpfrts_get_d (mpfr_t op); 
long    mpfrts_get_i(mpfr_t op);
int		mpfrts_get_prec (void); 
void 	mpfrts_set_prec (int dig);
void	mpfrts_add (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_add_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_sub (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void 	mpfrts_sub_i (mpfr_t *rop, mpfr_t op1, long int op2); 
void 	mpfrts_i_sub (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_mul (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_mul_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_div (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void 	mpfrts_div_i (mpfr_t *rop, mpfr_t op1, long int op2);
void	mpfrts_i_div (mpfr_t *rop, long int op1, mpfr_t op2);
void	mpfrts_pow (mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_pow_i (mpfr_t *rop, mpfr_t op1, long int op2);
void 	mpfrts_i_pow (mpfr_t *rop, unsigned long int op1, mpfr_t op2);
void	mpfrts_abs(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_neg(mpfr_t  *rop, mpfr_t op);
int		mpfrts_greater(mpfr_t op1, mpfr_t op2); 
int		mpfrts_greaterequal(mpfr_t op1, mpfr_t op2);
int		mpfrts_less(mpfr_t op1, mpfr_t op2); 
int		mpfrts_lessequal(mpfr_t op1, mpfr_t op2);
int		mpfrts_equal(mpfr_t op1, mpfr_t op2); 
void	mpfrts_sqrt(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log2(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_log10(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp2(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_exp10(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_cos(mpfr_t  *rop, mpfr_t op); 
void	mpfrts_sin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sin_cos(mpfr_t *rsin, mpfr_t *rcos, mpfr_t op); 
void	mpfrts_tan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sec(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csc(mpfr_t  *rop, mpfr_t op);
void	mpfrts_cot(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acos(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asin(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atan2(mpfr_t  *rop, mpfr_t op1, mpfr_t op2); 
void	mpfrts_cosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_tanh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_sech(mpfr_t  *rop, mpfr_t op);
void	mpfrts_csch(mpfr_t  *rop, mpfr_t op);
void	mpfrts_coth(mpfr_t  *rop, mpfr_t op);
void	mpfrts_acosh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_asinh(mpfr_t  *rop, mpfr_t op);
void	mpfrts_atanh(mpfr_t  *rop, mpfr_t op);
void 	mpfrts_write_var(mpfr_t op);
void 	mpfrts_write (char *c, mpfr_t op);
void 	mpfrts_fread (FILE *file, mpfr_t rop);
void 	mpfrts_fwrite (FILE *file, mpfr_t op, int prec); 
void	Array1MP_init(Array1MP *vec, long dim);
void	Array1MP_clear(Array1MP *vec, long dim);
void	Array2MP_init(Array2MP *vec, long rows, long columns);
void	Array2MP_clear(Array2MP *vec, long rows, long columns);
void	Array1MP_set(Array1MP rop, Array1MP op, long dim);
void	Array2MP_set(Array2MP rop, Array2MP op, long rows, long columns);
void	Array2MP_column_set(Array2MP rop, Array1MP op, long c, long dim);
void	Array2MP_row_set(Array2MP rop, Array1MP op, long r, long dim);

//mpfrPOL
void mp_pol_derivative(mpfr_t *pol, int grado, mpfr_t *dpol);
void mp_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_horner_der(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_classic_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_twosum(mpfr_t a, mpfr_t b, mpfr_t *x, mpfr_t *y);
void mp_split_real(mpfr_t a, mpfr_t *ah, mpfr_t *al);
void mp_twoproduct(mpfr_t a, mpfr_t b, mpfr_t *x, mpfr_t *y);
void mp_compensated_horner(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *eval);
void mp_ch_evaluate_poly_der(mpfr_t *pol, int grado, mpfr_t t, mpfr_t *vpol, mpfr_t *vder);
int  mp_one_zero(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *raiz);
int  mp_one_extremum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_one_maximum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_one_minimum(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *extr);
int  mp_brent(mpfr_t *pol, int grado, mpfr_t ta, mpfr_t tb, mpfr_t tol, mpfr_t *raiz);
int  mp_rtsafe(mpfr_t *pol, int grado, mpfr_t x1, mpfr_t x2, mpfr_t tol, mpfr_t *raiz);

//mpfrTODE
void mp_set_relminstep(mpfr_t val);
void mp_set_info_taylor(void);
void mp_unset_info_taylor(void);
void mp_str_info_taylor(void);
void mp_add_info_step(mpfr_t tstep);
int  mp_taylor_order(mpfr_t eps);
void mp_norm_inf_vec(mpfr_t *rop, int nvar, mpfr_t *y) ;
void mp_norm_inf(mpfr_t *rop, int n, int k, mpfr_t *coef, int MAX_ORDER);
void mp_compute_step(mpfr_t *rop, double *hant, mpfr_t tol, int n, int ord, mpfr_t *coef, int MAX_ORDER);
void mp_compute_tol (mpfr_t *tol, mpfr_t tolrel, mpfr_t tolabs, int nvar, mpfr_t *y); 
void mp_taylor_horner(int n, int ord, mpfr_t *coef, mpfr_t t, mpfr_t *x, int MAX_ORDER);
void mp_taylor_horner_der(int n, int ord, mpfr_t *coef, mpfr_t t, mpfr_t *x, int MAX_ORDER);
void mp_write_taylor_solution( int n, int j, mpfr_t tini, mpfr_t *x,  mpfr_t*** mat, FILE* fileout);	
int  mp_valid_step (MPLinkedFunction fcn, iteration_data *itd, mpfr_t *step, mpfr_t tip, mpfr_t eps, int nvar, int ncol, int order, mpfr_t *cvfd, mpfr_t *p, int MAX_ORDER);
int  mp_tides(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t *lt, int ntes, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout);
int  mp_tides_point(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t t0, mpfr_t tf, mpfr_t dt, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout);
int  mp_tides_list(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t *lt, int ntes, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout) ;
int  mp_tides_delta(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t t0, mpfr_t dt, int ntot, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout);
int  mp_tides_delta_ft(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t t0, mpfr_t tf, mpfr_t dt, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout);
int  mp_tides_kernel(MPLinkedFunction fcn, int *pdd, int nvar, int npar, int nfun, mpfr_t *x, mpfr_t *p, mpfr_t t0, mpfr_t dt,  mpfr_t *dtl, int ntot, mpfr_t tolrel, mpfr_t tolabs, mp_data_matrix *mat, FILE* fileout, mp_data_matrix *der);
long mp_number_of_columns(MPLinkedFunction fcn);
void init_mp_data_matrix(mp_data_matrix *dm, int r, int c);
void delete_mp_data_matrix(mp_data_matrix *dm);
int  getOrder (void);
int  getNsteps (void);

//mpfrEVENTS
void mp_init_event_list(mp_event_list *lista, int ncol);
void mp_add_to_event_list(mp_event_list *lista, mpfr_t *val);
void mp_event_list_to_array(mp_event_list lista, mp_data_matrix *array);

int  mp_tides_find_zeros(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout) ;
int  mp_tides_find_extrema(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout);
int  mp_tides_find_minimum(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout);
int  mp_tides_find_maximum(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout);
int  mp_tides_events(MPLinkedFunction fcn, int nvar, int npar,mpfr_t *x, mpfr_t *p, mpfr_t tini, mpfr_t tend, mpfr_t tol, int *numevents, mp_data_matrix *events, FILE* fileout,int evcase);

#endif


";


(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`Texts`"]

EndPackage[]

Null
