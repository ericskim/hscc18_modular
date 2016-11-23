(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`PWSTides: Symbolic Taylor Series Method and Power Series*)


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
(*PWSTides*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`PWSTides`",
				"MathTIDES`ODES`",
					"MathTIDES`LKFunctions`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
	PWS, 
	ZeroPWS,
	UnitPWS,
	OrderPWS,
	PWSeries,
	TSMSolve,
	TSMSolve$
}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`PWSTides`*"]
Clear @@ Names["MathTIDES`PWSTides`*"]


(* ::Section::Closed:: *)
(*Mensajes*)


Begin["`mess`"]


PWS::usage = "..." 
ZeroPWS::usage = "..."
UnitPWS::usage = "..."
OrderPWS::usage = "..."
PWSeries::usage = "..."
TSMSolve::usage = "..."
PWS::"divbyzero"="Order zero of divisor is equal to zero"

TSMSolve::"badvarg"="bad variables arguments"
TSMSolve::"badparg"="bad parameter arguments"



End[]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsection::Closed:: *)
(*Series de Taylor : PWS*)


ZeroPWS[n_]:= PWS @@ Table[0,{n+1}]

UnitPWS[n_]:= 
	Module[{ts},
		ts = ZeroPWS[n];
		ts[[1]] = 1;
		ts]

OrderPWS[t_PWS]:= Length[t]-1
OrderPWS[t_PWS, n_]:= PadRight[t,n+1]


PWSeries[ser_, {eps_Symbol, zero_,  n_Integer}]:=
	Module[{ls,ser2},
		ser2 = Normal[Series[ser,{eps,zero,n}]];
		ls = Map[(Coefficient[ser2,eps-zero,#] )&,Range[0,n]];
		PWS @@ Map[Simplify,ls]]
	
PWSeries[ser_List, {eps_Symbol, zero_,  n_Integer}]:=
	Map[PWSeries[#,{eps,zero,n}]&, ser]


horner[eps_,a_PWS]:=
	Module[{na, total},
		na = Reverse[List @@ a];
		total = na[[1]];
		na = Drop[na,1];
		Map[(total = (total * eps + #))&, na];
		total
	]


PWS[c__][t_?NumberQ]:= Expand[horner[t,PWS[c]]]
PWS[c__][t_]:= Collect[horner[t,PWS[c]], t]


PWS/: D[ser_PWS, x__] := Map[D[#,x]&, ser]
PWS/: Dt[ser_PWS, x__]:= Map[Dt[#,x]&, ser]


(* ::Subsection::Closed:: *)
(*Aritm\[EAcute]tica de series*)


(* ::Subsubsection::Closed:: *)
(*Operaciones en cada orden *)


nullOrderZeroQ[u_PWS]:=
	If[u[[1]] === 0, 
		Message[PWS::"divbyzero"]; Abort[],
		False]

add$[u_PWS, v_PWS, w_PWS, k_Integer]:= u[[k+1]]+v[[k+1]]

add$[u_PWS, v_  , w_PWS, 0]:= u[[1]]+v 
add$[u_PWS, v_  , w_PWS, k_Integer]:= u[[k+1]] 
add$[u_  , v_PWS, w_PWS, 0]:= u   +v[[1]]
add$[u_  , v_PWS, w_PWS, k_Integer]:= v[[k+1]]

mult$[u_PWS, v_PWS, w_PWS,  k_Integer]:= Sum[ u[[j+1]] v[[k-j+1]], {j,0,k}]
mult$[u_PWS, v_, w_PWS,  k_Integer]  := u[[k+1]] v
mult$[u_, v_PWS, w_PWS,  k_Integer]  := u v[[k+1]]

div$[u_PWS, v_PWS, w_PWS,  0]:= (u[[1]]/v[[1]])/; !nullOrderZeroQ[v]
div$[u_PWS, v_PWS, w_PWS,  k_Integer]:= (u[[k+1]]/v[[1]] - 
		Sum[ Binomial[k,j] v[[j+1]] w[[k-j+1]], {j,1,k}]/v[[1]])
div$[u_PWS, v_, w_PWS,  k_Integer]  := u[[k+1]] / v
div$[1,u__]:=inverse$[u]
div$[u_, v_PWS, w_PWS,  k_Integer]  := u inverse$[v,w,k]


inverse$[u_PWS, w_PWS,  0]:= (1/u[[1]])/; !nullOrderZeroQ[u]
inverse$[u_PWS, w_PWS,  k_Integer]:= - Sum[ u[[j+1]] w[[k-j+1]], {j,1,k}]/u[[1]]

sck$[u_PWS, v_PWS, k_Integer, signo_]:= 
	signo Sum[ (k-j) v[[j+1]] u[[k-j+1]], {j,0,k-1}]/k

sin$[u_PWS, c_PWS, w_PWS,  0] := Sin[u[[1]]];
sin$[u_PWS, c_PWS, w_PWS,  k_Integer]:= sck$[u, c, k, 1]

cos$[u_PWS, s_PWS, w_PWS, 0] := Cos[u[[1]]];
cos$[u_PWS, s_PWS, w_PWS, k_Integer]:= sck$[u, s, k, -1]

sinh$[u_PWS, c_PWS, w_PWS, 0] := Sinh[u[[1]]];
sinh$[u_PWS, c_PWS, w_PWS, k_Integer]:= sck$[u, c, k, 1]

cosh$[u_PWS, s_PWS, w_PWS, 0] := Cosh[u[[1]]];
cosh$[u_PWS, s_PWS, w_PWS, k_Integer]:= sck$[u, s, k, 1]


compl$[u_PWS, g_PWS, w_PWS, k_]:=
	((u[[k+1]]- (Sum[ j w[[j+1]] g[[k-j+1]], {j,1,k-1}]/k))/g[[1]])/; 
			!nullOrderZeroQ[g]


log$[u_PWS, w_PWS, 0]:= Log[u[[1]]]
log$[u_PWS, w_PWS, k_Integer]:= compl$[u, u, w, k]

asin$[u_PWS, g_PWS, w_PWS, 0]:= ArcSin[u[[1]]]
asin$[u_PWS, g_PWS, w_PWS, k_Integer]:= compl$[u, g, w, k]

acos$[u_PWS, g_PWS, w_PWS, 0]:= ArcCos[u[[1]]]
acos$[u_PWS, g_PWS, w_PWS, k_Integer]:= compl$[u, g, w, k]

atan$[u_PWS, g_PWS, w_PWS, 0]:= ArcTan[u[[1]]]
atan$[u_PWS, g_PWS, w_PWS, k_Integer]:= compl$[u, g, w, k]

asinh$[u_PWS, g_PWS, w_PWS, 0]:= ArcSinh[u[[1]]]
asinh$[u_PWS, g_PWS, w_PWS, k_Integer]:= compl$[u, g, w, k]

acosh$[u_PWS, g_PWS, w_PWS, 0]:= ArcCosh[u[[1]]]
acosh$[u_PWS, g_PWS, w_PWS, k_Integer]:= compl$[u, g, w, k]

atanh$[u_PWS, g_PWS, w_PWS, 0]:= ArcTanh[u[[1]]]
atanh$[u_PWS, g_PWS, w_PWS, k_Integer]:= compl$[u, g, w, k]

pow$[u_PWS, v_PWS, s_PWS, t_PWS, w_PWS, 0]:= u[[1]]^v[[1]]
pow$[u_PWS, v_PWS, s_PWS, t_PWS, w_PWS, k_Integer]:= 
	Sum[
		Sum[(i+1) (u[[i+2]] s[[j-i+1]] + v[[i+2]] t[[j-i+1]]),
			{i,0,j}] w[[k-j]],
				{j,0,k-1}]/k
				
pow$[u_PWS, a_, w_PWS, 0]:= u[[1]]^a
pow$[u_PWS, a_, w_PWS, k_Integer]:= 
	(Sum[(a (k-j) -j) w[[j+1]] u[[k-j+1]],{j,0,k-1}]/k)/u[[1]]/; !nullOrderZeroQ[u]

pow$[a_, u_PWS,  w_PWS, 0]:= a^u[[1]]
pow$[a_, u_PWS,  w_PWS, k_Integer]:= 
	(Sum[(k-j) w[[j+1]] u[[k-j+1]], {j,0,k-1}]/k) Log[a]

exp$[u_PWS,  w_PWS, 0]:= E^u[[1]]
exp$[u_PWS,  w_PWS, k_Integer]:= 
	(Sum[(k-j) w[[j+1]] u[[k-j+1]], {j,0,k-1}]/k) 

variables$[u_PWS, w_PWS, k_Integer]:= u[[k]]/k

variables$[u_, w_PWS, 1]:= u
variables$[u_, w_PWS, k_Integer]:= 0

der$[u_PWS, w_PWS, k_Integer]:= (k+1) u[[k+2]]


rADheads = {
	LKF$Plus->add$,
	LKF$Times-> mult$,
	LKF$Divide->div$,
	LKF$Sin->sin$,
	LKF$Cos->cos$,
	LKF$Tan->tan$,
	LKF$Sinh->sinh$,
	LKF$Cosh->cosh$,
	LKF$Tanh->tanh$,
	LKF$ArcSin->asin$,
	LKF$ArcCos->acos$,
	LKF$ArcTan->atan$,
	LKF$ArcSinh->asinh$,
	LKF$ArcCosh->acosh$,
	LKF$ArcTanh->atanh$,
	LKF$Power->pow$,
	LKF$Exp->exp$,
	LKF$Log->log$,
	LKF$Der->der$
}


(* ::Subsubsection::Closed:: *)
(*Funciones generales*)


(*
PWSeries[fun_,{eps_Symbol, zero_, n_Integer},
	vv$:{Rule[_Symbol, _]...}]:=
  Module[{var, vl, adf, tsfun},
		{var, vl} = If[vv$ == {},{{},{}}, Transpose[vv$/.Rule->List]];
		If[vl =!= {}, vl = PWSeries[vl,{eps,zero,n}]];
		var = Prepend[var,eps];
		vl  = Prepend[vl,PWS @@ Join[{0,1},Table[0, {n-1}]]];
		adf = ToTaylorLKFPar[fun, var, {}];  
		tsfun = PWSExpansion[adf, vl, {}, n]; 
		tsfun
	]
	
PWSExpansion[adf_LKFPar, var_List, par_List, n_Integer]:=
	Module[{npar,lpar,rpar, tpar, liter,rvar, riter, niter,
				item, itemk , fun, link, ks, kss},
		npar = Length[adf[[3]]];
		rpar = Map[(LKF$Par[#]->par[[#]])&, Range[adf[[2]]]];
		lpar = Map[LKF$Par, Range[adf[[2]]+1,adf[[2]]+npar]];
		lpar = lpar //.Take[RightIterationLKF[adf],npar]/.fromLKF$Link/.rpar;
		tpar = par;
		Map[AppendTo[tpar,#]&,lpar];
		
		
		ChangeIndexDer[];
		Map[depthDer,IterationLKF[adf]];
		
		liter = adf[[4]]/.{LKF$Constant[x_]->x}; 
		niter = Length[liter]; 
		rpar = Map[(LKF$Par[#]->tpar[[#]])&, Range[Length[tpar]]];
		rvar = Map[(LKF$Var[#]->var[[#]])&, Range[Length[var]]];
		riter = {};
		link = Table[PWS[],{niter}]; 
		riter := Map[(LKF$Link[#]->link[[#]])&,Range[niter]]; 
		Do[ item = liter[[i]]/.rADheads; 
			Map[AppendTo[item,#]&,{LKF$Link[i],ks}];
			kss = k - depthDer[LKF$Link[i]];
			If[kss >= 0, 
				itemk = item/.rpar/.rvar/.riter/.ks->kss; 
				AppendTo[link[[i]],itemk]],{k,0,n}, {i,niter} 
		];
		fun = adf[[5]]/.riter/.rvar/.rpar/cpar;
		Clear[link];
		If[Length[fun]===1, fun = fun[[1]]];
		Map[Simplify, fun]
	] 
*)


(* ::Subsubsection::Closed:: *)
(*Operaciones completas*)


PWS/: a_PWS + b_PWS := 
	Module[{an,bn,nord},
		nord = Min[OrderPWS[a], OrderPWS[b]];
		an = OrderPWS[a,nord];
		bn = OrderPWS[b,nord];
		PWS @@ ((List @@ an) + (List @@ bn))
	]

PWS/: a_PWS + b_ := a/.{a[[1]]->a[[1]] + b}


PWS/: Times[a_PWS, b_PWS] := 
	Module[{prod, nord},
		nord = Min[OrderPWS[a], OrderPWS[b]];
		prod = ZeroPWS[nord];
		Do[prod[[k+1]] = mult$[a, b, prod,k], {k,0,nord}];
		Map[Simplify,prod]
	]

PWS/: Times[a_PWS, b_] := Map[(b #)&,a] 


PWS/: Power[a_PWS, b_PWS] := 
	PWSeries[Power[a[t],b[t]], {t, 0, Min[OrderPWS[a], OrderPWS[b]]}]

PWS/: Power[a_PWS, b_] := 
	PWSeries[Power[a[t],b], {t, 0, OrderPWS[a]}]

PWS/: Power[a_, b_PWS] := 
	PWSeries[Power[a,b[t]], {t, 0, OrderPWS[b]}]



PWS/: Sin[a_PWS]:=PWSeries[Sin[a[t]], {t, 0, OrderPWS[a]}]
PWS/: Cos[a_PWS]:=PWSeries[Cos[a[t]], {t, 0, OrderPWS[a]}]
PWS/: Sec[a_PWS]:=PWSeries[Sec[a[t]], {t, 0, OrderPWS[a]}]
PWS/: Csc[a_PWS]:=PWSeries[Csc[a[t]], {t, 0, OrderPWS[a]}]

PWS/: Sinh[a_PWS]:=PWSeries[Sinh[a[t]], {t, 0, OrderPWS[a]}]
PWS/: Cosh[a_PWS]:=PWSeries[Cosh[a[t]], {t, 0, OrderPWS[a]}]
PWS/: Sech[a_PWS]:=PWSeries[Sech[a[t]], {t, 0, OrderPWS[a]}]
PWS/: Csch[a_PWS]:=PWSeries[Csch[a[t]], {t, 0, OrderPWS[a]}]

PWS/: ArcSin[a_PWS]:=PWSeries[ArcSin[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcCos[a_PWS]:=PWSeries[ArcCos[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcTan[a_PWS]:=PWSeries[ArcTan[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcSec[a_PWS]:=PWSeries[ArcSec[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcCsc[a_PWS]:=PWSeries[ArcCsc[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcCot[a_PWS]:=PWSeries[ArcCot[a[t]], {t, 0, OrderPWS[a]}]

PWS/: ArcSinh[a_PWS]:=PWSeries[ArcSinh[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcCosh[a_PWS]:=PWSeries[ArcCosh[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcTanh[a_PWS]:=PWSeries[ArcTanh[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcSech[a_PWS]:=PWSeries[ArcSech[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcCsch[a_PWS]:=PWSeries[ArcCsch[a[t]], {t, 0, OrderPWS[a]}]
PWS/: ArcCoth[a_PWS]:=PWSeries[ArcCoth[a[t]], {t, 0, OrderPWS[a]}]

PWS/: Log[a_PWS]:=PWSeries[Log[a[t]], {t, 0, OrderPWS[a]}]


(* ::Subsection::Closed:: *)
(*Integraci\[OAcute]n simb\[OAcute]lica de ODEs*)


TSMSolve[fode_FirstOrderODE$, {t_Symbol, t0_, n_Integer}, vv$:{Rule[_Symbol, _]...}]:=
	TSMSolve[fode, {t,t0,n}, vv$, {}]

TSMSolve[fode_FirstOrderODE$, {t_Symbol, t0_, n_Integer},
	vv$:{Rule[_Symbol, _]...}, pp$:{Rule[_Symbol, _]...}]:=
  Module[{fun, fvar, iter, rules, vlp ={},
  			var, vl, vts, par, pl, tsfun},
		{var, vl} = If[vv$ == {},{{},{}}, Transpose[vv$/.Rule->List]];
		{par, pl} = If[pp$ == {},{{},{}}, Transpose[pp$/.Rule->List]];
		fun = ToTaylorLKFPar[fode];
		If[Length[var] =!= fun[[1]]-1,
			Message[TSMSolve::"badvarg"]; Print[var]; Abort[]];
		If[Length[par] =!= fun[[2]],
			Message[TSMSolve::"badparg"]; Print[par]; Abort[]];

		vlp = Prepend[vl,t0]; 
		TSMSolve$[fun, vlp, pl, n]
	]
	
	
TSMSolve$[adf_LKFPar, var_List, par_List, n_Integer]:=
	Module[{npar,lpar,rpar, tpar, liter, rvar, nvar, riter, niter,
				item, itemk , fun, link, ks, kss, fin, sal},
		npar = Length[adf[[3]]];
		rpar = Map[(LKF$Par[#]->par[[#]])&, Range[adf[[2]]]];
		lpar = Map[LKF$Par, Range[adf[[2]]+1,adf[[2]]+npar]];
		lpar = lpar //.Take[RightIterationLKF[adf],npar]/.fromLKF$Link/.rpar;
		tpar = par;
		Map[AppendTo[tpar,#]&,lpar]; 
		tpar = tpar /. {LKF$Constant[x_]:>x}; 
	
		ChangeIndexDer[];
		Map[depthDer,IterationLKF[adf]];

	
		liter = adf[[4]]; 
		niter = Length[liter]; 

		nvar = {};
		Do[ AppendTo[nvar,PWS[var[[i]]]], {i,adf[[1]]}];
		nvar[[1]] = PWS @@ Join[{var[[1]],1},Table[0, {n-1}]];
	
		rpar = Map[(LKF$Par[#]->tpar[[#]])&, Range[Length[tpar]]];
		rvar := Map[(LKF$Var[#]->nvar[[#]])&, Range[Length[var]]];
		
		link = Table[PWS[],{niter}]; 
		riter := Map[(LKF$Link[#]->link[[#]])&,Range[niter]]; 
		
		Do[ 	
			If[k > 0,
				Do[ itemk = variables$[adf[[5,i]],nvar[[i+1]],ks]/.
						rpar/.rvar/.riter/.ks->k; 
					AppendTo[nvar[[i+1]],itemk],
					{i,adf[[1]]-1}]]; 
			Do[ item = liter[[i]]/.rADheads/.{LKF$Constant[x_]:>x}; 
				Map[AppendTo[item,#]&,{LKF$Link[i],ks}];
				kss = k - depthDer[LKF$Link[i]];
				If[kss >= 0, 
					itemk = item/.rpar/.rvar/.riter/.ks->kss; 
					AppendTo[link[[i]],itemk]],
				{i,niter} ],
		{k,0,n}]; 
		sal = Drop[nvar,1];
		fin = Map[#[[1]]&,Drop[adf[[5]],adf[[1]]-1]];
		Map[AppendTo[sal, link[[#]]]&, fin];
		sal = sal /. {LKF$Constant[x_]:>x}; 

		Clear[link];
		Clear[nvar];
		sal
	] 



(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`PWSTides`"]

EndPackage[]

Null
