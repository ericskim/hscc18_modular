(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`Iterations: Iterations for partial derivatives*)


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
(*Iterations*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`Iterations`",
				"MathTIDES`ODES`",
					"MathTIDES`LKFunctions`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
ListIndexFunDer,
CountDerivatives,
PreviousList,
PreviousIndexList,
IteratorsList,
IteratorsListStar,
CompleteIteratorsList,
CompleteIteratorsListStar,
SortListDer,
DerOutputList,
ListIndexFunLastOrder,
NewDerivativesList,
PartialDerivativesText
}


{NTuples, NTuplesLastOrder, OrderedNTuples, Kji, TIndex, VarOrder, AddDimensionList}


{MinusOnePosition,MinusOnePosition,MinusOneOptimum, ComputeBinomials, PreviousList}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`Iterations`*"]
Clear @@ Names["MathTIDES`Iterations`*"]


(* ::Section::Closed:: *)
(*Mensajes*)


Begin["`mess`"]


End[]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsubsection:: *)
(*n - tuplas*)


NTuples[nvar_,mord_]:=
	Select[Tuples[Range[0,mord], nvar],
		(Plus@@# <= mord)&]


NTuplesLastOrder[nvar_,mord_]:=
	Select[Tuples[Range[0,mord], nvar],
		(Plus@@# == mord)&]


OrderedNTuples[n_,mord_, Index_]:=
	Sort[ NTuples[n,mord], (Index[#1] <Index[#2])&]


Kji[j_, mu_List]:= 
	Sum[mu[[Length[mu]-i]],{i,0,j-1}]-1

TIndex[mu_List]:= 
	Sum[Binomial[Kji[j,mu]+j,j],{j,1,Length[mu]}]


ListIndexFunDer[nvar_,mord_]:=
	OrderedNTuples[nvar,mord, TIndex]


CountDerivatives[nvar_,mord_]:= 
	Binomial[mord+nvar,nvar]


ListIndexFunLastOrder[nvar_,mord_]:=
	Prepend[
		SortListDer[
			NTuplesLastOrder[nvar,mord]],Table[0,{nvar}]]


VarOrder[{mu__List}]:= 
	Module[{nv,  mord},
		nv = Map[Length,{mu}][[1]];
		mord = Max[Map[(Plus @@ #)&, {mu}]];
		{nv,mord}
	]


SortListDer[{mu__List}]:= 
	Module[{lt},
		lt = ListIndexFunDer[Sequence@@VarOrder[{mu}]];
		Sort[{mu}, (Position[lt,#1][[1,1]]<Position[lt,#2][[1,1]])&]
]/; Equal@@Map[Length,{mu}]


(* ::Subsubsection:: *)
(*Binomials *)


BinomialsInt[ord_]:= 
	Map[(BinTaylor[ord,#] = Binomial[ord,#];)&, 
		Range[0,ord]];


ComputeBinomials[max_]:=  
	( Map[BinomialsInt, Range[0,max]];)


DeleteBinomials[]:= Clear[BinTaylor]


BinomialProduct[li_List, lv_List]:=
	Inner[BinTaylor,li,lv,Times]


(* ::Subsubsection:: *)
(*Listas de anteriores*)


PreviousList[{mu__List}]:= 
	SortListDer[
		Union[
			Flatten[
				Map[PreviousList,{mu}],1]]]

PreviousList[mu_List]:= 
	Module[{nmu},
		nmu = Union[
			Flatten[
				Outer[List,Sequence @@ Map[Range[0,#]&, mu]],
					Length[mu]-1]];
		SortListDer[nmu]
	]


PreviousIndexList[mu_]:=
	Map[TIndex, PreviousList[mu]]


MinusOnePosition[mu_, pos_]:= 
	If[mu[[pos]] == 0 , {}, 
		{ReplacePart[mu,pos->mu[[pos]]-1]}][[1]]

FirstNonZeroMinimum[mu_List]:=
	Module[{munc, min},
		munc = Select[mu,(#!=0)&];
		min = Min[munc];
		Position[mu, min][[1]]
	]

MinusOneOptimum[mu_]:= 
	MinusOnePosition[mu,FirstNonZeroMinimum[mu][[1]]]


IteratorsList[mu_]:=
	Module[{la , lc, pb},
		la = Union[Flatten[Outer[List,
			Sequence @@ Map[Range[0,#]&, mu]],Length[mu]-1]];
		lc = Map[(mu-#)&, la];
		pb = Map[BinomialProduct [mu, #]&, la];
		{pb,Map[TIndex,la],Map[TIndex, lc]}]


IteratorsListStar[mu_]:=
	{{1},{0},{0}}/; mu === Table[0,{Length[mu]}]

IteratorsListStar[mu_]:=
	Module[{me, la , lc, pb},
		me = MinusOneOptimum[mu];
		la = Union[Flatten[Outer[List,
			Sequence @@ Map[Range[0,#]&, me]],Length[mu]-1]];
		lc = Map[(mu-#)&, la];
		pb = Map[BinomialProduct [me, #]&, la];
		{pb,Map[TIndex,la],Map[TIndex, lc]}]


IteratorsList[nvar_, ordmax_]:=
	Map[IteratorsList, ListIndexFunDer[nvar,ordmax]]

IteratorsList[{mu__List}]:= 
	Map[IteratorsList, {mu}]


IteratorsListStar[nvar_, ordmax_]:=
	Map[IteratorsListStar, ListIndexFunDer[nvar,ordmax]]

IteratorsListStar[{mu__List}]:= 
	Map[IteratorsListStar, {mu}]


(* ::Subsubsection:: *)
(*Auxiliares listas derivadas parciales*)


PositionVARPAR[var_List, par_List, der_]:=
	Module[{pv,pp, num},
		pv = Position[var, der];
		pp = Position[par,der];
		num = If[pv == {}, pp[[1,1]]+Length[var],pv[[1,1]]];
		num] 

DerivativesList[var_List, par_List, der_List]:=
	Module[{dv,dp, dx, ndv,ndp},
		dv = Intersection[var,der];
		dp =  Intersection[par,der];
		dx =Complement[der, Join[dv,dp]];
		If[Length[dx] != 0, Message[FirstOrderODE$::"badder"]; Abort[]];
		Map[PositionVARPAR[var,par,#]&,der]]


OnlyUntilList[{}]={}
OnlyUntilList[{a_List, n_, s_Symbol}]:= {a, s[Length[a], n]}
OnlyUntilList[{a_List, n_Integer}]:= {a, Until[Length[a],n]}
OnlyUntilList[{a_List, n_List}]:= {a, Only[Length[a], n]}


NewDerivativesList[var_List, par_List, pwrt_List]:=
	Module[{ldir, sal, dim}, 
		ldir = If[pwrt =!= {}, 
					DerivativesList[var,par,pwrt[[1]]]];
		sal = OnlyUntilList[If[pwrt == {}, {}, 
					pwrt/.{pwrt[[1]]->ldir}]];
		If[sal != {} && Head[sal[[2,2]]] == List, 
			dim = Union[Union[Map[Length,sal[[2,2]]]]];
			If[Length[dim] != 1 || sal[[2,1]] =!= dim[[1]], 
				Message[FirstOrderODE$::"badderord"]; Abort[]]];
		sal
	]


(* ::Subsection::Closed:: *)
(*Iteraciones generales*)


DerOutputList[Until[nvar_Integer, mord_Integer]]:=
	DerOutputList[nvar, mord]

DerOutputList[Only[nvar_Integer, mord_Integer]]:=
	DerOutputList[ListIndexFunLastOrder[nvar, mord]]

DerOutputList[Until[nvar_Integer, ord_List]]:=
	DerOutputList[ord, All]

DerOutputList[Only[nvar_Integer, ord_List]]:=
	DerOutputList[ord]


DerOutputList[nvar_, ordmax_]:=
	Module[{ ls,  lo, cd},
		cd = CountDerivatives[nvar,ordmax];
		ls = ListIndexFunDer[nvar,ordmax];
		lo = Map[TIndex, ListIndexFunDer[nvar,ordmax]];
		{cd, ls, lo}
	 ]

DerOutputList[{mu__List}, All]:=
	Module[{lmu, nmu,  omu},
		nmu = PreviousList[{mu}];
		lmu = Length[nmu];
		omu = Range[0,lmu-1];
		{lmu, nmu, omu}  ]

DerOutputList[{mu__List}]:=
	Module[{lmu, nmu, xmu, pmu, omu},
		nmu = PreviousList[{mu}];
		lmu = Length[nmu];
		xmu = {mu};
		PrependTo[xmu,Table[0,{Length[xmu[[1]]]}]];
		xmu = Union[xmu];
		pmu = SortListDer[xmu];
		omu = Map[Position[nmu, #][[1,1]]&,pmu];
		{lmu, pmu, omu-1}  ]


CompleteIteratorsList[Until[nvar_Integer, mord_Integer]]:=
	CompleteIteratorsList[nvar, mord]

CompleteIteratorsList[Only[nvar_Integer, mord_Integer]]:=
	CompleteIteratorsList[ListIndexFunLastOrder[nvar, mord]]

CompleteIteratorsList[Until[nvar_Integer, ord_List]]:=
	CompleteIteratorsList[ord]

CompleteIteratorsList[Only[nvar_Integer, ord_List]]:=
	CompleteIteratorsList[ord]


CompleteIteratorsListStar[Until[nvar_Integer, mord_Integer]]:=
	CompleteIteratorsListStar[nvar, mord]

CompleteIteratorsListStar[Only[nvar_Integer, mord_Integer]]:=
	CompleteIteratorsListStar[ListIndexFunLastOrder[nvar, mord]]

CompleteIteratorsListStar[Until[nvar_Integer, ord_List]]:=
	CompleteIteratorsListStar[ord]

CompleteIteratorsListStar[Only[nvar_Integer, ord_List]]:=
	CompleteIteratorsListStar[ord]


CompleteIteratorsList[nvar_, ordmax_]:=
	Module[{ todos, acum, coef, lv,lc, cd},
		cd = CountDerivatives[nvar,ordmax];
		ComputeBinomials[ordmax];
		todos = IteratorsList[nvar,ordmax];
		acum = Accumulate[Map[Length[#[[1]]]&, todos]];
		PrependTo[acum, 0];
		coef = Flatten[Map[#[[1]]&, todos]];
		lv = Flatten[Map[#[[2]]&, todos]];
		lc = Flatten[Map[#[[3]]&, todos]];
		{acum, coef, lv, lc}
	 ]


CompleteIteratorsListStar[nvar_, ordmax_]:=
	Module[{ todos, acum, coef, lv,lc, cd},
		cd = CountDerivatives[nvar,ordmax];
		ComputeBinomials[ordmax];
		todos = IteratorsListStar[nvar,ordmax]; 
		acum = Accumulate[Map[Length[#[[1]]]&, todos]];
		PrependTo[acum, 0];  
		coef = Flatten[Map[#[[1]]&, todos]]; 
		lv =  Flatten[Map[#[[2]]&, todos]];
		lc =  Flatten[Map[#[[3]]&, todos]]; 
		{acum, coef, lv, lc}
	 ]


CompleteIteratorsList[{mu__List}]:=
	Module[{nvar,ordmax, todos, acum, coef, lv,lc, cd, pl,mpl},
		{nvar,ordmax} = VarOrder[{mu}];
		pl = PreviousList[{mu}];
		cd = Length[pl];
		ComputeBinomials[ordmax];
		todos = IteratorsList[pl];
		acum = Accumulate[Map[Length[#[[1]]]&, todos]];
		PrependTo[acum, 0];
		coef = Flatten[Map[#[[1]]&, todos]];
		mpl = Map[TIndex, pl];
		rules = Inner[Rule, mpl, Range[0,Length[mpl]-1],List];
		lv = Flatten[Map[#[[2]]&, todos]]/.rules;
		lc = Flatten[Map[#[[3]]&, todos]]/.rules;
		{acum, coef, lv, lc}
	 ]


CompleteIteratorsListStar[{mu__List}]:=
	Module[{nvar,ordmax, todos, acum, coef, lv,lc, cd, pl,mpl},
		{nvar,ordmax} = VarOrder[{mu}];
		pl = PreviousList[{mu}];
		cd = Length[pl];
		ComputeBinomials[ordmax];
		todos = IteratorsListStar[pl];
		acum = Accumulate[Map[Length[#[[1]]]&, todos]];
		PrependTo[acum, 0];
		coef = Flatten[Map[#[[1]]&, todos]];
		mpl = Map[TIndex, pl];
		rules = Inner[Rule, mpl, Range[0,Length[mpl]-1],List];
		lv = Flatten[Map[#[[2]]&, todos]]/.rules;
		lc = Flatten[Map[#[[3]]&, todos]]/.rules;
		{acum, coef, lv, lc}
	 ]


(* ::Subsubsection:: *)
(*Output Text*)


PartialDerivativesText[var_List, par_List, pwrt_List, name_String]:=
	Module[{ndl,np, pp,oder,tdl,cil,cils, std, ndlist, stn, texto = ""}, 
		ndl = NewDerivativesList[var,par,pwrt];
		{np, pp,oder} = If[ndl[[1]]=={}, 
			{0,{},Until[ Length[var],0]}, {ndl[[2,1]],ndl[[1]],ndl[[2]]}];
		tdl   =  DerOutputList[oder];
		cil   =  AddDimensionList[CompleteIteratorsList[oder]];
		cils =   AddDimensionList[CompleteIteratorsListStar[oder]];
		std = Flatten[{AddDimensionList[pp],tdl[[1]], cil, cils}];
		stn = Flatten[tdl[[2]]];
		ndlist = Join[std,stn];
		texto = texto <> "int   "<>name <> "_PDData[]  = " <> ToString[ndlist] <> ";\n";
		texto
	]

AddDimensionList[{lx__List}]:= Map[AddDimensionList,{lx}]
AddDimensionList[{}]:= {0}
AddDimensionList[lx_List]:= Prepend[lx,Length[lx]]


(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`Iterations`"]

EndPackage[]

Null
