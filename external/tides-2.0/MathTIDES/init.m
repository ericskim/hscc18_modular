(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES: Taylor Integrator of Differential EquationS*)


(* ::Text:: *)
(*Version : 2.00*)


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


(* ::Title:: *)
(*init*)


(* ::Section::Closed:: *)
(*Iniciaci\[OAcute]n*)


(* ::Subsection::Closed:: *)
(*Fichero: LKFunctions*)


DeclarePackage["MathTIDES`LKFunctions`", 
	{   "LKF","LKFC", "ToLKF","ToTaylorLKF","ExtractConstants",
		"LKFPar","LKFCPar", "ToLKFPar","ToTaylorLKFPar","IterationLKF",
		"ToLKFC","ToTaylorLKFC","ToLKFCPar","ToTaylorLKFCPar",
		"RightIterationLKF","LeftIterationLKF","LieDer",
		"NumberOfVariables", "NumberOfParameters", "NumberOfFunctions",
		"NumberOfLinks", "LinksVariables", "LinksFunctions",
		"strFunC","strFunF","strFunCmp","StringCNumber","StringFNumber",
		"ListCTextConstants","ListFTextConstants","ListDoubleConstants",
		"TextMPConstants","TextDPConstants","TextMPClearConstants",
		"Double", "Multiple"}];



(* ::Subsection::Closed:: *)
(*Fichero: Odes*)


DeclarePackage["MathTIDES`ODES`", 
	{   "FirstOrderODE$", "FirstOrderODE","HamiltonianToODE",
		"PotentialToODE","NthOrderODE", "Screen","Points", "Delta","Only","Until"}];


(* ::Subsection::Closed:: *)
(*Fichero: Iterations*)


DeclarePackage["MathTIDES`Iterations`", 
	{   "ListIndexFunDer","CountDerivatives","PreviousList",
		"PreviousIndexList", "IteratorsList","IteratorsListStar",
		"CompleteIteratorsList","CompleteIteratorsListStar",
		"SortListDer", "DerOutputList","ListIndexFunLastOrder",
		"NewDerivativesList","PartialDerivativesText"}];



(* ::Subsection::Closed:: *)
(*Fichero: MinimalCode*)


DeclarePackage["MathTIDES`MinimalCode`", 
	{"MinCCText","MinCHText","DriverMinC","MinFFText","DriverMinFortran",
		"gmecopyrightC","gmecopyrightF", "gmecopyrightDrC", "gmecopyrightDrF"}];



(* ::Subsection::Closed:: *)
(*Fichero: StandardCode*)


DeclarePackage["MathTIDES`StandardCode`", 
	{"StandardCText", "StandardHText", "DriverStdC", "DriverEventC"}];



(* ::Subsection::Closed:: *)
(*Fichero : Texts*)


DeclarePackage["MathTIDES`Texts`",
	{ "mincgen", "minhgen","minfgen","headstdDP","headstdMP"}];


(* ::Subsection::Closed:: *)
(*Fichero: Codes*)


DeclarePackage["MathTIDES`Codes`",
	{ "TSMCodeFiles", "CodeFiles", "PrecisionDigits","MinTIDES", 
	"Driver","ODEFiles","TIDESFiles","MpfrTIDES",
	"ParametersValue", "InitialConditions", "IntegrationPoints", 
	"Output", "DataMatrix", "Factor1","Factor2", "OutputCoefficients",
	"Factor3","MaxStepRatio","MinStepRatio","MaxIterationsNumber",
	"OrderIncrement","MinOrder","MaxOrder",	"RelativeTolerance",
	"AbsoluteTolerance", "DefectErrorControl", "AddFunctions", "AddPartials",
	"Optimization","KahanSummation","CompensatedHorner","EventsNumber",
	"FindZeros","FindExtrema","FindMinima","FindMaxima","EventTolerance",
	"WriteMinTIDESCFiles", "WriteMinTIDESFFiles", "WriteDPTIDESHFiles","WriteMPTIDESHFiles"}];



(* ::Subsection::Closed:: *)
(*Fichero : PWSTides*)


DeclarePackage["MathTIDES`PWSTides`",
	{ "PWS", "ZeroPWS","UnitPWS","OrderPWS","PWSeries", "TSMSolve"}];


(* ::Subsection::Closed:: *)
(*Mensaje Inicial y s\[IAcute]mbolos en el contexto general*)


mathTIDESVersion$ = "2.00"


mathTIDESInit$ = "    MathTIDES "<>mathTIDESVersion$<>"\n" <>
"    MathTIDES    is   part   of   the   TIDES   project.\n"<>
"    Copyright(C)  2010  Abad, A.,  Barrio, R.,  Blesa, F.  and  Rodriguez, M.\n";


FrameBox[RowBox[{StyleBox[mathTIDESInit$,FontFamily->"Geneva"],
	StyleBox[Hyperlink["http://gme.unizar.es/software/tides"]]}],
	Background->LightYellow]//DisplayForm
