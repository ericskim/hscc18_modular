(* ::Package:: *)

(* ::Title:: *)
(*MathTIDES`StandardCode:  C Standard TIDES Code*)


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
(*StandardCode*)


(* ::Section::Closed:: *)
(*Contexto y diccionario*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


BeginPackage["MathTIDES`StandardCode`",
				"MathTIDES`Iterations`",
					"MathTIDES`ODES`",
						"MathTIDES`Texts`",
							"MathTIDES`LKFunctions`",
								"MathTIDES`MinimalCode`"]


(* ::Subsection::Closed:: *)
(*S\[IAcute]mbolos*)


{
StandardCText, 
StandardHText, 
DriverStdC, 
DriverEventC,
SinCosLKFList$, 
SinCoshLKFList$
}


(* ::Subsection::Closed:: *)
(*Protecci\[OAcute]n*)


Unprotect @@ Names["MathTIDES`StandardCode`*"]
Clear @@ Names["MathTIDES`StandardCode`*"]


(* ::Section::Closed:: *)
(*Mensajes*)


Begin["`mess`"]


End[]


(* ::Section::Closed:: *)
(*C\[OAcute]digo*)


(* ::Subsection::Closed:: *)
(*Comienzo*)


Begin["`code`"]


(* ::Subsection::Closed:: *)
(*Constantes e Iteraciones*)


ListOfPartials[fun_FirstOrderODE$]:=
	If[fun[[5]] =={}, {}, fun[[5,1]]]

NumberOfPartials[fun_FirstOrderODE$]:= 
	Length[ListOfPartials[fun]]


$ParametersValue = Null 
$InitialConditions = Null 
$IntInt = Null 
$FileOutput = False 
$FileCoefficients = False ;
$DataMatrix = False 
$Factor1 = Null 
$Factor2 = Null 
$Factor3 = Null 
$MaxStepRatio = Null 
$MinStepRatio = Null 
$MaxIterationsNumber = Null 
$ExcessOrder = Null 
$MinOrder = Null 
$MaxOrder = Null
$RelativeTolerance = Null
$AbsoluteTolerance = Null
$StepSizeEstimator = False
$KahanSummation = True
$CompensatedHorner = False
$EventTolerance = Null
$EventVector = False
$MpfrTIDES = False


CleanMParams[]:=
(
	$ParametersValue = Null; 
	$InitialConditions = Null ;
	$IntInt = Null ;
	$FileOutput = False ;
	$FileCoefficients = False ;
	$DataMatrix = False ;
	$Factor1 = Null ;
	$Factor2 = Null ;
	$Factor3 = Null ;
	$MaxStepRatio = Null ;
	$MinStepRatio = Null ;
	$MaxIterationsNumber = Null ;
	$ExcessOrder = Null ;
	$MinOrder = Null; 
	$MaxOrder = Null;
	$RelativeTolerance = Null;
	$AbsoluteTolerance = Null;
	$StepSizeEstimator = False;
	$KahanSummation = True;
	$CompensatedHorner = False;
	$EventTolerance = Null;
	$EventVector = False;
	$MpfrTIDES = False;
)


StdIP[x_, 3999999999]:= Null

StdIP[Null,dp_]:= Null;
		
StdIP[{x_,Delta[dt_],Points[n_?Positive]},dp_]:=
	Module[{st0,std},
		If[dp>16,  
			st0 = strNUMBERC[x,dp];
			std = strNUMBERC[dt,dp],
			st0 = StringNumber[x];
			std =  StringNumber[dt]];
		{$TDN, st0, std, ToString[n]}
	]

StdIP[{x_,Points[n_?Positive],Delta[dt_]},dp_]:=
	StdIP[{x,Delta[dt] ,Points[n]},dp]

StdIP[{x_,y_,Points[n_?Positive]},dp_]:=
	Module[{st0,std},
		If[dp>16,  
			st0 = strNUMBERC[x,dp];
			std = strNUMBERC[y,dp], 
			st0 = StringNumber[x];
			std =  StringNumber[y]];
		{$TTN, st0, std, ToString[n]}
	]

StdIP[{x_,y_,Delta[dt_]}, dp_]:=
	Module[{st0,std,stf},
		If[dp>16,  
			st0 = strNUMBERC[x,dp];
			std = strNUMBERC[dt,dp];
			stf = strNUMBERC[y,dp],
			st0 = StringNumber[x];
			std =  StringNumber[dt];
			stf =  StringNumber[y]];
		{$TDT, st0, std, stf}
	]
	
StdIP[{x_,y_},dp_]:= 
	StdIP[{x,y,Points[1]},dp]

StdIP[x_List, dp_]:=
	Module[{stl},
		If[dp>16,  
			stl = Map[strNUMBERC[#,dp]&, x],
			stl = Map[StringNumber[#]&, x]];
		Flatten[{$TTT, stl}]
	]

StdIP[x_, dp_]:= Null



ValueMParams[params_, dp_]:=
(
	$ParametersValue = params[[1]]; 
	$InitialConditions = params[[2]] ;
	$IntInt = StdIP[params[[3]], dp] ;
	$FileOutput = params[[4]] ;
	$FileCoefficients = params[[5]] ;
	$DataMatrix = params[[6]] ;
	$Factor1 = params[[7]] ;
	$Factor2 = params[[8]] ;
	$Factor3 = params[[9]] ;
	$MaxStepRatio = params[[10]] ;
	$MinStepRatio = params[[11]] ;
	$MaxIterationsNumber = params[[12]] ;
	$ExcessOrder = params[[13]] ;
	$MinOrder = params[[14]]; 
	$MaxOrder = params[[15]];
	$RelativeTolerance = params[[16]];
	$AbsoluteTolerance = params[[17]];
	$StepSizeEstimator = params[[18]];
	$KahanSummation = params[[19]];
	$CompensatedHorner = params[[20]];
	$EventTolerance = params[[21]];
	$EventVector = params[[6]];
	$MpfrTIDES = params[[22]];
)


(* ::Subsection::Closed:: *)
(*C\[OAcute]digo C*)


(* ::Subsubsection::Closed:: *)
(*Repeticiones codigo C*)


separador =
"/*********************************************************************************************/"


headFunction1CDP[name_?StringQ]:=
	"long  "<> name <> 
	"(iteration_data *itd, double t, double v[], double p[], int ORDER, double *cvfd)\n{\n"

headFunction1CMP[name_?StringQ]:=
	"long  "<> name <> 
	"(iteration_data *itd, mpfr_t t, mpfr_t v[], mpfr_t p[], int ORDER, mpfr_t *cvfd)\n{\n"


numColumns[]:=
	"\tif(ORDER < 0) return NUM_COLUMNS;\n\n" 


setFirstC2 := 
	"\tfor(i=0;  i<=ORDER; i++) {";


setLast1CDP = 
	"\n\t}\n\n"<>
	"\twrite_dp_solution();\n"
setLast1CMP = 
	"\n\t}\n\n"<>
	"\twrite_mp_solution();\n"

setLast2C = 
	"\n\treturn NUM_COLUMNS;\n}\n\n";


includehDPT[]:=  "#include \"dp_tides.h\"\n" 
includehMPT[]:=  "#include \"mp_tides.h\"\n" 

includehC[name_?StringQ]:=  "#include \""<> name <>".h\"\n" 

includeSystemhC[name_?StringQ]:= "#include <"<> name <>">\n"


EndMessage[name_?StringQ]:= 
	"Files \"" <> name <> ".h\" and \"" <> name <> ".c\" written on directory \""<> Directory[]<> "\"."


(* ::Subsubsection::Closed:: *)
(*Listas y constantes en C*)


ConstantStringC[tipo_?StringQ, name_?StringQ, n_Integer]:=
	Module[{texto = "\tstatic "},
		texto = texto <> tipo <> " " <> name;
		texto = texto <>  " = " <> ToString[n]<> ";\n";
		texto]


ListIntegersToStringC[tipo_?StringQ, name_?StringQ,{}]:=
	Module[{texto="\tstatic "},
		texto = texto <> tipo <>" "; 
		texto = texto <> name <>"["<>ToString[1]<>"] = {0};\n";
		texto]


ListIntegersToStringC[tipo_?StringQ, name_?StringQ,mu_List]:=
	Module[{texto="\tstatic ",lmu},
		lmu = Length[mu];
		texto = texto <> tipo <>" "; 
		texto = texto <> name <>"["<>ToString[lmu]<>"] = {";
		Do[texto = texto <> ToString[mu[[i]]]<>",", {i,lmu-1}];
		texto=  texto <> ToString[mu[[lmu]]]<> "};\n";
		texto]


ListStringsToStringC[name_?StringQ,mu_List]:=
	Module[{texto="\tstatic char* ",lmu},
		lmu = Length[mu];
		texto = texto <> name <>"["<>ToString[lmu]<>"] = {";
		Do[texto = texto <> "\""<>mu[[i]]<> "\""<>",", {i,lmu-1}];
		texto= texto <> "\""<>mu[[lmu]]<> "\""<>"};\n";
		texto]


(* ::Subsubsection::Closed:: *)
(*Variables e iteraciones C*)


StringNumber[x_]:= StringCNumber[x]


strNUMBERC[x_,  3999999999]:= 
	"\"****\""

strNUMBERC[x_,  nd_]:= 
	"\""<> ToString[CForm[N[x,nd]]]<>"\""


strTayCDB[LKF$Plus]:= "double_add_t";
strTayCDB[LKF$Minus]:= "double_sub_t";
strTayCDB[LKF$Times]:= "double_mul_t";
strTayCDB[LKF$Inv]:= "double_inv_t";
strTayCDB[LKF$Divide]:= "double_div_t";
strTayCDB[LKF$DivideCV]:= "double_div_t_cv";
strTayCDB[LKF$DivideVC]:= "double_div_t_vc";
strTayCDB[LKF$Power]:= "double_pow_t";
strTayCDB[LKF$Sin]:= "double_sin_t";
strTayCDB[LKF$Cos]:= "double_cos_t";
strTayCDB[LKF$SinCos]:= "double_sin_cos_t";
strTayCDB[LKF$Tan]:= "double_tan_t";
strTayCDB[LKF$Sinh]:= "double_sinh_t";
strTayCDB[LKF$Cosh]:= "double_cosh_t";
strTayCDB[LKF$SinhCosh]:= "double_sinh_cosh_t";
strTayCDB[LKF$Tanh]:= "double_tanh_t";
strTayCDB[LKF$ArcSin]:= "double_asin_t";
strTayCDB[LKF$ArcCos]:= "double_acos_t";
strTayCDB[LKF$ArcTan]:= "double_atan_t";
strTayCDB[LKF$ArcSinh]:= "double_asinh_t";
strTayCDB[LKF$ArcCosh]:= "double_acosh_t";
strTayCDB[LKF$ArcTanh]:= "double_atanh_t";
strTayCDB[LKF$Log]:= "double_log_t";
strTayCDB[LKF$Exp]:= "double_exp_t";
strTayCDB[LKF$Der]:= "double_der_t";
strTayCDB[Var$]:= "double_var_t";
strTayCDB[Var$C]:= "double_var_t_cc";

strTayCMP[LKF$Plus]:= "mpfrts_add_t";
strTayCMP[LKF$Minus]:= "mpfrts_sub_t";
strTayCMP[LKF$Times]:= "mpfrts_mul_t";
strTayCMP[LKF$Inv]:= "mpfrts_inv_t";
strTayCMP[LKF$Divide]:= "mpfrts_div_t";
strTayCMP[LKF$DivideCV]:= "mpfrts_div_t_cv";
strTayCMP[LKF$DivideVC]:= "mpfrts_div_t_vc";
strTayCMP[LKF$Power]:= "mpfrts_pow_t";
strTayCMP[LKF$Sin]:= "mpfrts_sin_t";
strTayCMP[LKF$Cos]:= "mpfrts_cos_t";
strTayCMP[LKF$SinCos]:= "mpfrts_sin_cos_t";
strTayCMP[LKF$Tan]:= "mpfrts_tan_t";
strTayCMP[LKF$Sinh]:= "mpfrts_sinh_t";
strTayCMP[LKF$Cosh]:= "mpfrts_cosh_t";
strTayCMP[LKF$SinhCosh]:= "mpfrts_sinh_cosh_t";
strTayCMP[LKF$Tanh]:= "mpfrts_tanh_t";
strTayCMP[LKF$ArcSin]:= "mpfrts_asin_t";
strTayCMP[LKF$ArcCos]:= "mpfrts_acos_t";
strTayCMP[LKF$ArcTan]:= "mpfrts_atan_t";
strTayCMP[LKF$ArcSinh]:= "mpfrts_asinh_t";
strTayCMP[LKF$ArcCosh]:= "mpfrts_acosh_t";
strTayCMP[LKF$ArcTanh]:= "mpfrts_atanh_t";
strTayCMP[LKF$Log]:= "mpfrts_log_t";
strTayCMP[LKF$Exp]:= "mpfrts_exp_t";
strTayCMP[LKF$Der]:= "mpfrts_der_t";
strTayCMP[Var$]:= "mpfrts_var_t";
strTayCMP[Var$C]:= "mpfrts_var_t_cc";



strOBJC[LKF$Var[n_]]:= "var[" <> ToString[n-1] <> "]";
strOBJC[LKF$Par[n_]]:= "par[" <> ToString[n-1] <> "]";
strOBJC[LKF$Link[n_]]:= "link[" <> ToString[n-1] <> "]";
strOBJC[LKF$Const[n_]]:= "ct[" <> ToString[n-1] <> "]";


restaDepthDer[i_]:= 
	Module[{dpt, text},
		dpt = depthDer[LKF$Link[i]];
		If[dpt === 0, 
			text = "",
			text = "-" <> ToString[dpt] ];
		text
	]
strLITERC[ LKF$Sin[x_,LKF$Link[j_]], i_ ]:= strLITERCSinCos[i,j,x]
strLITERC[ LKF$Cos[x_,LKF$Link[j_]], i_ ]:= strLITERCSinCos[j,i,x]
strLITERCSinCos[i_,j_,x_]:=
	Module[{prev, sal},
		prev = MemberQ[SinCosLKFList$, {i,j}];
		If[prev , 
			sal = ""; 
			Drop[SinCosLKFList$, Position[SinCosLKFList$,{i,j}][[1]]],
			sal = strLITERCSinCosPrint[i,j,x]; 
			AppendTo[SinCosLKFList$,{i,j}]];
		sal
	]
strLITERCSinCosPrint[i_,j_,x_]:=
	"\n\t\t" <> strTayC[LKF$SinCos] <> "(itd, " <> strOBJC[x] <> "," <> strOBJC[LKF$Link[i]] <>
	"," <> strOBJC[LKF$Link[j]] <> ",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Sinh[x_,LKF$Link[j_]], i_ ]:= strLITERCSinCosh[i,j,x]
strLITERC[ LKF$Cosh[x_,LKF$Link[j_]], i_ ]:= strLITERCSinCosh[j,i,x]
strLITERCSinCosh[i_,j_,x_]:=
	Module[{prev, sal},
		prev = MemberQ[SinCoshLKFList$, {i,j}];
		If[prev , 
			sal = ""; 
			Drop[SinCoshLKFList$, Position[SinCoshLKFList$,{i,j}][[1]]],
			sal = strLITERCSinCoshPrint[i,j,x]; 
			AppendTo[SinCoshLKFList$,{i,j}]];
		sal
	]
strLITERCSinCoshPrint[i_,j_,x_]:=
	"\n\t\t" <> strTayC[LKF$SinhCosh] <> "(itd, " <> strOBJC[x] <> "," <> strOBJC[LKF$Link[i]] <>
	"," <> strOBJC[LKF$Link[j]] <> ",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Divide[1,z_], i_ ]:= 
	"\n\t\t" <> strTayC[LKF$Inv] <> "(itd, " <> strOBJC[z] <> "," <> strOBJC[LKF$Link[i]] <> 
	",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Divide[x_?numberADQ,z_], i_ ]:= 
	"\n\t\t" <> strTayC[LKF$DivideCV] <> "(itd, " <> strOBJC[x] <> "," <> strOBJC[z] <> ","<> strOBJC[LKF$Link[i]] <> 
	",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Divide[x_,z_?numberADQ], i_ ]:= 
	"\n\t\t" <> strTayC[LKF$DivideVC] <> "(itd, " <> strOBJC[x] <> "," <> strOBJC[z] <> ","<> strOBJC[LKF$Link[i]] <> 
	",i" <> restaDepthDer[i] <> ");";	

strLITERC[ LKF$Power[E, y_], i_ ]:= 
	"\n\t\t" <> strTayC[LKF$Exp] <> "(itd, "  <>  StringJoin @@ Map[(strOBJC[#] <> ",")&,{y}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";

strLITERC[ LKF$Power[x_, y_?numberADQ], i_ ]:= 
	"\n\t\t" <> strTayC[LKF$Power] <> "_cc" <> "(itd, "  <> strOBJC[x] <> "," <> 
	StringJoin @@ Map[(strOBJC[#] <> ",")&,{y}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";

strLITERC[ h_[x_?numberADQ, y__], i_ ]:= 
	"\n\t\t" <> strTayC[h] <> "_cc" <> "(itd, "  <> strOBJC[x] <> "," <>
	StringJoin @@ Map[(strOBJC[#] <> ",")&,{y}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";

strLITERC[ LKF$Minus[y_,x_?numberADQ], i_ ]:= 
	strLITERC[ LKF$Plus[y,-x], i ]

strLITERC[ h_[y_,x_?numberADQ, k___], i_ ]:= 
	"\n\t\t" <> strTayC[h] <> "_cc" <> "(itd, "  <> strOBJC[x] <> "," <>
	StringJoin @@ Map[(strOBJC[#] <> ",")&,{y,k}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";

strLITERC[ h_[y___], i_ ]:= 
	"\n\t\t" <> strTayC[h] <> "(itd, "  <> 
	StringJoin @@ Map[(strOBJC[#] <> ",")&,{y}] <>
	strOBJC[LKF$Link[i]] <> ",i" <> restaDepthDer[i] <> ");";
	
strVITERC[it_?numberADQ,i_]:=
	"\n\t\t"<> strTayC[Var$C]<>"(itd, " <> strOBJC[it] <> ",var[" <> ToString[i] <> "], i);"
	
strVITERC[it_,i_]:=
	"\n\t\t"<> strTayC[Var$]<>"(itd, " <> strOBJC[it] <> ",var[" <> ToString[i] <> "], i);"


textVARSC[lf_]:=
	StringJoin @@ Map[strVITERC[lf[[#]],#]&, Range[Length[lf]]];
	
textITERSC[it_]:=
	StringJoin @@ Map[strLITERC[it[[#]],#]&, Range[Length[it]]];


(* ::Subsubsection::Closed:: *)
(*Cadena C*)


StandardCText[name_?StringQ, oldfun_, n_Integer]:=
	Module[{texto = "",lkf,lkv,iter, dgt},
		$digitsNUMTIDES = n;
		SinCosLKFList$ = {};
		SinCoshLKFList$ = {};
		lkf = ToTaylorLKF[oldfun];
		lkf = ExtractConstants[lkf];
		iter =  IterationLKF[lkf[[4]]];
		lkv  =  LinksVariables[lkf];
		ChangeIndexDer[];
		Map[depthDer,iter];

		If[n > 16, strTayC = strTayCMP, strTayC = strTayCDB];
		texto = texto <> gmecopyrightC[];
		texto = texto <> If[n > 16, includehMPT[], includehDPT[]];
		texto = texto <> includehC[name]<>"\n\n";
		texto = texto <> If[n > 16, headFunction1CMP[name], headFunction1CDP[name]];
		texto = texto <> "\n\tint i;\n";
		texto = texto <> ConstantStringC["int  " , "VARIABLES       ", NumberOfVariables[lkf]];
		texto = texto <> ConstantStringC["int  " , "PARAMETERS      ", NumberOfParameters[lkf]];
		texto = texto <> ConstantStringC["int  " , "FUNCTIONS       ", NumberOfFunctions[lkf]];
		texto = texto <> ConstantStringC["int  " , "LINKS           ", NumberOfLinks[lkf]];
		texto = texto <> ListIntegersToStringC["int  ", "POS_FUNCTIONS",LinksFunctions[lkf]];
		texto = texto <> "\n";
		texto = texto <> If[n > 16, "\tinitialize_mp_case();\n\n", "\tinitialize_dp_case();\n\n"];

		If[n > 16, 
			texto = texto <> TextMPConstants[lkf[[3]]]<>"\n",
			texto = texto <> TextDPConstants[lkf[[3]]]<>"\n"];

		texto =  texto <> setFirstC2;
		texto =  texto <> textVARSC[lkv] ;
		texto =  texto <> textITERSC[lkf[[4]]] ;

		texto = texto <> If[n > 16, setLast1CMP, setLast1CDP];
		If[n > 16,  
			texto = texto <> "\tclear_vpl();\n";
			texto =  texto <> TextMPClearConstants[lkf[[3]]];
		];
		texto =  texto <> setLast2C;
		texto		
	]


StandardHText[name_?StringQ, n_Integer]:=
	Module[{texto =""}, 
		texto = texto <> gmecopyrightC[];

		texto = texto <> "\n#ifndef "<> name <>"_tides_h\n";
		texto = texto <> "#define "<> name <>"_tides_h\n\n";
		If[n > 16,
			texto = texto <> "long  "<> name <> 
				"(iteration_data *itd, mpfr_t t, mpfr_t v[], mpfr_t p[], int ORDER, mpfr_t *cvfd);\n\n",
			texto = texto <> "long  "<> name <> 
				"(iteration_data *itd, double t, double v[], double p[], int ORDER, double *cvfd);\n\n"
		 ];
		texto = texto <> "#endif\n\n\n";
		texto
	]



(* ::Subsubsection::Closed:: *)
(*Drivers*)


TextDate[]:=
	Module[{dt, texto, mes,min, tmin},
	mes = {"January ", "February ", "March ", 
	"April ","May ", "June ", "July ", "August ", 
	"September ", "October ","November ", "December "};
	dt = Date[];
	texto = "This file has been created by MathTIDES ("<> Global`mathTIDESVersion$ <> ") ";
	texto = texto <> mes[[dt[[2]]]] <>ToString[dt[[3]]] <> ", ";
	texto = texto <> ToString[dt[[1]]]<>", ";
	min = Round[dt[[5]]];
	tmin = If[min < 10, "0"<>ToString[min],  ToString[min]];
	texto = texto <> ToString[dt[[4]]] <> ":" <>tmin;
	texto
]


NameDataMatrix[name_]:= 
	If[Head[$DataMatrix] === String, $DataMatrix, name<> "_DataMatrix"]

NameEventVector[name_]:= 
	If[Head[$EventVector] === String, $EventVector, name<> "_EventsMatrix"]


bloqueDBIncludes[name_]:= 
	"#include \"dp_tides.h\"\n" <>
	"#include \""<> name <> ".h\""

bloqueMPIncludes[name_]:= 
	"#include \"mpfr.h\"\n"<> 
	"#include \"mp_tides.h\"\n" <>
	"#include \""<> name <> ".h\""

bloqueMPFIncludes[name_]:= 
	bloqueMPIncludes[name]


bloquePDCExtern[name_, fun_]:= 
	If[fun[[5]] === {}, "", 
		"\n\n" <> PartialDerivativesText[fun[[3]], fun[[4]], fun[[5]], name]]



bloqueDBMainBegin[]:= "\n\nint main() {"

bloqueMPMainBegin[n_]:=
	Module[{texto = "\n\nint main() {"},
			texto = texto <> "\n\n\tint i;";
			texto = texto <> "\n\n/* --- SET PRECISION  ------------ */";
		If[n === 3999999999, 
			texto = texto <> "\n\tset_precision_digits(****);",
			texto = texto <> "\n\tset_precision_digits(" <> ToString[n] <> ");" ];
		texto]

bloqueMPFMainBegin[n_]:=bloqueMPMainBegin[n]



bloqueMainEnd[]:= 
	Module[{texto="\n\n/* --- END  ---------------------- */"},
		If[Head[$FileOutput] === String,
			texto = texto  <>  "\n\tfclose(fd);"];
		texto = texto <> "\n\treturn 0;";
		texto = texto <> "\n}\n\n\n";
		texto]


bloqueExterns[]:=
	Module[{texto=""},
		If[$Factor1 =!= Null || $Factor2 =!= Null || $Factor3 =!= Null ||
			$MaxStepRatio =!= Null || $MinStepRatio =!= Null ||
			$MaxIterationsNumber =!= Null ||$ExcessOrder =!= Null  ||
			$MinOrder =!= Null ||$MaxOrder =!= Null ||$StepSizeEstimator === True ||
			$KahanSummation =!= True || $CompensatedHorner =!= False, 
				texto = texto <> "\n\n/* --- CONSTANTS OF THE METHOD  -- */"];
		If[$Factor1 =!= Null,  texto = texto <> 
			"\n\tfac1 = "<>StringNumber[$Factor1]<>";"];
		If[$Factor2 =!= Null,  texto = texto <> "\n\tfac2 = " 
			StringNumber[$Factor2]<>";"];
		If[$Factor3 =!= Null,  texto = texto  
			"\n\tfac3 = "<>StringNumber[$Factor3]<>";"];
		If[$MaxStepRatio =!= Null,  texto = texto <> 
			"\n\trmaxstep = "<>StringNumber[$MaxStepRatio]<>";"];
		If[$MinStepRatio =!= Null,  texto = texto <> 
			"\n\trminstep = "<>StringNumber[$MinStepRatio]<>";"];
		If[$MaxIterationsNumber =!= Null,  texto = texto <> 
			"\n\tnitermax = "<>ToString[$MaxIterationsNumber]<>";"];
		If[$ExcessOrder =!= Null,  texto  = texto <> 
			"\n\tnordinc = "<>ToString[$ExcessOrder]<>";"];
		If[$MinOrder =!= Null,  texto = texto <> 
			"\n\tminord = "<>ToString[$MinOrder]<>";"];
		If[$StepSizeEstimator === True,  texto = texto <> 
			"\n\tdefect_error_control = 1;"];
		If[$KahanSummation === False,  texto = texto <> 
			"\n\tkahan_summation = 0;"];
		If[$CompensatedHorner === True,  texto = texto <> 
			"\n\tcompensated_horner = 1;"];
		texto
		]



bloqueDBParams[p_]:=
	Module[{texto = ""},
		texto = texto <> "\n\n/* --- PARAMETERS  --------------- */";
		texto = texto <> "\n\tint npar = " <> If[p > 0, ToString[p], "0"] <> ";" ;
		If[p > 0,
			texto = texto <> "\n\tdouble p[npar];";
			Map[(texto = texto <> "\n\tp[" <> ToString[#] <>
				 "] = "<> 
				If[$ParametersValue === Null, "*****", StringNumber[$ParametersValue[[#+1]]]] <>
				" ; ")&, Range[0,p-1]]];
		texto]


bloqueMPParams[p_,n_]:=
	Module[{texto = ""},
		texto = texto <> "\n\n/* --- PARAMETERS  --------------- */";
		texto = texto <> "\n\tint npar = " <> If[p > 0, ToString[p], "0"] <> ";" ;
		If[p > 0,
			texto = texto <> "\n\tmpfr_t p[npar];";
			texto = texto <> "\n\tfor(i=0; i<npar; i++) mpfrts_init(&p[i]);";
			Map[(texto = texto <> "\n\tmpfrts_set_str(&p[" <> ToString[#] <> "], "<> 
				If[$ParametersValue === Null, "\"*****\"", strNUMBERC[$ParametersValue[[#+1]],n]] <>
				"); ")&, Range[0,p-1]]];
		texto]

bloqueMPFParams[p_,n_]:=
	Module[{texto = ""},
		texto = texto <> "\n\n/* --- PARAMETERS  --------------- */";
		texto = texto <> "\n\tint npar = " <> If[p > 0, ToString[p], "0"] <> ";" ;
		If[p > 0,
			texto = texto <> "\n\tmpfr_t p[npar];";
			texto = texto <> "\n\tfor(i=0; i<npar; i++) mpfr_init2(p[i], TIDES_PREC);";
			Map[(texto = texto <> "\n\tmpfr_set_str(p[" <> ToString[#] <> "], "<> 
				If[$ParametersValue === Null, "\"*****\"", strNUMBERC[$ParametersValue[[#+1]],n]] <>
				", 10, TIDES_RND); ")&, Range[0,p-1]]];
		texto]


bloqueDBVars[v_]:=
	Module[{texto = ""},
		texto = texto <> "\n\n/* --- VARIABLES   --------------- */";
		texto = texto <> "\n\tint nvar = " <> ToString[v] <> ";" ;
		texto = texto <> "\n\tdouble v[nvar];";
		Map[(texto = texto <>  "\n\tv[" <> ToString[#] <>
				 "] = "<> 
				If[$InitialConditions === Null, "*****", StringNumber[$InitialConditions[[#+1]]]] <>
				" ; ")&, Range[0,v-1]];
		texto]

bloqueMPVars[v_,n_]:=
	Module[{texto = ""},
		texto = texto <> "\n\n/* --- VARIABLES   --------------- */";
		texto = texto <> "\n\tint nvar = " <> ToString[v] <> ";" ;
		texto = texto <> "\n\tmpfr_t v[nvar];";
		texto = texto <> "\n\tfor(i=0; i<nvar; i++) mpfrts_init(&v[i]);";
		Map[(texto = texto <>  "\n\tmpfrts_set_str(&v[" <> ToString[#] <> "], "<> 
				If[$InitialConditions === Null, "\"*****\"", strNUMBERC[$InitialConditions[[#+1]],n]] <>
				"); ")&, Range[0,v-1]];
		texto
	]

bloqueMPFVars[v_,n_]:=
	Module[{texto = ""},
		texto = texto <> "\n\n/* --- VARIABLES   --------------- */";
		texto = texto <> "\n\tint nvar = " <> ToString[v] <> ";" ;
		texto = texto <> "\n\tmpfr_t v[nvar];";
		texto = texto <> "\n\tfor(i=0; i<nvar; i++) mpfr_init2(v[i], TIDES_PREC);";
		Map[(texto = texto <>  "\n\tmpfr_set_str(v[" <> ToString[#] <> "], "<> 
				If[$InitialConditions === Null, "\"*****\"", strNUMBERC[$InitialConditions[[#+1]],n]] <>
				", 10, TIDES_RND); ")&, Range[0,v-1]];
		texto
	]


bloqueFuns[nfun_]:=
	"\n\n/* --- NUMBER OF FUNCTIONS   ----- */" <> 
	"\n\tint nfun = " <> ToString[nfun] <> ";" ; 



bloqueNEVents[nevents_]:=
	"\n\n/* --- NUMBER OF EVENTS   -------- */" <> 
	"\n\tint nevents = " <> ToString[nevents] <> ";" ; 



bloqueDBTols[]:= 
	Module[{texto="", trel, tabs},
		If[$RelativeTolerance === Null && $AbsoluteTolerance === Null, 
			tabs = trel = N[10^-16]];
		If[$RelativeTolerance === Null && $AbsoluteTolerance =!= Null, 
			tabs = trel = N[$AbsoluteTolerance,16]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance === Null, 
			tabs = trel = N[$RelativeTolerance,16]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance =!= Null, 
			trel = N[$RelativeTolerance,16]; tabs = N[$AbsoluteTolerance,16]];
		texto = texto <> "\n\n/* --- TOLERANCES  --------------- */";
		texto = texto <> "\n\tdouble tolrel = "<> StringNumber[trel]<>" ;";
		texto = texto <> "\n\tdouble tolabs = "<> StringNumber[tabs]<>" ;";
		texto
	]

bloqueDBETols[]:= 
	Module[{texto="", tol},
		If[$EventTolerance === Null,  
			tol = N[10^-16],
			tol = N[$EventTolerance,16]];
		texto = texto <> "\n\n/* --- TOLERANCE  ---------------- */";
		texto = texto <> "\n\tdouble tol = "<> StringNumber[tol]<>" ;"; 
		texto
	]


bloqueMPTols[n_]:= 
	Module[{texto="", trel="", tabs=""},
		If[$RelativeTolerance === Null && $AbsoluteTolerance === Null, 
			If[n === 3999999999, 
				tabs = trel = strNUMBERC[0,n],
				tabs = trel = strNUMBERC[10^-(n-1),n]]];
		If[$RelativeTolerance === Null && $AbsoluteTolerance =!= Null, 
			tabs = trel = strNUMBERC[$AbsoluteTolerance,n]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance === Null, 
			tabs = trel = strNUMBERC[$RelativeTolerance,n]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance =!= Null, 
			trel = strNUMBERC[$RelativeTolerance,n]; 
			tabs = strNUMBERC[$AbsoluteTolerance,n]];
		texto = texto <> "\n\n/* --- TOLERANCES  --------------- */";
		texto = texto <> "\n\tmpfr_t tolrel, tolabs;";
		texto = texto <> "\n\tmpfrts_init(&tolrel); "; 
		texto = texto <> "\n\tmpfrts_init(&tolabs); "; 
		texto = texto <> "\n\tmpfrts_set_str(&tolrel, "<> trel <>");";
		texto = texto <> "\n\tmpfrts_set_str(&tolabs, "<> tabs <>");";
		texto
	]

bloqueMPETols[n_]:= 
	Module[{texto="", tol=""},
		If[$EventTolerance === Null , 
			If[n === 3999999999, 
				tol = strNUMBERC[0,n],
				tol = strNUMBERC[10^-(n-1),n]],
			tol = strNUMBERC[$EventTolerance,n]];
		texto = texto <> "\n\n/* --- TOLERANCE  ---------------- */";
		texto = texto <> "\n\tmpfr_t tol;";
		texto = texto <> "\n\tmpfrts_init(&tol); "; 
		texto = texto <> "\n\tmpfrts_set_str(&tol, "<> tol <>");";
		texto
	]

bloqueMPFTols[n_]:= 
	Module[{texto="", trel="", tabs=""},
		If[$RelativeTolerance === Null && $AbsoluteTolerance === Null, 
			If[n === 3999999999, 
				tabs = trel = strNUMBERC[0,n],
				tabs = trel = strNUMBERC[10^-(n-1),n]]];
		If[$RelativeTolerance === Null && $AbsoluteTolerance =!= Null, 
			tabs = trel = strNUMBERC[$AbsoluteTolerance,n]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance === Null, 
			tabs = trel = strNUMBERC[$RelativeTolerance,n]];
		If[$RelativeTolerance =!= Null && $AbsoluteTolerance =!= Null, 
			trel = strNUMBERC[$RelativeTolerance,n]; 
			tabs = strNUMBERC[$AbsoluteTolerance,n]];
		texto = texto <> "\n\n/* --- TOLERANCES  --------------- */";
		texto = texto <> "\n\tmpfr_t tolrel, tolabs;";
		texto = texto <> "\n\tmpfr_init2(tolrel, TIDES_PREC); "; 
		texto = texto <> "\n\tmpfr_init2(tolabs, TIDES_PREC); "; 
		texto = texto <> "\n\tmpfr_set_str(tolrel, "<> trel <>", 10, TIDES_RND);";
		texto = texto <> "\n\tmpfr_set_str(tolabs, "<> tabs <>", 10, TIDES_RND);";
		texto
	]

bloqueMPEFTols[n_]:= 
	Module[{texto="", tol=""},
		If[$EventTolerance === Null , 
			If[n === 3999999999, 
				tol = strNUMBERC[0,n],
				tol = strNUMBERC[10^-(n-1),n]],
			tol = strNUMBERC[$EventTolerance,n]];
		texto = texto <> "\n\n/* --- TOLERANCE  ---------------- */";
		texto = texto <> "\n\tmpfr_t tol;";
		texto = texto <> "\n\tmpfr_init2(tol, TIDES_PREC); "; 
		texto = texto <> "\n\tmpfr_set_str(tol, "<> tol <>", 10, TIDES_RND);";
		texto
	]



bloqueDB$IntPt[{$TDN, t0_String, dt_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tdouble tini = " <> t0 <> ";";
		texto = texto <> "\n\tdouble dt   = " <> dt <> ";";
		texto = texto <> "\n\tint    nipt = " <> n <> ";";
		texto
	]

bloqueMP$IntPt[{$TDN, t0_String, dt_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, dt; " ;
		texto = texto <> "\n\tmpfrts_init(&tini); " ;
		texto = texto <> "\n\tmpfrts_init(&dt); " ;
		texto = texto <> "\n\tmpfrts_set_str(&tini, "<> t0 <>");" ;
		texto = texto <> "\n\tmpfrts_set_str(&dt, "<> dt <>");" ;
		texto = texto <> "\n\tint  nipt  = " <> n <> ";";
		texto
	]

bloqueMPF$IntPt[{$TDN, t0_String, dt_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, dt; " ;
		texto = texto <> "\n\tmpfr_init2(tini, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_init2(dt, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_set_str(tini, "<> t0 <>", 10, TIDES_RND);" ;
		texto = texto <> "\n\tmpfr_set_str(dt, "<> dt <>", 10, TIDES_RND);" ;
		texto = texto <> "\n\tint  nipt  = " <> n <> ";";
		texto
	]

bloqueDB$IntPt[{$TDT, t0_String, dt_String, tf_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tdouble tini = " <> t0 <> ";";
		texto = texto <> "\n\tdouble dt   = " <> dt <> ";";
		texto = texto <> "\n\tdouble tend = " <> tf <> ";";
		texto = texto <> "\n\tint    nipt = (int) floor ((tend-tini)/dt);";
		texto
	]

bloqueMP$IntPt[{$TDT, t0_String, dt_String, tf_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, dt, tend; " ;
		texto = texto <> "\n\tmpfrts_init(&tini); " ;
		texto = texto <> "\n\tmpfrts_init(&dt); " ;
		texto = texto <> "\n\tmpfrts_init(&tend); " ;
		texto = texto <> "\n\tmpfrts_set_str(&tini, "<> t0 <>");" ;
		texto = texto <> "\n\tmpfrts_set_str(&dt, "<> dt <>");" ;
		texto = texto <> "\n\tmpfrts_set_str(&tend, "<> tf <>");" ;
		texto = texto <> "\n\tmpfr_t aux; " ;
		texto = texto <> "\n\tmpfrts_init(&aux); " ;
		texto = texto <> "\n\tmpfrts_set(&aux,tend); " ;
		texto = texto <> "\n\tmpfrts_sub(&aux,aux,tini); " ;
		texto = texto <> "\n\tmpfrts_div(&aux,aux,dt); " ;
		texto = texto <> "\n\tdouble auxd = mpfrts_get_d(aux); " ;
		texto = texto <> "\n\tint    nipt = (int) floor(auxd);";
		texto
	]

bloqueMPF$IntPt[{$TDT, t0_String, dt_String, tf_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, dt, tend; " ;
		texto = texto <> "\n\tmpfr_init2(tini, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_init2(dt, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_init2(tend, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_set_str(tini, "<> t0 <>", 10, TIDES_RND);" ;
		texto = texto <> "\n\tmpfr_set_str(dt, "<> dt <>", 10, TIDES_RND);" ;
		texto = texto <> "\n\tmpfr_set_str(tend, "<> tf <>", 10, TIDES_RND);" ;
		texto = texto <> "\n\tmpfr_t aux; " ;
		texto = texto <> "\n\tmpfr_init2(aux, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_set(aux,tend, TIDES_RND); " ;
		texto = texto <> "\n\tmpfr_sub(aux,aux,tini, TIDES_RND); " ;
		texto = texto <> "\n\tmpfr_div(aux,aux,dt, TIDES_RND); " ;
		texto = texto <> "\n\tdouble auxd = mpfr_get_d(aux, TIDES_RND); " ;
		texto = texto <> "\n\tint    nipt = (int) floor(auxd);";
		texto
	]

bloqueDB$IntPt[{$TTN, t0_String, tf_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tdouble tini = " <> t0 <> ";";
		texto = texto <> "\n\tdouble tend = " <> tf <> ";";
		texto = texto <> "\n\tint    nipt = " <> n <> ";";
		texto = texto <> "\n\tdouble dt = (tend - tini)/nipt ;";
		texto
	]

bloqueMP$IntPt[{$TTN, t0_String, tf_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, dt, tend; " ;
		texto = texto <> "\n\tmpfrts_init(&tini); " ;
		texto = texto <> "\n\tmpfrts_init(&dt); " ;
		texto = texto <> "\n\tmpfrts_init(&tend); " ;
		texto = texto <> "\n\tmpfrts_set_str(&tini, "<> t0 <>");" ;
		texto = texto <> "\n\tmpfrts_set_str(&tend, "<> tf <>");" ;
		texto = texto <> "\n\tint  nipt  = " <> n <> ";";
		texto = texto <> "\n\tmpfrts_set(&dt,tend); " ;
		texto = texto <> "\n\tmpfrts_sub(&dt,dt,tini); " ;
		texto = texto <> "\n\tmpfrts_div_i(&dt,dt,nipt); " ;
		texto
	]

bloqueMPF$IntPt[{$TTN, t0_String, tf_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, dt, tend; " ;
		texto = texto <> "\n\tmpfr_init2(tini, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_init2(dt, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_init2(tend, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_set_str(tini, "<> t0 <>", 10, TIDES_RND);" ;
		texto = texto <> "\n\tmpfr_set_str(tend, "<> tf <>", 10, TIDES_RND);" ;
		texto = texto <> "\n\tint  nipt  = " <> n <> ";";
		texto = texto <> "\n\tmpfr_set(dt,tend, TIDES_RND); " ;
		texto = texto <> "\n\tmpfr_sub(dt,dt,tini, TIDES_RND); " ;
		texto = texto <> "\n\tmpfr_div_si(dt,dt,nipt, TIDES_RND); " ;
		texto
	]

bloqueDB$IntPt[{$TTT, tt__String}]:=
	Module[{n , ltt, texto=""},
		ltt = {tt};
		n =  Length[ltt];
		texto = texto <> "\n\tint ntot  = " <> ToString[n]<> ";";
		texto = texto <> "\n\tdouble lt[ntot] ;";
		Map[(texto = texto <> "\n\tlt["<>ToString[#]<>"] = "
			<> ltt[[#+1]]<>" ;")&, Range[0,n-1]];
		texto
	]

bloqueMP$IntPt[{$TTT, tt__String}]:=
	Module[{n , ltt , texto=""},
		ltt = {tt};
		n =  Length[ltt];
		texto = texto <> "\n\tint ilt, ntot  = " <> ToString[n] <> ";";
		texto = texto <> "\n\tmpfr_t lt[ntot] ;";
		texto = texto <> "\n\tfor(ilt = 0; ilt < ntot; ilt++) mpfrts_init(&lt[ilt]);";
		Map[(texto = texto <> "\n\tmpfrts_set_str(&lt["<>ToString[#]<>"], "
			<> ltt[[#+1]]<>" );")&, Range[0,n-1]];
		texto
	]

bloqueMPF$IntPt[{$TTT, tt__String}]:=
	Module[{n , ltt , texto=""},
		ltt = {tt};
		n =  Length[ltt];
		texto = texto <> "\n\tint ilt, ntot  = " <> ToString[n] <> ";";
		texto = texto <> "\n\tmpfr_t lt[ntot] ;";
		texto = texto <> "\n\tfor(ilt = 0; ilt < ntot; ilt++) mpfr_init2(lt[ilt], TIDES_PREC);";
		Map[(texto = texto <> "\n\tmpfr_set_str(lt["<>ToString[#]<>"], "
			<> ltt[[#+1]]<>", 10, TIDES_RND);")&, Range[0,n-1]];
		texto
	]

bloqueDB$IntPt[Null]:=
	Module[{texto = ""},
		texto = texto <> "\n\tdouble tini = *****;";
		texto = texto <> "\n\tdouble dt   = *****;";
		texto = texto <> "\n\tint    nipt = *****;";
		texto
	]

bloqueMP$IntPt[Null]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, dt; " ;
		texto = texto <> "\n\tmpfrts_init(&tini); " ;
		texto = texto <> "\n\tmpfrts_init(&dt); " ;
		texto = texto <> "\n\tmpfrts_set_str(&tini, \"*****\");" ;
		texto = texto <> "\n\tmpfrts_set_str(&dt, \"*****\");" ;
		texto = texto <> "\n\tint  nipt  = *****;";
		texto
	]

bloqueMPF$IntPt[Null]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, dt; " ;
		texto = texto <> "\n\tmpfr_init2(tini, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_init2(dt, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_set_str(tini, \"*****\", 10, TIDES_RND);" ;
		texto = texto <> "\n\tmpfr_set_str(dt, \"*****\", 10, TIDES_RND);" ;
		texto = texto <> "\n\tint  nipt  = *****;";
		texto
	]



bloqueIntPoints[x_,n_]:=
	Module[{texto=""},
		texto = texto <> "\n\n/* --- INTEGRATION POINTS   ------ */";
		texto = texto <>If[n >16, 
			If[$MpfrTIDES, bloqueMP$IntPt[x], bloqueMPF$IntPt[x]] , 
			bloqueDB$IntPt[x]];
		texto
	]


bloqueDBE$IntPt[{$TTN, t0_String, tf_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tdouble tini = " <> t0 <> ";";
		texto = texto <> "\n\tdouble tend = " <> tf <> ";";
		texto
	]

bloqueMPE$IntPt[{$TTN, t0_String, tf_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, tend; " ;
		texto = texto <> "\n\tmpfrts_init(&tini); " ;
		texto = texto <> "\n\tmpfrts_init(&tend); " ;
		texto = texto <> "\n\tmpfrts_set_str(&tini, "<> t0 <>");" ;
		texto = texto <> "\n\tmpfrts_set_str(&tend, "<> tf <>");" ;
		texto
	]

bloqueMPEF$IntPt[{$TTN, t0_String, tf_String, n_String}]:=
	Module[{texto=""},
		texto = texto <> "\n\tmpfr_t tini, tend; " ;
		texto = texto <> "\n\tmpfr_init2(tini, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_init2(tend, TIDES_PREC); " ;
		texto = texto <> "\n\tmpfr_set_str(tini, "<> t0 <>", 10, TIDES_RND);" ;
		texto = texto <> "\n\tmpfr_set_str(tend, "<> tf <>", 10, TIDES_RND);" ;
		texto
	]

bloqueEIntPoints[x_,n_]:=
	Module[{texto=""},
		texto = texto <> "\n\n/* --- INTEGRATION POINTS   ------ */";
		texto = texto <>If[n >16, 
			If[$MpfrTIDES, bloqueMPE$IntPt[x], bloqueMPEF$IntPt[x]], 
			bloqueDBE$IntPt[x]];
		texto
	]


bloquePDC[name_, fun_]:= 
	Module[{texto=""},	
		If[fun[[5]] != {}, 
			texto = texto <> "\n\n/* - DECLARE PARTIAL DERIVATIVES - */";
			texto = texto <> "\n\tdeclare_partial_derivatives("<> name <> "_PDData);"];
		texto	
	]


bloqueDBOutput[name_]:= 
	Module[{texto=""},	
		texto = texto <> "\n\n/* --- OUTPUT  ------------------- */";
		If[$FileOutput === Screen , 
			texto = texto <> "\n\tFILE* fd = stdout;"];
		If[Head[$FileOutput] === String, 
			texto = texto <> "\n\tFILE* fd =" ;
			texto = texto <> " fopen(\""<>$FileOutput<>"\", \"w\");"];
		If[$DataMatrix =!= False, 
			texto =texto <> "\n\tdp_data_matrix "<> NameDataMatrix[name] <> ";"];
		texto	
	]

bloqueDBEOutput[name_,nevents_]:= 
	Module[{texto=""},	
		texto = texto <> "\n\n/* --- OUTPUT  ------------------- */";
		If[$FileOutput === Screen , 
			texto = texto <> "\n\tFILE* fd = stdout;"];
		If[Head[$FileOutput] === String, 
			texto = texto <> "\n\tFILE* fd =" ;
			texto = texto <> " fopen(\""<>$FileOutput<>"\", \"w\");"];
		If[$EventVector =!= False, 
			texto =texto <> "\n\tdp_data_matrix "<> NameEventVector[name] <> ";"];
		texto	
	]



bloqueMPOutput[name_]:= 
	Module[{texto=""},	
		texto = texto <> "\n\n/* --- OUTPUT  ------------------- */";
		If[$FileOutput === Screen , 
			texto = texto <> "\n\tFILE* fd = stdout;"];
		If[Head[$FileOutput] === String, 
			texto = texto <> "\n\tFILE* fd =" ;
			texto = texto <> " fopen(\""<>$FileOutput<>"\", \"w\");"];
		If[$DataMatrix =!= False, 
			texto =texto <>"\n\tmp_data_matrix "<> NameDataMatrix[name] <> ";"];
		texto	
	]

bloqueMPEOutput[name_,nevents_]:= 
	Module[{texto=""},	
		texto = texto <> "\n\n/* --- OUTPUT  ------------------- */";
		If[$FileOutput === Screen , 
			texto = texto <> "\n\tFILE* fd = stdout;"];
		If[Head[$FileOutput] === String, 
			texto = texto <> "\n\tFILE* fd =" ;
			texto = texto <> " fopen(\""<>$FileOutput<>"\", \"w\");"];
		If[$EventVector =!= False, 
			texto =texto <> "\n\tmp_data_matrix "<> NameEventVector[name] <> ";"];
		texto	
	]


bloqueMPFOutput[name_]:=bloqueMPOutput[name]
bloqueMPEFOutput[name_,nevents_]:= bloqueMPEOutput[name,nevents]


bloqueDBIntCall[p_,name_, fun_]:= 
	Module[{texto="", tdm = "", tfo, ptext, namepdd },
		namepdd = If[fun[[5]] != {}, name <> "_PDData", "NULL"];	
		tfo = If[$FileOutput === Screen || Head[$FileOutput] === String, "fd", "NULL"];
		tdm = If[$DataMatrix =!= False, "&"<>NameDataMatrix[name], "NULL"];
		ptext = If[p === 0, "NULL, ", "p, "];
		texto = texto <> "\n\n/* --- INTEGRATOR  --------------- */";
		If[$IntInt === Null || $IntInt[[1]] =!= $TTT, 
			texto = texto <>"\n\tdp_tides_delta("<> name <>", "<> namepdd <>
					", nvar, npar, nfun, v, "<> ptext<>"tini, dt, nipt, tolrel, tolabs, ",
			texto = texto <>"\n\tdp_tides_list("<> name <>  ", "<> namepdd <>
					", nvar, npar, nfun, v, "<> ptext<>"lt, ntot, tolrel, tolabs, "
			];
		texto = texto <> tdm <> ", "<> tfo <> ");";
		texto
	]

bloqueMPIntCall[p_,name_, fun_]:= 
	Module[{texto="", tdm = "", tfo, ptext, namepdd },	
		namepdd = If[fun[[5]] != {}, name <> "_PDData", "NULL"];	
		tfo = If[$FileOutput === Screen || Head[$FileOutput] === String, "fd", "NULL"];
		tdm = If[$DataMatrix =!= False, "&"<>NameDataMatrix[name], "NULL"];
		ptext = If[p === 0, "NULL, ", "p, "];
		texto = texto <> "\n\n/* --- INTEGRATOR  --------------- */";
		If[$IntInt === Null || $IntInt[[1]] =!= $TTT, 
			texto = texto <>"\n\tmp_tides_delta("<> name <>", "<> namepdd <>
					", nvar, npar, nfun, v, "<> ptext<>"tini, dt, nipt, tolrel, tolabs, ",
			texto = texto <>"\n\tmp_tides_list("<> name <>  ", "<> namepdd <>
					", nvar, npar, nfun, v, "<> ptext<>"lt, ntot, tolrel, tolabs, "];
		texto = texto <> tdm <> ", "<> tfo <> ");";
		texto
	]

bloqueMPFIntCall[p_,name_,fun_]:= bloqueMPIntCall[p,name,fun]


bloqueDBEventCall[name_,p_,fevent_]:= 
	Module[{texto="", head},	
		texto = texto <> "\n\n/* --- EVENT FINDER  ------------- */";
		head = Switch[fevent, 
			1, "dp_tides_find_zeros(",
			2, "dp_tides_find_extrema(",
			3, "dp_tides_find_minimum(",
			4, "dp_tides_find_maximum("];
		texto = texto <> "\n\t" <> head <> name <>", nvar, npar, v, ";
		texto = texto <> If[p === 0, "NULL, ", "p, "] ;
		texto = texto <> "tini, tend, tol, &nevents, "; 
		texto = texto <> 
			If[$EventVector === False, "NULL", "&"<>NameEventVector[name]] <> ", ";
		texto = texto <> 
			If[$FileOutput === False, "NULL", "fd"] <> ");\n";
		texto
	]

bloqueMPEventCall[name_,p_,fevent_]:= 
	Module[{texto="", head},	
		texto = texto <> "\n\n/* --- EVENT FINDER  ------------- */";
		head = Switch[fevent, 
			1, "mp_tides_find_zeros(",
			2, "mp_tides_find_extrema(",
			3, "mp_tides_find_minimum(",
			4, "mp_tides_find_maximum("];
		texto = texto <> "\n\t" <> head <> name <>", nvar, npar, v, ";
		texto = texto <> If[p === 0, "NULL, ", "p, "] ;
		texto = texto <> "tini, tend, tol, &nevents, "; 
		texto = texto <> 
			If[$EventVector === False, "NULL", "&"<>NameEventVector[name]] <> ", ";
		texto = texto <> 
			If[$FileOutput === False, "NULL", "fd"] <> ");\n";
		texto
	]

bloqueMPEventFCall[name_,p_,fevent_]:= bloqueMPEventCall[name,p,fevent]


DriverDB[name_?StringQ, v_, p_, nf_, fun_]:= 
	Module[{texto = ""},
		texto = gmecopyrightDrC["dp"];
		texto = texto <> bloqueDBIncludes[name]; 
		texto = texto <> bloquePDCExtern[name,fun];
		texto = texto <> bloqueDBMainBegin[];
		texto = texto <> bloqueExterns[];
		texto = texto <> bloqueDBParams[p];
		texto = texto <> bloqueDBVars[v];
		texto = texto <> bloqueFuns[nf];
		texto = texto <> bloqueDBTols[];
		texto = texto <> bloqueIntPoints[$IntInt, 16];
		texto = texto <> bloqueDBOutput[name];
		texto = texto <> bloqueDBIntCall[p,name,fun];
		texto = texto <> bloqueMainEnd[];
		texto
	]


DriverMP[name_?StringQ, n_, v_, p_, nf_, fun_, mpfrt_]:= 
	If[mpfrt, DriverMPFR[name,n,v,p, nf,fun], DriverMPTD[name,n,v,p, nf,fun]]

DriverMPFR[name_?StringQ, n_, v_, p_, nf_, fun_]:= 
	Module[{texto = ""},
		texto = gmecopyrightDrC["mp"];
		texto = texto <> bloqueMPIncludes[name];
		texto = texto <> bloquePDCExtern[name,fun];
		texto = texto <> bloqueMPMainBegin[n];
		texto = texto <> bloqueExterns[];
		texto = texto <> bloqueMPParams[p,n];
		texto = texto <> bloqueMPVars[v,n];
		texto = texto <> bloqueFuns[nf];
		texto = texto <> bloqueMPTols[n];
		texto = texto <> bloqueIntPoints[$IntInt, n];
		texto = texto <> bloqueMPOutput[name];
		texto = texto <> bloqueMPIntCall[p,name,fun];
		texto = texto <> bloqueMainEnd[];
		texto
	]

DriverMPTD[name_?StringQ, n_, v_, p_, nf_, fun_]:= 
	Module[{texto = ""},
		texto = gmecopyrightDrC["mp"];
		texto = texto <> bloqueMPFIncludes[name];
		texto = texto <> bloquePDCExtern[name,fun];
		texto = texto <> bloqueMPFMainBegin[n];
		texto = texto <> bloqueExterns[];
		texto = texto <> bloqueMPFParams[p,n];
		texto = texto <> bloqueMPFVars[v,n];
		texto = texto <> bloqueFuns[nf];
		texto = texto <> bloqueMPFTols[n];
		texto = texto <> bloqueIntPoints[$IntInt, n];
		texto = texto <> bloqueMPFOutput[name];
		texto = texto <> bloqueMPFIntCall[p,name,fun];
		texto = texto <> bloqueMainEnd[];
		texto
	]



DriverStdC[name_?StringQ, n_, fun_, params_List]:=
	Module[{v, p, np, nf, texto},
		ValueMParams[params,n]; 
		adf = ToTaylorLKF[fun];
		v = adf[[1]]-1;
		p = adf[[2]];
		nf = NumberOfFunctions[adf];
		texto = If[n >16, DriverMP[name,n,v,p, nf,fun,$MpfrTIDES], DriverDB[name,v,p, nf,fun]];
		CleanMParams[];
		texto
	]




DriverEventDB[name_?StringQ, v_, p_, nf_,nevents_, fevent_]:= 
	Module[{texto = ""},
		texto = gmecopyrightDrC["dp"];
		texto = texto <> bloqueDBIncludes[name];
		texto = texto <> bloqueDBMainBegin[];
		texto = texto <> bloqueExterns[];
		texto = texto <> bloqueDBParams[p];
		texto = texto <> bloqueDBVars[v];
		texto = texto <> bloqueNEVents[nevents];
		texto = texto <> bloqueDBETols[];
		texto = texto <> bloqueEIntPoints[$IntInt, 16];
		texto = texto <> bloqueDBEOutput[name,nevents];
		texto = texto <> bloqueDBEventCall[name,p,fevent];
		texto = texto <> bloqueMainEnd[];
		texto
	]


DriverEventMP[name_?StringQ, n_, v_, p_, nf_,nevents_, fevent_, mpfrt_]:= 
	If[mpfrt, 
		DriverEventMPFR[name,n,v,p, nf,nevents,fevent], 
		DriverEventMPTD[name,n,v,p, nf,nevents,fevent]]

DriverEventMPFR[name_?StringQ, n_, v_, p_, nf_,nevents_, fevent_]:= 
	Module[{texto = ""},
		texto = gmecopyrightDrC["mp"];
		texto = texto <> bloqueMPIncludes[name];
		texto = texto <> bloqueMPMainBegin[n];
		texto = texto <> bloqueExterns[];
		texto = texto <> bloqueMPParams[p,n];
		texto = texto <> bloqueMPVars[v,n];
		texto = texto <> bloqueNEVents[nevents];
		texto = texto <> bloqueMPETols[n];
		texto = texto <> bloqueEIntPoints[$IntInt, n];
		texto = texto <> bloqueMPEOutput[name,nevents];
		texto = texto <> bloqueMPEventCall[name,p,fevent];
		texto = texto <> bloqueMainEnd[];
		texto
	]

DriverEventMPTD[name_?StringQ, n_, v_, p_, nf_,nevents_, fevent_]:= 
	Module[{texto = ""},
		texto = gmecopyrightDrC["mp"];
		texto = texto <> bloqueMPFIncludes[name];
		texto = texto <> bloqueMPFMainBegin[n];
		texto = texto <> bloqueExterns[];
		texto = texto <> bloqueMPFParams[p,n];
		texto = texto <> bloqueMPFVars[v,n];
		texto = texto <> bloqueNEVents[nevents];
		texto = texto <> bloqueMPEFTols[n];
		texto = texto <> bloqueEIntPoints[$IntInt, n];
		texto = texto <> bloqueMPEFOutput[name,nevents];
		texto = texto <> bloqueMPEventFCall[name,p,fevent];
		texto = texto <> bloqueMainEnd[];
		texto
	]


DriverEventC[name_?StringQ, n_, fun_, params_List, nevents_, fevent_ ]:=
	Module[{v, p, np, nf, texto},
		ValueMParams[params,n]; 
		adf = ToTaylorLKF[fun];
		v = adf[[1]]-1;
		p = adf[[2]];
		nf = NumberOfFunctions[adf];
		texto = If[n >16, 
			DriverEventMP[name,n,v,p, nf,nevents,fevent, $MpfrTIDES], 
			DriverEventDB[name,v,p, nf,nevents,fevent]];
		CleanMParams[];
		texto
	]




(* ::Subsection::Closed:: *)
(*Final*)


End[]


(* ::Section::Closed:: *)
(*Final*)


Protect @@ Names["MathTIDES`StandardCode`"]

EndPackage[]

Null
