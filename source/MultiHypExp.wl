(* ::Package:: *)

(* ::Section:: *)
(*Version details*)


(* ::Text:: *)
(*v2.0 :: KdF and MKdFs can be computed with this version*)
(*v2.1 :: one variable functions are included*)
(*v2.2 :: three variable functions are included.*)
(*v2.3 :: comments added*)
(*v2.7 :: three variables cases are added FA,FB, FD, FN ,FS*)
(*v2.8 :: special cases of F4 are added*)
(*v3.0	:: New command "ReduceFunction" is added	*)


(* ::Section:: *)
(*Begin Package*)


BeginPackage["MultiHypExp`"]


Print["MultiHypExp v1.0
Author : Souvik Bera (souvikbera@iisc.ac.in)\n\n"];


(* ::Section:: *)
(*Public Commands*)


SeriesExpand::usage=
"This command finds the series expansion of multivariate hypergeometric functions. It takes the input in two different forms.\n
\[FilledSquare] SeriesExpand[ Indices_List, exp, Var_List, e_Symbol, p_Integer]
Finds the first 'p_Integer' coefficients of the series expansion of 'exp' with respect to 'e_Symbol',
where 'Var_List' = list of variables and 'Indices_List' = list of summation indices.
The 'exp' must be series representation of a multivariate hypergeometric function.\n
For instance\n SeriesExpand[{m},\!\(\*FractionBox[\(\(\\\ \)\(Pochhammer[e, m] Pochhammer[e, m] \*SuperscriptBox[\(x\), \(m\)]\)\), \(Pochhammer[1 + e, m] \(m!\)\)]\),{x},e,3]
yields the series expansion of\!\(\*SubscriptBox[\(\\\ \), \(2\)]\)\!\(\*SubscriptBox[\(F\), \(1\)]\)(e,e,1+e,x) up to O(e^3)\n
\[FilledSquare] SeriesExpand[FunctionName,Parameter_List,Var_List,e_Symbol,p_Integer]
Finds the first 'p_Integer' coefficients of 'FunctionName' series with respect to 'e_Symbol',
where 'Var_List' = list of variables and 'Parameter_List' = list of Pochhammer paramters.
For the expansion of one variable\!\(\*SubscriptBox[\(\\\ \), \(p\)]\)\!\(\*SubscriptBox[\(F\), \(p - 1\)]\) functions, the 'FunctionName' input can be omitted
and the 'Parameter_List'  must be given in the following form.
'Parameter_List' = {Upper_Pochhammer_Parameter_List, Lower_Pochhammer_Parameter_List}\n
For instance\n
SeriesExpand[F1, {e,e,e,1+e} , {x,y} , e ,3]
finds the first 3 series coefficients of Appell \!\(\*SubscriptBox[\(F\), \(1\)]\)(e,e,e,1+e,x,y)\n
SeriesExpand[ { {e,e,e} , {1+e,1+e} }  , {x} , e , 4]
yields the series expansion of\!\(\*SubscriptBox[\(\\\ \), \(3\)]\)\!\(\*SubscriptBox[\(F\), \(2\)]\)(e,e,e,1+e,1+e,x) up to O(e^4).

\[FivePointedStar] The available 'FunctionName' are given below.
Double variable series : F1, F2 ,F3, F4, G1, G2, G3, H1, H2, H3, H4, H6 and H7
Three variable series : FA3, FB3, FD3, FK3, FM3, FN3 and FS3

\[Bullet] Note that, the series expansion of Appell \!\(\*SubscriptBox[\(F\), \(4\)]\) and Horn \!\(\*SubscriptBox[\(H\), \(1\)]\) can be found with certain restrictions of its Pochhammer parameters
";



ReduceFunction::usage=
" This command finds the reduction formula of a multivariate hypergeometric functions.
\[FilledSquare] ReduceFunction[FunctionName, Parameter_List, Var_List]
Finds the reduction formula of the function 'FunctionName',
where 'Var_List' = list of variables and 'Parameter_List' = list of Pochhammer paramters, which should be non-negative integers.
For Instance\n
ReduceFunction[F2, {1, 1, 1, 2, 2}, {x, y}]
yields the reduction formula of Appell \!\(\*SubscriptBox[\(F\), \(2\)]\)(1,1,1;2,2;x,y) in term of multiple polylogs.

\[FivePointedStar] The available 'FunctionName' are : F1, F2, F3, F4, FD3, FS3
";



G3replacements::usage=
"The SeriesExpand command yields the series expansion of G3 function in terms of Global`w and Global`z variables.
The transformation rule from Global`w, Global`z to x,y can be found using this command
G3replacements[x,y]";


(*ToCANONICA::usage="";
SolutionVector::usage="";
Solution::usage="";
ToHYPERDIREInput::usage="";
FromHYPERDIREOutput::usage="";

infotoseries::usage="";
info::usage="";
generalser::usage="";
DiffOp::usage="";
Operators::usage="";
PDEGenerator::usage="";
ToPfaffSystem::usage="";
ToPfaffSystem2::usage="";
ShiftedSeries::usage="shifts the series";
SolvePfaffSystem::usage="";
SolutionVectorGeneral::usage="";
ApplyOperator::usage="";
BoundaryCondition::usage="";*)


(* ::Section:: *)
(*Calling other packages*)


(*LaunchKernels[4];*)


(*SetSharedVariable[Global`$PolyLogPath,Global`$HPLPath];*)


(* ::Text:: *)
(*Calling CANONICA in main kernel*)


AppendTo[System`$Path,Global`$CANONICAPath];
AppendTo[System`$Path,Global`$HYPERDIREPath];
Get["AppellF1F4.m"]//Quiet;
Get["HornFunctions.m"]//Quiet;
Get["FdFunction.m"]//Quiet;
Get["FsFunction.m"]//Quiet;
Print["\n"];
Get["CANONICA.m"]//Quiet;
(*Get["RISC/HolonomicFunctions`"]//Quiet;*)
$ComputeParallel=True;
(*SET RISC path*)


(* ::Text:: *)
(*Calling PolyLogTools in all kernels*)


Off[General::shdw];ParallelEvaluate[Off[General::shdw]];
Off[ParallelCombine::mopar1];ParallelEvaluate[Off[ParallelCombine::mopar1]];


(* ::Section:: *)
(*Private` begins*)


Begin["`Private`"]


(* ::Section::Closed:: *)
(*ReduceFunction*)


F1ser[x_,y_]={1,0,-G[0,1,x]-G[0,1,y],G[1,y] G[0,1,x]-G[1,y] G[0,y,x]+G[0,0,1,y]+G[0,1,1,x]+G[0,1,1,y]+G[0,y,1,x],G[0,1,y] G[0,y,x]-G[0,1,x] G[1,1,y]+G[0,y,x] G[1,1,y]-G[1,y] G[0,1,1,x]-G[1,y] G[0,y,1,x]+2 G[1,y] G[0,y,y,x]-G[0,0,0,1,y]-G[0,0,1,1,y]-G[0,1,1,1,x]-G[0,1,1,1,y]+G[0,y,0,1,x]-G[0,y,1,1,x]-2 G[0,y,y,1,x],-G[0,y,x] G[0,0,1,y]+G[1,1,y] G[0,1,1,x]-G[0,y,x] G[0,1,1,y]+G[1,1,y] G[0,y,1,x]-2 G[0,1,y] G[0,y,y,x]-2 G[1,1,y] G[0,y,y,x]+G[0,1,x] G[1,1,1,y]-G[0,y,x] G[1,1,1,y]+G[1,y] G[0,1,1,1,x]-G[1,y] G[0,y,0,1,x]+G[1,y] G[0,y,0,y,x]+G[1,y] G[0,y,1,1,x]+2 G[1,y] G[0,y,y,1,x]-4 G[1,y] G[0,y,y,y,x]+G[0,0,0,0,1,y]+G[0,0,0,1,1,y]+G[0,0,1,1,1,y]+G[0,1,1,1,1,x]+G[0,1,1,1,1,y]-G[0,y,0,1,1,x]-G[0,y,0,y,1,x]+G[0,y,1,1,1,x]-2 G[0,y,y,0,1,x]+2 G[0,y,y,1,1,x]+4 G[0,y,y,y,1,x]};
F2ser[x_,y_]={1,0,-G[0,1,x]-G[0,1,y],G[1,y] G[0,1,x]+G[0,0,1,x]+G[0,0,1,y]+G[0,1,1,y]+G[0,1,1-y,x],-G[0,1,x] G[1,1,y]-G[1,y] G[0,0,1,x]-G[1,y] G[0,1,1-y,x]-G[0,0,0,1,x]-G[0,0,0,1,y]-G[0,0,1,1,y]-G[0,0,1,1-y,x]-G[0,1,1,1,y]-G[0,1,1-y,1-y,x],G[1,1,y] G[0,0,1,x]+G[1,1,y] G[0,1,1-y,x]+G[0,1,x] G[1,1,1,y]+G[1,y] G[0,0,0,1,x]+G[1,y] G[0,0,1,1-y,x]+G[1,y] G[0,1,1-y,1-y,x]+G[0,0,0,0,1,x]+G[0,0,0,0,1,y]+G[0,0,0,1,1,y]+G[0,0,0,1,1-y,x]+G[0,0,1,1,1,y]+G[0,0,1,1-y,1-y,x]+G[0,1,1,1,1,y]+G[0,1,1-y,1-y,1-y,x]};
F3ser[x_,y_]={1,0,-G[0,1,x]-G[0,1,y],G[0,0,1,x]+G[0,0,1,y]+G[0,1,1,x]+G[0,1,1,y],G[0,1,x] G[0,1,y]-G[1,y] G[0,0,y/(-1+y),x]+G[1,y] G[0,1,y/(-1+y),x]-G[0,0,0,1,y]-G[0,0,1,1,x]-G[0,0,1,1,y]-G[0,0,y/(-1+y),1,x]-G[0,1,0,1,x]-G[0,1,1,1,x]-G[0,1,1,1,y]+G[0,1,y/(-1+y),1,x],-G[0,1,y] G[0,0,1,x]-G[0,1,x] G[0,0,1,y]+G[0,1,y] G[0,0,y/(-1+y),x]+G[1,1,y] G[0,0,y/(-1+y),x]-G[0,1,y] G[0,1,1,x]-G[0,1,x] G[0,1,1,y]-G[0,1,y] G[0,1,y/(-1+y),x]-G[1,1,y] G[0,1,y/(-1+y),x]-G[1,y] G[0,0,1,y/(-1+y),x]+3 G[1,y] G[0,0,y/(-1+y),y/(-1+y),x]+2 G[1,y] G[0,1,0,y/(-1+y),x]-G[1,y] G[0,1,1,y/(-1+y),x]-3 G[1,y] G[0,1,y/(-1+y),y/(-1+y),x]+G[0,0,0,0,1,y]+G[0,0,0,1,1,y]+G[0,0,1,0,1,x]+G[0,0,1,1,1,x]+G[0,0,1,1,1,y]-G[0,0,1,y/(-1+y),1,x]-2 G[0,0,y/(-1+y),0,1,x]+G[0,0,y/(-1+y),1,1,x]+3 G[0,0,y/(-1+y),y/(-1+y),1,x]-G[0,1,0,0,1,x]+G[0,1,0,1,1,x]+2 G[0,1,0,y/(-1+y),1,x]+G[0,1,1,0,1,x]+G[0,1,1,1,1,x]+G[0,1,1,1,1,y]-G[0,1,1,y/(-1+y),1,x]+2 G[0,1,y/(-1+y),0,1,x]-G[0,1,y/(-1+y),1,1,x]-3 G[0,1,y/(-1+y),y/(-1+y),1,x]};
(*F4ser[x_,y_]={1/(1-(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)),-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]/(-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]/(-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)),-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2/(2 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y))))-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y))-G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2/(2 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y))),-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^3/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y))))-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(2 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2)/(2 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^3/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y))),-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^4/(24 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y))))-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^3 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2)/(4 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^3)/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^4/(24 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))};*)
FD3ser[x_,y_,z_]={1,0,-G[0,1,x]-G[0,1,y]-G[0,1,z],G[1,y] G[0,1,x]+G[1,z] G[0,1,x]+G[1,z] G[0,1,y]-G[1,y] G[0,y,x]-G[1,z] G[0,z,x]-G[1,z] G[0,z,y]-G[0,0,1,x]+G[0,0,1,z]+G[0,1,1,x]+G[0,1,1,y]+G[0,1,1,z]+G[0,y,1,x]+G[0,z,1,x]+G[0,z,1,y],-G[1,y] G[1,z] G[0,1,x]+G[1,y] G[1,z] G[0,y,x]-G[1,z] G[z,y] G[0,y,x]+G[1,z] G[z,y] G[0,z,x]+G[0,1,y] G[0,z,x]+G[0,1,z] G[0,z,x]+G[0,1,z] G[0,z,y]-G[0,1,x] G[1,1,y]+G[0,y,x] G[1,1,y]-G[0,1,x] G[1,1,z]-G[0,1,y] G[1,1,z]+G[0,z,x] G[1,1,z]+G[0,z,y] G[1,1,z]+G[0,y,x] G[z,1,y]-G[0,z,x] G[z,1,y]+G[1,y] G[0,0,1,x]+G[1,z] G[0,0,1,x]-G[1,y] G[0,0,y,x]-G[1,z] G[0,0,z,x]-G[1,y] G[0,1,1,x]-G[1,z] G[0,1,1,x]-G[1,z] G[0,1,1,y]-G[1,y] G[0,y,1,x]-G[1,z] G[0,y,1,x]+2 G[1,y] G[0,y,y,x]+G[1,z] G[0,y,z,x]-G[1,y] G[0,z,1,x]-G[1,z] G[0,z,1,x]-G[1,z] G[0,z,1,y]+G[1,y] G[0,z,y,x]+2 G[1,z] G[0,z,z,x]+2 G[1,z] G[0,z,z,y]-G[0,0,0,1,x]-G[0,0,0,1,z]+G[0,0,1,1,x]-G[0,0,1,1,z]+G[0,0,y,1,x]+G[0,0,z,1,x]-G[0,1,1,1,x]-G[0,1,1,1,y]-G[0,1,1,1,z]+2 G[0,y,0,1,x]-G[0,y,1,1,x]-2 G[0,y,y,1,x]-G[0,y,z,1,x]+2 G[0,z,0,1,x]+G[0,z,0,1,y]-G[0,z,1,1,x]-G[0,z,1,1,y]-G[0,z,y,1,x]-2 G[0,z,z,1,x]-2 G[0,z,z,1,y]};
FS3ser[x_,y_,z_]={1,0,-G[0,1,x]-G[0,1,y]-G[0,1,z],G[1,z] G[0,1,y]-G[1,z] G[0,z,y]+G[0,0,1,x]+G[0,0,1,z]+G[0,1,1,x]+G[0,1,1,y]+G[0,1,1,z]+G[0,z,1,y],G[0,1,x] G[0,1,y]+G[0,1,x] G[0,1,z]+G[0,1,z] G[0,z,y]-G[0,1,y] G[1,1,z]+G[0,z,y] G[1,1,z]-G[1,y] G[0,0,y/(-1+y),x]-G[1,z] G[0,0,z/(-1+z),x]-G[1,z] G[0,1,1,y]+G[1,y] G[0,1,y/(-1+y),x]+G[1,z] G[0,1,z/(-1+z),x]-G[1,z] G[0,z,1,y]+2 G[1,z] G[0,z,z,y]+G[0,0,0,1,x]-G[0,0,0,1,z]-G[0,0,1,1,x]-G[0,0,1,1,z]-G[0,0,y/(-1+y),1,x]-G[0,0,z/(-1+z),1,x]-2 G[0,1,0,1,x]-G[0,1,1,1,x]-G[0,1,1,1,y]-G[0,1,1,1,z]+G[0,1,y/(-1+y),1,x]+G[0,1,z/(-1+z),1,x]+G[0,z,0,1,y]-G[0,z,1,1,y]-2 G[0,z,z,1,y]};
F4ser[x_,y_]={(1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2])/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)]),((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)])/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(2 (-1+x+y) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2]),((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2)/(12 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2)/(12 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(2 (-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-(2 (-1+x+y) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+(5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(8 (-1+x+y) G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+((1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]),((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^3)/(36 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(12 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2)/(12 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^3)/(36 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-((-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-(2 (-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-((-1+x+y) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(2 G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(2 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+(5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(4 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+1/Sqrt[x^2+(-1+y)^2-2 x (1+y)] (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(8 (-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-(8 (-1+x+y) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(8 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(4 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 2 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-(5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/Sqrt[x^2+(-1+y)^2-2 x (1+y)]+(20 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(8 G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(20 (-1+x+y) G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(32 (-1+x+y) G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(4 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]),((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^4)/(144 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^3 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(36 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2)/(24 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^3)/(36 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y+5 Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^4)/(144 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^3 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(18 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^3 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(18 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-((-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^3 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(9 Sqrt[-4 x y+(-1+x+y)^2])-((-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-((-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-((-1+x+y) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^3 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(9 Sqrt[-4 x y+(-1+x+y)^2])+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(2 G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+1/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+(5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(4 G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(4 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+1/Sqrt[x^2+(-1+y)^2-2 x (1+y)] (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+1/Sqrt[x^2+(-1+y)^2-2 x (1+y)] (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(6 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(2 G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(2 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(4 G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))-(4 (-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-(8 (-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-(4 (-1+x+y) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(8 G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(8 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(8 G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(4 G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(4 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(4 G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y))+(16 G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+1/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)]^2 G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(6 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]^2 G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 2 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 2 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 8 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/Sqrt[x^2+(-1+y)^2-2 x (1+y)] (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-(5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/Sqrt[x^2+(-1+y)^2-2 x (1+y)]-(5 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/Sqrt[x^2+(-1+y)^2-2 x (1+y)]+(4 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y))+1/Sqrt[x^2+(-1+y)^2-2 x (1+y)] 3 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(20 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(20 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(16 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))-1/Sqrt[x^2+(-1+y)^2-2 x (1+y)] 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(8 G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(8 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(16 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(20 (-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(20 (-1+x+y) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(20 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 10 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(16 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,0,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(2 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(32 (-1+x+y) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])-(32 (-1+x+y) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(32 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 16 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)]+(4 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(4 G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(4 G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(8 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(4 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 8 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 8 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 2 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1+x-y-Sqrt[-4 x y+(-1+x+y)^2])/(2 x)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,(1-x+y-Sqrt[-4 x y+(-1+x+y)^2])/(2 y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 2 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)] G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+(15 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/Sqrt[x^2+(-1+y)^2-2 x (1+y)]-(20 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/Sqrt[x^2+(-1+y)^2-2 x (1+y)]-(50 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(80 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[0,1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(40 G[1,0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(20 (-1+x+y) G[1,0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/Sqrt[-4 x y+(-1+x+y)^2]+(8 G[1,0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(80 (-1+x+y) G[1,0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(16 G[1,0,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(8 G[1,1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(80 (-1+x+y) G[1,1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+G[1,1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(128 (-1+x+y) G[1,1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)])/(3 Sqrt[-4 x y+(-1+x+y)^2])+(4 G[1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(16 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3-(3 (-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2)/(4 x y))+(4 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))+(4 G[1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 (-1+(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y)))-(20 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+(4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 8 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+(8 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])+((1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)])/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)])-1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) 4 (-1+x+y+Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),0,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]+1/(3 Sqrt[x^2+(-1+y)^2-2 x (1+y)]) (1-x-y-Sqrt[(-1+x)^2-2 (1+x) y+y^2]) G[(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])^2/(4 x y),1,(-1+x+y+Sqrt[-4 x y+(-1+x+y)^2])/(2 y)]};


(*ToSeriesData[exp_,e_Symbol,p0_Integer]:=SeriesData[e,0,{exp},0,p0,1];
ToSeriesData[exp_,e_Symbol,p0_Integer]/;Head[exp]===SeriesData:=exp;*)
ReduceFunction[FunctionName_,ParameterList_List,Variables_List]/;And@@(#>0&/@ParameterList):=Module[
{changeparameterlist,operators,hyperdireop, formula, opapplied,x,y,z,result,e,temp,ind,N0,
a,b,c,d},
If[And@@(IntegerQ/@ParameterList),Nothing[],Print["Parameters must be integers"];Abort[];];

(*consistency conditions*)
If[Or[ToString[FunctionName]==="F1"&&Length[ParameterList]===4,
ToString[FunctionName]==="F2"&&Length[ParameterList]===5,
ToString[FunctionName]==="F3"&&Length[ParameterList]===5,
ToString[FunctionName]==="F4"&&Length[ParameterList]===4,
ToString[FunctionName]==="FD3"&&Length[ParameterList]===5,
ToString[FunctionName]==="FS3"&&Length[ParameterList]===6
],Nothing[],Print["wrong parameter list"];Abort[];];

Which[(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~   F1    ~~~~~~~~~~~~~~~~~~~~~~~*)
ToString[FunctionName]==="F1",
changeparameterlist = {0,0,0,1}-ParameterList;
hyperdireop = FullSimplify/@First[F1IndexChange[changeparameterlist,Join[ParameterList+e,{x,y}]]];

ind=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(hyperdireop));
	N0=Abs[Min[ind[[All,4]]]];

operators = {#&, x PolyLogTools`DG[#,x]&,y PolyLogTools`DG[#,y]&};
formula=F1ser[x,y][[1;;N0+1]]/.MultiHypExp`Private`G-> PolyLogTools`G;,
(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~   F2    ~~~~~~~~~~~~~~~~~~~~~~~*)
ToString[FunctionName]==="F2",
changeparameterlist = {0,0,0,1,1}-ParameterList;
hyperdireop = FullSimplify/@First[F2IndexChange[changeparameterlist,Join[ParameterList+e,{x,y}]]];
ind=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(hyperdireop));
	N0=Abs[Min[ind[[All,4]]]];
operators = {#&, x PolyLogTools`DG[#,x]&,y PolyLogTools`DG[#,y]&, x y PolyLogTools`DG[PolyLogTools`DG[#,x],y]&};
formula=F2ser[x,y][[1;;N0+1]]/.MultiHypExp`Private`G-> PolyLogTools`G;,
(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~   F3    ~~~~~~~~~~~~~~~~~~~~~~~*)
ToString[FunctionName]==="F3",
changeparameterlist = {0,0,0,0,1}-ParameterList;
hyperdireop = FullSimplify/@First[F3IndexChange[changeparameterlist,Join[ParameterList+e,{x,y}]]];
ind=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(hyperdireop));
	N0=Abs[Min[ind[[All,4]]]];
operators = {#&, x PolyLogTools`DG[#,x]&,y PolyLogTools`DG[#,y]&, x y PolyLogTools`DG[PolyLogTools`DG[#,x],y]&};
formula=F3ser[x,y][[1;;N0+1]]/.MultiHypExp`Private`G-> PolyLogTools`G;,
(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~   F4   ~~~~~~~~~~~~~~~~~~~~~~~*)
ToString[FunctionName]==="F4",
changeparameterlist = {0,0,0,0}-ParameterList;
(*Print[changeparameterlist];*)
hyperdireop =Simplify/@First[F4IndexChange[changeparameterlist,Join[{a,b,c,d},{x,y}]]]/.MapThread[Rule,{{a,b,c,d},ParameterList+{e,2e,3e,2e}}];
Print["hyperdire done"];
ind=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(hyperdireop));
(*Print[ind];*)
N0=Abs[Min[ind[[All,4]]]];
operators = {#&, x PolyLogTools`DG[#,x]&,y PolyLogTools`DG[#,y]&, x y PolyLogTools`DG[PolyLogTools`DG[#,x],y]&};
formula=F4ser[x,y][[1;;N0+1]]/.MultiHypExp`Private`G-> PolyLogTools`G;,
(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~  FD3   ~~~~~~~~~~~~~~~~~~~~~~~*)
ToString[FunctionName]==="FD3"&&Length[Variables]===3,
changeparameterlist = {0,0,0,0,1}-ParameterList;
hyperdireop = FullSimplify/@First[FdIndexChange[{changeparameterlist[[1]],changeparameterlist[[2;;-2]],changeparameterlist[[-1]]},Join[{ParameterList[[1]],ParameterList[[2;;-2]],ParameterList[[-1]]}+e,{{x,y,z}}]]];
ind=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(hyperdireop));
	N0=Abs[Min[ind[[All,4]]]];
operators = {#&, x PolyLogTools`DG[#,x]&,y PolyLogTools`DG[#,y]&,z PolyLogTools`DG[#,z]&};
formula=FD3ser[x,y,z][[1;;N0+1]]/.MultiHypExp`Private`G-> PolyLogTools`G;,
(*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~  FS3   ~~~~~~~~~~~~~~~~~~~~~~~*)
ToString[FunctionName]==="FS3"&&Length[Variables]===3,
changeparameterlist = {0,0,0,0,0,1}-ParameterList;
hyperdireop = Simplify/@First[FsIndexChange[changeparameterlist,Join[ParameterList+e,{x,y,z}]]];
ind=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(hyperdireop));
	N0=Abs[Min[ind[[All,4]]]];
operators = {#&, x PolyLogTools`DG[#,x]&,y PolyLogTools`DG[#,y]&,z PolyLogTools`DG[#,z]&, x y PolyLogTools`DG[PolyLogTools`DG[#,x],y]&, x z PolyLogTools`DG[PolyLogTools`DG[#,x],z]&};
formula=FS3ser[x,y,z][[1;;N0+1]]/.MultiHypExp`Private`G-> PolyLogTools`G;


];

(*common part*)
DistributeDefinitions[x,y,z,e];
SetSharedVariable[operators,formula,hyperdireop,opapplied];
opapplied = ParallelTable[If[hyperdireop[[i]]===0,0,hyperdireop[[i]](operators[[i]]@(Table[e^(i-1),{i,1,Length[formula]}]formula))],{i,1,Length[operators]}];
result =Table[PolyLogTools`GCoefficientSimplify[Plus@@opapplied[[i]]],{i,1,Length[operators]}];
result = Normal[Series[Total[result],{e,0,0}]];
Return[result/.MapThread[Rule,{If[Length[Variables]===2,{x,y},{x,y,z}],Variables}]]]


(* ::Section::Closed:: *)
(*Common commands 	*)


(* ::Text:: *)
(*info :: This command converts the information of the parameters in a series into a list form.*)
(*infotoseries :: This command does the opposite work of the command info. It takes the o/p format of info command and convert into a series*)
(*generalser :: given a info_list it finds the general expression of a series*)


infotoseries[indices_List,infolist_(*info o/p*)]:=Module[{pre,serinfo,var,num,denom},
{pre,{serinfo,var}}=infolist; 
(*Print[serinfo];*)
num=Times@@(Pochhammer@@{#[[1]],(*Plus@@*)(#[[2]] . indices)}&/@serinfo[[1]]);
(*Print[num];*)
denom=Times@@(Pochhammer@@{#[[1]],(*Plus@@*)(#[[2]] . indices)}&/@serinfo[[2]]);
Return[pre num Times@@Power[var,indices]/(denom Times@@(Factorial[#]&/@indices) )];
];


info[indices_List,exp_]:=Module[{prefactor,poch,list11},
prefactor= exp/.((#-> 0)&/@indices);

poch=(#/(#/.Table[indices[[i]]-> 0,{i,1,Length[indices]}]))&/@List@@If[Head[exp]=!=Plus,Flatten[{exp}],Expand[exp]];
(*poch=poch/.Factorial[m_]\[Rule] Pochhammer[1,m];*);
(*Print[{"info",poch}];*)
list11=Table[{{SortBy[Flatten[List@@Which[Numerator[#]===1,{},Head[Numerator[#]]=!=Times,{Numerator[#]},True,Numerator[#]]/.Power[x_,n_]:> Table[x,{i,1,n}]]/.Pochhammer[z_,m_]:>{z, Plus@@(CoefficientRules[m,indices]/.Rule-> Times)},Last],
			SortBy[Flatten[List@@Which[Denominator[#]===1,{},Head[Denominator[#]]=!=Times,{Denominator[#]},True,Denominator[#]]/.Power[x_,n_]:> Table[x,{i,1,n}]]/.Pochhammer[z_,m_]:>{z, Plus@@(CoefficientRules[m,indices]/.Rule-> Times)},Last]}&@(poch[[i]]/(poch[[i]]/.Pochhammer[z_,m_]-> 1)),Table[((#/.Pochhammer[z_,m_]-> 1/.Table[If[k===j,indices[[k]]-> 1,indices[[k]]-> 0],{k,1,Length[indices]}]))&@poch[[i]],
{j,1,Length[indices]}]}
,{i,1,Length[poch]}];
(*Print[list11];*)
Return[Prepend[list11,prefactor]];
];


generalser[indices_,infolist_]:=Module[{pre,num,denom,num1,denom1,var},
{pre,{{num,denom},var}}=infolist;
denom=Select[denom,#[[1]]=!=1&];
num1=Table[{{ToExpression[ToString[Global`aa]<>ToString[i]],num[[All,2]][[i]]},ToExpression[ToString[Global`aa]<>ToString[i]]-> num[[All,1]][[i]]},{i,1,Length[num[[All,2]]]}];
denom1=Table[{{ToExpression[ToString[Global`bb]<>ToString[i]],denom[[All,2]][[i]]},ToExpression[ToString[Global`bb]<>ToString[i]]-> denom[[All,1]][[i]]},{i,1,Length[denom[[All,2]]]}];
Return[{infotoseries[indices,{pre,{{num1[[All,1]],denom1[[All,1]]},var}}],Join[num1[[All,2]],denom1[[All,2]]]}];];


(* ::Section::Closed:: *)
(*ShiftedSeries*)


(* ::Text:: *)
(*ShiftedSeries : Given a series representation, it finds the associated series that can be expanded in Taylor series*)


parameterlist[exp_,eps_Symbol]:={Coefficient[exp,eps,0],Coefficient[exp,eps,1]};
shiftparameter[exp_List,eps_Symbol]:=Module[{},Return[Total[{1,eps}Which[exp[[2]]===0,exp,exp[[1]]<=0,{1,exp[[2]]},True,exp]]]];
shiftparameter2[exp_List,eps_Symbol]:=Module[{},Return[Total[{1,eps}Which[exp[[2]]===0,exp,exp[[1]]>0,{0,exp[[2]]},True,exp]]]];
shiftedser[indices_List,exp_,e_Symbol]:=Module[{infoser,num,denom},
infoser  = info[indices,exp];
num={shiftparameter2[parameterlist[#[[1]],e],e],#[[2]]}&/@infoser[[2,1,1]];
denom={shiftparameter[parameterlist[#[[1]],e],e],#[[2]]}&/@infoser[[2,1,2]];
Return[{infoser[[1]],{{num,denom},infoser[[2,-1]]}}];
];


ShiftedSeries[indices_, exp_, e_Symbol]:=Module[{serinfo},
  serinfo = shiftedser[indices, exp, e];
  Return[infotoseries[indices, serinfo]];
  ];


(* ::Section::Closed:: *)
(*Operators	*)


(* ::Text:: *)
(*Operators : It finds the differential operators that relates the two series.*)


DiffOp[indices_,exp_,e_]:=Module[{num,denom,seriesinfo,th,fromnum,fromdenom,all},
(*{num,denom}=Through[{Numerator,Denominator}[exp]]/.Pochhammer[a_,m_]\[Rule] Gamma[a+m]/.Factorial[m_]\[Rule] Gamma[m+1];*)
th=Table[Subscript[Global`t, i],{i,1,Length[indices]}];
seriesinfo=info[indices,exp][[2,1]];
(*Print[seriesinfo];*)
num={seriesinfo[[1,All,1]],seriesinfo[[1,All,2]]};
denom={seriesinfo[[2,All,1]],seriesinfo[[2,All,2]]};
fromnum=Flatten[Table[If[FreeQ[num[[1,i]],e],Nothing,Table[(num[[2,i]] . th)/(num[[1,i]]-j)+ 1,{j,1,(num[[1,i]]-e)/.e-> 0}]],{i,1,Length[num[[2]]]}]];
fromdenom=Flatten[Table[If[FreeQ[denom[[1,i]],e],Nothing,Table[(denom[[2,i]] . th)/(denom[[1,i]]+j-1)+ 1,{j,1,(e+1-denom[[1,i]])/.e-> 0}]],{i,1,Length[denom[[2]]]}]];
all=Join[fromnum,fromdenom];
Return[If[all==={},1,NonCommutativeMultiply@@all]];
];


Operators[indices_List,exp_,variables_List,e_Symbol]:=Module[{infoser,newinfoser,pde,newser,genser,op,op2,gb,theta,reducedop,result,temp1},
infoser = info[indices,exp];
newinfoser = shiftedser[indices,exp,e];
genser=generalser[indices,newinfoser];
(*Print[genser];*)
pde= PDEGenerator[indices,variables,genser[[1]]]/.genser[[2]];

gb=RISC`HolonomicFunctions`OreGroebnerBasis[pde,Last[Flatten[{RISC`HolonomicFunctions`OreAlgebra[pde]}]]];
(*Print[RISC`HolonomicFunctions`UnderTheStaircase[gb]];*)

theta[x_] := x RISC`HolonomicFunctions`Der[x];
op=DiffOp[indices,exp,e];
(*Print[op];*)
(*If[op ===1, Return[op]];*)
op2 = RISC`HolonomicFunctions`ToOrePolynomial[op/.Table[Subscript[Global`t, i]-> theta[variables[[i]]],{i,1,Length[indices]}],
RISC`HolonomicFunctions`OreAlgebra@@Join[variables,RISC`HolonomicFunctions`Der/@variables]];
(*Print[op2];*)
reducedop=RISC`HolonomicFunctions`OreReduce[op2,gb(*,RISC`HolonomicFunctions`OreAlgebra@@Join[variables,RISC`HolonomicFunctions`Der/@variables]*)];
result=Collect[reducedop,RISC`HolonomicFunctions`UnderTheStaircase[gb],Simplify];
temp1=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(result[[1,All,1]]));
(*Print[temp1];*)
Print["Differential operator has pole in e of order : ", Abs[Min[temp1[[All,4]]]]];

(*Print["The result will be correct upto ",SeriesData[e,0,List[],1,6- Abs[Min[temp1[[All,4]]]],1]];*)

Return[{Abs[Min[temp1[[All,4]]]],result}];];


(* ::Section::Closed:: *)
(*ToPfaff system	*)


(* ::Text:: *)
(*PDEGenerator ::  Given a series representation, it finds the PDE associated with it*)
(*ToPfaffSystem :: Given a series representation, it finds the Pfaff system associated with it*)


PDEGenerator[indices_List,var_List,series_]:=Module[{serinfo,replacements,AA,AA1,AA2,AA3,AA4,th},

(*<<RISC`HolonomicFunctions`  is needed *)

(*AA is the ratio of A_{m+1,n}/A_{m,n}...etc*)
(*serinfo=info[indices,series0];
{series,replacements}=generalser[indices,serinfo];*)
AA=Table[Through[{Numerator,Denominator}[FullSimplify[(series/.Table[var[[i]]-> 1,{i,1,Length[var]}]/.indices[[j]]-> indices[[j]]+1)/(series/.Table[var[[i]]-> 1,{i,1,Length[var]}])]]],{j,1,Length[indices]}];
(*th is Euler theta operators*)
(*Print[AA//Factor];*)
th=Table[var[[i]] RISC`HolonomicFunctions`Der[var[[i]]],{i,1,Length[indices]}];
AA1=AA/.Table[indices[[i]]->  th[[i]],{i,1,Length[indices]}];
AA2=Table[NonCommutativeMultiply@@@AA1[[i]],{i,1,Length[AA1]}];
AA3=Table[{RISC`HolonomicFunctions`ToOrePolynomial[AA2[[i,1]]],RISC`HolonomicFunctions`ToOrePolynomial[AA2[[i,2]]**(1 /var[[i]])]},{i,1,Length[AA2]}];
(*AA4 are the PDEs*)
AA4=Table[AA3[[i,2]]-AA3[[i,1]],{i,1,Length[AA3]}];
(*Print[AA4];*)

Return[AA4];];


ToPfaffSystem[indices_List,variables_List,exp_](*/;Length[variables]===2*):=Module[
{x,y,z,A,A1,vec,genser,replacements,infoser,pde,gb,mat,fff,mat2,elements,basis,mat10,mat11,mat12},
infoser=info[indices,exp];
(*Print[infoser];*)
{genser,replacements}=generalser[indices,infoser];
(*Print[{genser,replacements}];*)
pde = PDEGenerator[indices,variables,genser];
(*Print["pde"];*)
gb = RISC`HolonomicFunctions`OreGroebnerBasis[pde,Last[Flatten[{RISC`HolonomicFunctions`OreAlgebra[pde]}]]];
(*Print["gb",Length[gb]];
Print[gb];*)
elements=RISC`HolonomicFunctions`UnderTheStaircase[gb];
If[elements===Infinity,Print["Holonomic Rank is infinity"];Abort[];];

basis=Table[Times@@Power[variables,elements[[i]][[1,1,2]]]** elements[[i]],{i,1,Length[elements]}];
mat = Table[RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@(Collect[#,RISC`HolonomicFunctions`UnderTheStaircase[gb],Simplify]&/@(RISC`HolonomicFunctions`OreReduce[#,gb]&/@Thread[RISC`HolonomicFunctions`Der[variables[[i]]]**basis])),{i,1,Length[variables]}]
/.MapThread[Rule,{RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@elements,(RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@elements)/Table[Times@@Power[variables,elements[[i]][[1,1,2]]],{i,1,Length[elements]}]}];
(*mat = Table[ApplyOreOperator[#,fff@@variables]&/@(Collect[#,UnderTheStaircase[gb],Simplify]&/@(OreReduce[#,gb]&/@Thread[Der[variables[[i]]]**basis])),{i,1,Length[variables]}];*)
mat2 =Table[ Table[Apart//@Factor//@Simplify[Coefficient[Together[mat[[j]][[iii]]],RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@RISC`HolonomicFunctions`UnderTheStaircase[gb]]],{iii,1,Length[mat[[j]]]}],{j,1,Length[mat]}];
Print["basis vector is ",basis];
Return[{indices,exp,mat2/.replacements,basis}];


	
];


ToPfaffSystem2[indices_List,variables_List,exp_](*/;Length[variables]===2*):=Module[{genser,replacements,infoser,pde,gb,mat,fff,mat2,elements,basis,basis0},
infoser=info[indices,exp];
{genser,replacements}=generalser[indices,infoser];
pde = PDEGenerator[indices,variables,genser];
gb = RISC`HolonomicFunctions`OreGroebnerBasis[pde,Last[Flatten[{RISC`HolonomicFunctions`OreAlgebra[pde]}]]]/.replacements;
elements=RISC`HolonomicFunctions`UnderTheStaircase[gb];
If[elements===Infinity,Print["Holonomic Rank is infinity"];Abort[];];
basis=Table[Times@@Power[variables,elements[[i]][[1,1,2]]]** elements[[i]],{i,1,Length[elements]}];
Print["basis vector is ",basis];
mat = Table[RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@(Collect[#,RISC`HolonomicFunctions`UnderTheStaircase[gb],Simplify]&/@(RISC`HolonomicFunctions`OreReduce[#,gb]&/@Thread[RISC`HolonomicFunctions`Der[variables[[i]]]**basis])),{i,1,Length[variables]}]
/.MapThread[Rule,{RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@elements,(RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@elements)/Table[Times@@Power[variables,elements[[i]][[1,1,2]]],{i,1,Length[elements]}]}];
(*mat = Table[ApplyOreOperator[#,fff@@variables]&/@(Collect[#,UnderTheStaircase[gb],Simplify]&/@(OreReduce[#,gb]&/@Thread[Der[variables[[i]]]**basis])),{i,1,Length[variables]}];*)
mat2 =Table[ Table[Apart//@Factor//@Simplify[Coefficient[Together[mat[[j]][[iii]]],RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@RISC`HolonomicFunctions`UnderTheStaircase[gb]]],{iii,1,Length[mat[[j]]]}],{j,1,Length[mat]}];
Return[{indices,exp,mat2,basis}];

	
];


ToPfaffSystem3[indices_List,variables_List,exp_](*/;Length[variables]===2*):=Module[
{x,y,z,A,A1,vec,genser,replacements,infoser,pde,gb,mat,fff,mat2,elements,basis,mat10,mat11,mat12},
infoser=info[indices,exp];
{genser,replacements}=generalser[indices,infoser];
pde = PDEGenerator[indices,variables,genser];
Print["pde is generated"];
gb = RISC`HolonomicFunctions`OreGroebnerBasis[pde,Last[Flatten[{RISC`HolonomicFunctions`OreAlgebra[pde]}]]];
Print["gb is calculated"];
elements=RISC`HolonomicFunctions`UnderTheStaircase[gb];
If[elements===Infinity,Print["Holonomic Rank is infinity"];Abort[];];

If[Length[variables]===3&&Length[elements]<= 8,
{x,y,z}=variables;
vec={1,x,y,z,x y,y z,z x ,x y z}{1, RISC`HolonomicFunctions`Der[x],RISC`HolonomicFunctions`Der[y],RISC`HolonomicFunctions`Der[z],
RISC`HolonomicFunctions`Der[y] RISC`HolonomicFunctions`Der[x], RISC`HolonomicFunctions`Der[y] RISC`HolonomicFunctions`Der[z],
RISC`HolonomicFunctions`Der[z] RISC`HolonomicFunctions`Der[x],
RISC`HolonomicFunctions`Der[y] RISC`HolonomicFunctions`Der[x] RISC`HolonomicFunctions`Der[z]};



basis=RISC`HolonomicFunctions`ToOrePolynomial[vec[[1;;Length[elements]]]];
A = RISC`HolonomicFunctions`ApplyOreOperator[Collect[RISC`HolonomicFunctions`OreReduce[RISC`HolonomicFunctions`ToOrePolynomial[basis],gb],elements,Simplify],fff@@variables];
(*Print["A"];*)
A1= Table[Apart//@Factor//@FullSimplify[Coefficient[Together[A[[iii]]],RISC`HolonomicFunctions`ApplyOreOperator[elements,fff[x,y,z]]]],{iii,1,Length[A]}];];

If[Det[A1]===0,

basis=Table[Times@@Power[variables,elements[[i]][[1,1,2]]]** elements[[i]],{i,1,Length[elements]}];
mat = Table[RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@(Collect[#,RISC`HolonomicFunctions`UnderTheStaircase[gb],Simplify]&/@(RISC`HolonomicFunctions`OreReduce[#,gb]&/@Thread[RISC`HolonomicFunctions`Der[variables[[i]]]**basis])),{i,1,Length[variables]}]
/.MapThread[Rule,{RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@elements,(RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@elements)/Table[Times@@Power[variables,elements[[i]][[1,1,2]]],{i,1,Length[elements]}]}];
(*mat = Table[ApplyOreOperator[#,fff@@variables]&/@(Collect[#,UnderTheStaircase[gb],Simplify]&/@(OreReduce[#,gb]&/@Thread[Der[variables[[i]]]**basis])),{i,1,Length[variables]}];*)
mat2 =Table[ Table[Apart//@Factor//@Simplify[Coefficient[Together[mat[[j]][[iii]]],RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@RISC`HolonomicFunctions`UnderTheStaircase[gb]]],{iii,1,Length[mat[[j]]]}],{j,1,Length[mat]}];
Print["basis vector is ",basis];
Return[mat2/.replacements];

,
basis=Table[Times@@Power[variables,elements[[i]][[1,1,2]]]** elements[[i]],{i,1,Length[elements]}];
Print["basis vector is ",basis];
(*mat = Table[RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@(Collect[#,elements,Simplify]&/@(RISC`HolonomicFunctions`OreReduce[#,gb]&/@Thread[RISC`HolonomicFunctions`Der[variables[[i]]]**basis])),{i,1,Length[variables]}];*)

mat10= Table[Thread[RISC`HolonomicFunctions`Der[variables[[i]]]**basis],{i,1,Length[variables]}];
(*Print[1];
Print["mat10"];*)
SetSharedVariable[gb,mat10];
(*mat11 = RISC`HolonomicFunctions`OreReduce[#,gb]&/@mat10;*)
mat11 = ParallelMap[RISC`HolonomicFunctions`OreReduce[#,gb]&,mat10];
(*Print[2];
Print["mat11"];*)
mat12=Table[Collect[#,elements,Simplify]&/@(mat11[[i]]),{i,1,Length[mat11]}];
(*Print[3];
Print["mat12"];*)
mat = RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@mat12;
(*Print["mat"];*)
mat2 =Table[ Table[Apart//@Factor//@Simplify[Coefficient[Together[mat[[j]][[iii]]],RISC`HolonomicFunctions`ApplyOreOperator[#,fff@@variables]&/@elements]],{iii,1,Length[mat[[j]]]}],{j,1,Length[mat]}];
(*Print["mat2"];*)
(*Print[mat2];*)
SetSharedVariable[mat2,A1];
Return[ParallelMap[Simplify,mat2/.replacements]];

];




	
];


(* ::Section::Closed:: *)
(*SolvePfaffSystem*)


(* ::Text:: *)
(*SolvePfaffSystem :: given the Pfaff system, it finds the solution of it.*)
(*SolutionVectorGeneral :: There are three version of this command depending on the number of variables*)


(*BoundaryCondition[indices_List,exp_,Tmatrix_,basis_,variables_,e_Symbol]:=Module[{m0=5,sum,temp,temp2,temp3},
sum=Sum@@Join[{exp},Table[{indices[[i]],0,m0},{i,1,Length[indices]}]];
temp=(Inverse[Tmatrix].Table[RISC`HolonomicFunctions`ApplyOreOperator[basis[[i]],sum],{i,1,Length[basis]}]);
SetSharedVariable[temp];
temp2 = ParallelMap[Simplify,temp];
temp3=Limit[temp2,MapThread[Rule,{variables,Table[0,{i,1,Length[variables]}]}]];
Return[Simplify/@temp3];
];*)


BoundaryCondition[indices_List,exp_,Tmatrix_,basis_,variables_,e_Symbol]:=Module[{m0=5,sum,temp,temp2,temp3},
(*Print[FullForm[basis]];*)
sum=Sum@@Join[{exp},Table[{indices[[i]],0,m0},{i,1,Length[indices]}]];
(*Print[sum];*)
temp=(Inverse[Tmatrix] . Table[RISC`HolonomicFunctions`ApplyOreOperator[basis[[i]],sum],{i,1,Length[basis]}]);
SetSharedVariable[temp];
temp2 = ParallelMap[Simplify,temp];
(*Print[temp2];*)
(*temp3=temp2/.MapThread[Rule,{variables,Table[0,{i,1,Length[variables]}]}];*)
temp3=Limit[#,MapThread[Rule,{variables,Table[0,{i,1,Length[variables]}]}]]&/@temp2;
Return[Simplify/@temp3];
];


Solution[exp1_List(*ToCANONICA output*),exp2_List(*SolutionVector output*),e_Symbol,p0_Integer(*o/p is shown upto e^p0*),m0_Integer(*internal calculations are done upto e^m0*)]:=Module[
{result2, lowestorder, result1},
{lowestorder,result1}=exp2;
result2 = Plus@@( Table[e^i,{i,lowestorder,Length[result1]+lowestorder-1}]result1);
Return[ToSeriesData[#,e,p0]&@(Series[First[exp1[[1]] . result2],{e,0,p0}])];

];


SolvePfaffSystem[Pfaff0_(*ToPfaffSystem o/p*),variables_List,e_Symbol,p0_Integer]:=Module[{indices,tocanonica,sol1,sol2,m0=p0,Pfaff,basis,exp,BC},
{indices,exp,Pfaff,basis}=Pfaff0;
tocanonica=Quiet[CANONICA`TransformDiagonalBlock[Pfaff/.e->CANONICA`eps,variables]]/.CANONICA`eps-> e;
Print["CANONICA results are obtained"];
If[tocanonica===False,Abort[]];

BC=BoundaryCondition[indices,exp,tocanonica[[1]],basis,variables,e];
(*Print[BC];*)
(*Print[{tocanonica,BC,variables,e,m0}];*)
sol1=SolutionVectorGeneral[tocanonica,BC,variables,e,m0];
sol2 = Solution[tocanonica,sol1,e,m0,m0];
Return[sol2]];


(*For 1 variable*)
SolutionVectorGeneral[exp_List(*ToCANONICA output*),BC0_, variables_List, e_Symbol, p0_Integer]/;Length[variables]===1:=Module[{m0,BC,de,x,y,Omega,result1,i1,i2,i3,i4,i5,i6,i7,i8,i9,lowestorder},
{x}=variables;
SetSharedVariable[x];
Omega = exp[[2]];

(*BC = ToSeriesData[#,e,p0+1]&/@(Series[#,{e,0,p0}]&/@Limit[Simplify[Inverse[exp[[1]]].Join[{1},Table[0,{i,1,Length[exp[[2,1]]]-1}]]],MapThread[Rule,{variables,{0,0}}]]);			
*)
BC=ToSeriesData[#,e,p0+1]&/@(Series[#,{e,0,p0}]&/@BC0);	
(*Print[BC];*)
lowestorder=Min[BC[[All,4]]];
(*m0=Abs[lowestorder];*)
m0=p0;
(*Print[m0];*)
(*Print[lowestorder];*)
(*setting the boundary condition*)
i4=SeriesCoefficient[#,lowestorder]&/@BC;
(*i9 =If[ToString[FunctionName]==="F2",SeriesCoefficient[#,lowestorder]&/@BC];*)
(*i9 =If[ToString[FunctionName]==="F2",If[#==={},0,#[[1]]]&/@((ToSeriesData[#,e,p0]&/@BC)[[All,3]])(*SeriesCoefficient[#,0]&/@BC*)];*)
(*result1 =If[ToString[FunctionName]==="F2",{i9}];*)
result1={i4};
(*Print[i4];*)


(* i1 to i9 are steps for solving the DE*)
Table[
PrintTemporary["Solving "<>ToString[i]<>" -th order"];
de =(# . i4)&/@Omega/.e-> 1;
SetSharedVariable[de];
i1 = ParallelMap[PolyLogTools`GCoefficientSimplify,de[[1]]];
SetSharedVariable[i1];
i2=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{x,y}]&,If[Head[i1[[i]]]===Plus,List@@(i1[[i]]),{i1[[i]]}]]],{i,1,Length[i1]}];
SetSharedVariable[i2];
i3 = ParallelMap[PolyLogTools`GIntegrate[#,x]&,i2];
SetSharedVariable[i3];
i4=  ParallelMap[PolyLogTools`GCoefficientSimplify,i3];
SetSharedVariable[i4];
i4=i4+(SeriesCoefficient[#,lowestorder+i]&/@BC);
result1 = AppendTo[result1,i4];,
{i,1,m0}];
Return[{lowestorder,result1/.SeriesCoefficient[a_,n_Integer]/;FreeQ[a,e]&&n>=1-> 0}];
];


(*For 2 variables*)
SolutionVectorGeneral[exp_List(*ToCANONICA output*),BC0_, variables_List, e_Symbol, p0_Integer]/;Length[variables]===2:=Module[{m0,BC,de,x,y,Omega,result1,i1,i2,i3,i4,i5,i6,i7,i8,i9,lowestorder},
{x,y}=variables;
SetSharedVariable[x,y];
Omega = exp[[2]];

(*BC = ToSeriesData[#,e,p0+1]&/@(Series[#,{e,0,p0}]&/@Limit[Simplify[Inverse[exp[[1]]].Join[{1},Table[0,{i,1,Length[exp[[2,1]]]-1}]]],MapThread[Rule,{variables,{0,0}}]]);			
*)
BC=ToSeriesData[#,e,p0+1]&/@(Series[#,{e,0,p0}]&/@BC0);	
(*Print[BC];*)
lowestorder=Min[BC[[All,4]]];
(*m0=Abs[lowestorder];
Print[m0];*)
m0=p0;	
(*Print[lowestorder];*)
(*setting the boundary condition*)
i9=SeriesCoefficient[#,lowestorder]&/@BC;
(*i9 =If[ToString[FunctionName]==="F2",SeriesCoefficient[#,lowestorder]&/@BC];*)
(*i9 =If[ToString[FunctionName]==="F2",If[#==={},0,#[[1]]]&/@((ToSeriesData[#,e,p0]&/@BC)[[All,3]])(*SeriesCoefficient[#,0]&/@BC*)];*)
(*result1 =If[ToString[FunctionName]==="F2",{i9}];*)
result1={i9};
(*Print[i9];*)


(* i1 to i9 are steps for solving the DE*)
Table[
PrintTemporary["Solving "<>ToString[i]<>" -th order"];
de =(# . i9)&/@Omega/.e-> 1;
SetSharedVariable[de];
i1 = ParallelMap[PolyLogTools`GCoefficientSimplify,de[[1]]];
SetSharedVariable[i1];
i2=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{x,y}]&,If[Head[i1[[i]]]===Plus,List@@(i1[[i]]),{i1[[i]]}]]],{i,1,Length[i1]}];
SetSharedVariable[i2];
i3 = ParallelMap[PolyLogTools`GIntegrate[#,x]&,i2];
SetSharedVariable[i3];
i4=  ParallelMap[PolyLogTools`GCoefficientSimplify,i3];
SetSharedVariable[i4];
(*Print[i3];*)
i5 =de[[2]]- ParallelMap[PolyLogTools`DG[#,y]&,i4];
SetSharedVariable[i5];
i6 = ParallelMap[PolyLogTools`GCoefficientFullSimplify,i5];
SetSharedVariable[i6];
i7=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{y,x}]&,If[Head[i6[[i]]]===Plus,List@@(i6[[i]]),{i6[[i]]}]]],{i,1,Length[i6]}];
SetSharedVariable[i7];
i8 = ParallelMap[PolyLogTools`GIntegrate[#,y]&,i7];
SetSharedVariable[i8];
i9= ParallelMap[PolyLogTools`GCoefficientFullSimplify,(i4+i8)];
SetSharedVariable[i9];
i9=i9+(SeriesCoefficient[#,lowestorder+i]&/@BC);
result1 = AppendTo[result1,i9];,
{i,1,m0}];
Return[{lowestorder,result1/.SeriesCoefficient[a_,n_Integer]/;FreeQ[a,e]&&n>=1-> 0}];
];


(*For three variables*)
SolutionVectorGeneral[exp_List(*ToCANONICA output*),BC0_, variables_List, e_Symbol, p0_Integer]/;Length[variables]===3:=Module[{m0,de,x,y,z,Omega,result1,i1,i2,i3,i4,i5,i6,
i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,BC,lowestorder},
{x,y,z}=variables;
SetSharedVariable[x,y,z];
Omega = exp[[2]];

(*BC = ToSeriesData[#,e,p0+1]&/@(Series[#,{e,0,p0}]&/@Which[
ToString[FunctionName]==="Fd3",Normal[Series[#,{x,0,0}]]&/@(Normal[Series[#,{y,0,0}]]&/@(Normal[Series[#,{z,0,0}]]&/@(Simplify[Inverse[exp[[1]]] . {1,0,0,0}])))]);	*)		
(*Print[Normal[Series[#,{x,0,0}]]&/@(Normal[Series[#,{y,0,0}]]&/@(Normal[Series[#,{z,0,0}]]&/@(Simplify[Inverse[exp[[1]]] . {1,0,0,0}])))];	*)
(*Print[BC];*)
BC=ToSeriesData[#,e,p0+1]&/@(Series[#,{e,0,p0}]&/@BC0);
lowestorder=Min[BC[[All,4]]];
(*m0=Abs[lowestorder];*)
m0=p0;	
(*Print[lowestorder]*);
(*setting the boundary condition*)
i16=SeriesCoefficient[#,lowestorder]&/@BC;
(*i9 =If[ToString[FunctionName]==="F2",SeriesCoefficient[#,lowestorder]&/@BC];*)
(*i9 =If[ToString[FunctionName]==="F2",If[#==={},0,#[[1]]]&/@((ToSeriesData[#,e,p0]&/@BC)[[All,3]])(*SeriesCoefficient[#,0]&/@BC*)];*)
(*result1 =If[ToString[FunctionName]==="F2",{i9}];*)
result1={i16};
(*Print[i16];*)



(* i1 to i9 are steps for solving the DE*)
Table[
PrintTemporary["Solving "<>ToString[i]<>" -th order"];
de =(# . i16)&/@Omega/.e-> 1;
SetSharedVariable[de];
i1 = ParallelMap[PolyLogTools`GCoefficientSimplify,de[[1]]];
SetSharedVariable[i1];
i2=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{x,y,z}]&,If[Head[i1[[i]]]===Plus,List@@(i1[[i]]),{i1[[i]]}]]],{i,1,Length[i1]}];
SetSharedVariable[i2];
(*Print[i2];*)
i3 = ParallelMap[PolyLogTools`GIntegrate[#,x]&,i2];
SetSharedVariable[i3];
i4=  ParallelMap[PolyLogTools`GCoefficientSimplify,i3];
SetSharedVariable[i4];
(*Print[i3];*)
i5 =de[[2]]- ParallelMap[PolyLogTools`DG[#,y]&,i4];
SetSharedVariable[i5];
i6 = ParallelMap[PolyLogTools`GCoefficientFullSimplify,i5];
SetSharedVariable[i6];
i7=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{y,z,x}]&,If[Head[i6[[i]]]===Plus,List@@(i6[[i]]),{i6[[i]]}]]],{i,1,Length[i6]}];
SetSharedVariable[i7];
i8 = ParallelMap[PolyLogTools`GIntegrate[#,y]&,i7];
SetSharedVariable[i8];
i9= ParallelMap[PolyLogTools`GCoefficientFullSimplify,i8];
SetSharedVariable[i9];
i10 =  ParallelMap[PolyLogTools`DG[#,z]&,i4];
SetSharedVariable[i10];
i11 =  ParallelMap[PolyLogTools`DG[#,z]&,i9];
SetSharedVariable[i11];
i12 =de[[3]]-i10-i11;
SetSharedVariable[i12];
i13 = ParallelMap[PolyLogTools`GCoefficientFullSimplify,i12];
SetSharedVariable[i13];
i14=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{z,x,y}]&,If[Head[i13[[i]]]===Plus,List@@(i13[[i]]),{i13[[i]]}]]],{i,1,Length[i13]}];
SetSharedVariable[i14];
i15 = ParallelMap[PolyLogTools`GIntegrate[#,z]&,i14];
SetSharedVariable[i15];
i16=  ParallelMap[PolyLogTools`GCoefficientSimplify,i15+i4+i9];
i16=i16+(SeriesCoefficient[#,lowestorder+i]&/@BC);
result1 = AppendTo[result1,i16];,
{i,1,m0}];
Return[{lowestorder,result1/.SeriesCoefficient[a_,n_Integer]/;FreeQ[a,e]&&n>=1-> 0}];
];


(* ::Section::Closed:: *)
(*ApplyOperators*)


(* ::Text:: *)
(*ApplyOperator :: apply the differential operator on the series expansion of the secondary series*)


ToPolyLogOps[op_(*Operators o/p*)]:={op[[1,All,1]],(Times@@Power[op[[2,1]],#]&/@(op[[1,All,2]])) /. 1-> DGOp /.Power[RISC`HolonomicFunctions`Der[x_],n_.]:>DGOp[x,n]/.Times->Composition }


DGOp[exp_]:=exp;
DGOp[x_Symbol,n_Integer][exp_]:= Nest[PolyLogTools`DG[#,x]&,exp,n];
DistributeDefinitions[DGOp];


ApplyOperator[op0_(*ToPolyLogOps o/p*),exp_(*SolvePfaffSystem o/p*),e_Symbol,p0_Integer]:=Module[{temp,optemp,temp2,temp3,op,N0},
{N0,op}=op0;
temp = Table[e^i,{i,0,Length[exp[[3]]]-1}]exp[[3]];
optemp=ToPolyLogOps[op];
SetSharedVariable[optemp];
temp2=optemp[[1]]Table[ParallelMap[optemp[[2,i]],temp],{i,1,Length[optemp[[2]]]}];
SetSharedVariable[temp2];
(*Print[temp2];*)
temp3 = ParallelMap[PolyLogTools`GCoefficientSimplify,Total[Flatten[temp2]]];
Return[ToSeriesData[#,e,p0]&@Series[temp3,{e,0,p0-N0}]];
]


(* ::Section:: *)
(*GeneralSeriesExpand*)


SeriesExpand[{{a_,b_},{c_}},{x_},e_Symbol,p_Integer]/;p>=1&&p<=6&&MemberQ[{a,b},c]:=If[b===c,Series[(1-x)^-a,{e,0,p}],Series[(1-x)^-b,{e,0,p}]];


SeriesExpand[{ParameterList1_List,ParameterList2_List},variables_List,e_Symbol,p_Integer]/;And[Length[variables]===1,p>=1,p<=6]:=Module[
{m},

Return[SeriesExpand[{m},(Times @@ (Pochhammer[#, m] & /@ ParameterList1) variables[[1]]^m)/(Times @@ (Pochhammer[#, m] & /@ ParameterList2) m!),variables,e,p]]
]


SeriesExpand[indices_List,exp_,variables_,e_Symbol,p_Integer]/;p>=1&&p<=6:=Module[
{step1,step2,step3,step4,step5,p0,N0,ind},

(*p0=If[p>6,6,p];*)
step1 = ShiftedSeries[indices,exp,e];
PrintTemporary["Shifted Series is obtained"];
step2 = ToPfaffSystem[indices,variables,step1];
PrintTemporary["Pfaff system is obtained"];
(*Print[step2];*)
step4 = Operators[indices,exp,variables,e];
PrintTemporary["Operators are obtained"];
(*Print[step4];*)
(*Print[step4[[2]]];*)
ind=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(step4[[2,1,All,1]]));
(*Print[ind];
Print[Abs[Min[#[[-3]]&/@ind]]];*)
(*N0 is same as lowest order*)
(*N0=Abs[Min[ind[[All,4]]]];
p0=N0+p-1;*)

N0=Abs[Min[#[[-3]]&/@ind]];
p0=N0+p;
(*If[p0===0,p0=1];*)

(*Print[{N0,p0}];*)


(*N0=Abs[Min[#[[-3]]&/@ind]];
p0=N0+p-1;
If[p0===0,p0=1];
If[p0>6,p0=6];*)
(*Print[p0];*)
(*Print[{N0,p0}];*)
(*Here make changes*)


step3 = SolvePfaffSystem[step2,variables,e,p0];
PrintTemporary["Solution of the Pfaff system is obtained"];
(*Print[step3];*)
step5 = ApplyOperator[step4,step3,e,p0];
(*Print[step5];*)
step5 = ToSeriesData[#,e,p]&@ReplacePart[step5,-2-> step5[[-3]]+ p];
(*Print[step5//FullForm];*)
PrintTemporary["Operators are applied"];
SetSharedVariable[step5];
Return[ReplacePart[step5,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,step5[[3]]]]];

]



SeriesExpand[FunctionName_,ParameterList_List,variables_List,e_Symbol,p_Integer]/;And[Length[variables]===3,p>=1,p<=6,
Or[ToString[FunctionName]==="FD3",ToString[FunctionName]==="FN3",ToString[FunctionName]==="FA3",ToString[FunctionName]==="FB3",
ToString[FunctionName]==="FK3",ToString[FunctionName]==="FS3",ToString[FunctionName]==="FM3"]]:=Module[
{step1,step2,step3,step4,step5,p0,mm,nn,pp,indices,exp,x,y,z,vec,infoser,genser,
replacements,gb,pde,elements,newparameter,ind,N0},
(*p0=If[p>6,6,p];*)
indices = {mm,nn,pp};
{x,y,z}=variables;
vec={1,x,y,z,x y,y z,z x ,x y z}{1, RISC`HolonomicFunctions`Der[x],RISC`HolonomicFunctions`Der[y],RISC`HolonomicFunctions`Der[z],
RISC`HolonomicFunctions`Der[y] RISC`HolonomicFunctions`Der[x], RISC`HolonomicFunctions`Der[y] RISC`HolonomicFunctions`Der[z],
RISC`HolonomicFunctions`Der[z] RISC`HolonomicFunctions`Der[x],
RISC`HolonomicFunctions`Der[y] RISC`HolonomicFunctions`Der[x] RISC`HolonomicFunctions`Der[z]};

Which[
ToString[FunctionName]==="FD3",If[Length[ParameterList]=!=5,Print["Length of parameter list should be 5"];Abort[];];exp = FD3ser[ParameterList,variables,indices],
ToString[FunctionName]==="FN3",If[Length[ParameterList]=!=7,Print["Length of parameter list should be 7"];Abort[];];exp = FN3ser[ParameterList,variables,indices],
ToString[FunctionName]==="FA3",If[Length[ParameterList]=!=7,Print["Length of parameter list should be 7"];Abort[];];exp = FA3ser[ParameterList,variables,indices],
ToString[FunctionName]==="FB3",If[Length[ParameterList]=!=7,Print["Length of parameter list should be 7"];Abort[];];exp = FB3ser[ParameterList,variables,indices],
ToString[FunctionName]==="FK3",If[Length[ParameterList]=!=7,Print["Length of parameter list should be 7"];Abort[];];exp = FK3ser[ParameterList,variables,indices],
ToString[FunctionName]==="FS3",If[Length[ParameterList]=!=6,Print["Length of parameter list should be 6"];Abort[];];exp = FS3ser[ParameterList,variables,indices],
ToString[FunctionName]==="FM3",If[Length[ParameterList]=!=6,Print["Length of parameter list should be 6"];Abort[];];exp = FM3ser[ParameterList,variables,indices]];
infoser=info[indices,exp];

{genser,replacements}=generalser[indices,infoser];
pde = PDEGenerator[indices,variables,genser];
gb = RISC`HolonomicFunctions`OreGroebnerBasis[pde,Last[Flatten[{RISC`HolonomicFunctions`OreAlgebra[pde]}]]];
elements=RISC`HolonomicFunctions`UnderTheStaircase[gb];

step1 = ShiftedSeries[indices,exp,e];
PrintTemporary["Shifted Series is obtained"];
newparameter=(Join@@(info[indices,step1][[2,1]]))[[All,1]];

step2 = {indices,step1,Which[
ToString[FunctionName]==="FD3",
{mat1FD3[newparameter,variables],mat2FD3[newparameter,variables],mat3FD3[newparameter,variables]}
,ToString[FunctionName]==="FN3",
{mat1FN3[newparameter,variables],mat2FN3[newparameter,variables],mat3FN3[newparameter,variables]} 
,ToString[FunctionName]==="FA3",
{mat1FA3[newparameter,variables],mat2FA3[newparameter,variables],mat3FA3[newparameter,variables]}
,ToString[FunctionName]==="FB3",
{mat1FB3[newparameter,variables],mat2FB3[newparameter,variables],mat3FB3[newparameter,variables]} 
,ToString[FunctionName]==="FS3",
{mat1FS3[newparameter,variables],mat2FS3[newparameter,variables],mat3FS3[newparameter,variables]} 
,ToString[FunctionName]==="FK3",
{mat1FK3[newparameter,variables],mat2FK3[newparameter,variables],mat3FK3[newparameter,variables]}
,ToString[FunctionName]==="FM3",
{mat1FM3[newparameter,variables],mat2FM3[newparameter,variables],mat3FM3[newparameter,variables]} ],
RISC`HolonomicFunctions`ToOrePolynomial[vec[[1;;Length[elements]]]]};


(*step2 = ToPfaffSystem[indices,variables,step1];*)
PrintTemporary["Pfaff system is obtained"];
step4 = Operators[indices,exp,variables,e];



ind=ToSeriesData[#,e,1]&/@(Series[#,{e,0,0}]&/@(step4[[2,1,All,1]]));
(*Print[ind];*)
N0=Abs[Min[#[[-3]]&/@ind]];
p0=N0+p-1;
If[p0===0,p0=1];
(*Return[step4];*)

(*ind=ToSeriesData[#,eps,1]&/@(Series[#,{eps,0,0}]&/@(step2[[1]]));
(*Print[ind[[All,4]]];*)

(*N0 is same as lowest order*)
N0=Abs[Min[ind[[All,4]]]];
m0=N0+p0-1;*)

PrintTemporary["Operators are obtained"];
step3 = SolvePfaffSystem[step2,variables,e,p0];
(*Return[step3];*)


PrintTemporary["Solution of the Pfaff system is obtained"];
step5 = ApplyOperator[step4,step3,e,p0];
step5 = ToSeriesData[#,e,p]&@ReplacePart[step5,-2->  step5[[-3]]+ p];
PrintTemporary["Operators are applied"];
SetSharedVariable[step5];
Return[ReplacePart[step5,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,step5[[3]]]]];

]


(* ::Section:: *)
(*Main commands for two variables*)


ToCANONICA[FunctionName_,ParameterList_List,variables_List,e_Symbol]:=Module[{Mat, var,canonicaresult},
Mat = Which[
ToString[FunctionName]==="F1",{mat1F1[ParameterList,variables],mat2F1[ParameterList,variables]},
ToString[FunctionName]==="F2",{mat1F2[ParameterList,variables],mat2F2[ParameterList,variables]},
ToString[FunctionName]==="F3",{mat1F3[ParameterList,variables],mat2F3[ParameterList,variables]},
ToString[FunctionName]==="H2",{mat1H2[ParameterList,variables],mat2H2[ParameterList,variables]},
ToString[FunctionName]==="FD3",{mat1FD3[ParameterList,variables],mat2FD3[ParameterList,variables],mat3FD3[ParameterList,variables]},
True,Return[]]/.e-> CANONICA`eps;
var = Variables[variables];

(*Print[Mat];*)
(*If[CheckIntegrability[Mat,var], Print[True];,Print[False];];*)

canonicaresult = Quiet[CANONICA`TransformDiagonalBlock[Mat, var]]/.CANONICA`eps-> e;

Return[canonicaresult];
];


ToSeriesData[exp_,e_Symbol,p0_Integer]:=SeriesData[e,0,{exp},0,p0,1];
ToSeriesData[exp_,e_Symbol,p0_Integer]/;Head[exp]===SeriesData:=exp;


(*For 2 variables*)
SolutionVector[FunctionName_, exp_List(*ToCANONICA output*), variables_List, e_Symbol, p0_Integer]/;Length[variables]===2:=Module[{de,x,y,Omega,result1,i1,i2,i3,i4,i5,i6,i7,i8,i9,BC,lowestorder},
{x,y}=variables;
SetSharedVariable[x,y];
Omega = exp[[2]];

BC = ToSeriesData[#,e,p0+1]&/@(Series[#,{e,0,p0}]&/@Which[Or[ToString[FunctionName]==="F2",ToString[FunctionName]==="F3",ToString[FunctionName]==="H2"],
Limit[Simplify[Inverse[exp[[1]]] . {1,0,0,0}],MapThread[Rule,{variables,{0,0}}]],
ToString[FunctionName]==="F1",Limit[Simplify[Inverse[exp[[1]]] . {1,0,0}],MapThread[Rule,{variables,{0,0}}]]]);			
	
(*Print[BC];*)
lowestorder=Min[BC[[All,4]]];
	
(*Print[lowestorder];*)
(*setting the boundary condition*)
i9=SeriesCoefficient[#,lowestorder]&/@BC;
(*i9 =If[ToString[FunctionName]==="F2",SeriesCoefficient[#,lowestorder]&/@BC];*)
(*i9 =If[ToString[FunctionName]==="F2",If[#==={},0,#[[1]]]&/@((ToSeriesData[#,e,p0]&/@BC)[[All,3]])(*SeriesCoefficient[#,0]&/@BC*)];*)
(*result1 =If[ToString[FunctionName]==="F2",{i9}];*)
result1={i9};
(*Print[i9];*)


(* i1 to i9 are steps for solving the DE*)
Table[
de =(# . i9)&/@Omega/.e-> 1;
SetSharedVariable[de];
i1 = ParallelMap[PolyLogTools`GCoefficientSimplify,de[[1]]];
SetSharedVariable[i1];
i2=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{x,y}]&,If[Head[i1[[i]]]===Plus,List@@(i1[[i]]),{i1[[i]]}]]],{i,1,Length[i1]}];
SetSharedVariable[i2];
i3 = ParallelMap[PolyLogTools`GIntegrate[#,x]&,i2];
SetSharedVariable[i3];
i4=  ParallelMap[PolyLogTools`GCoefficientSimplify,i3];
SetSharedVariable[i4];
(*Print[i3];*)
i5 =de[[2]]- ParallelMap[PolyLogTools`DG[#,y]&,i4];
SetSharedVariable[i5];
i6 = ParallelMap[PolyLogTools`GCoefficientFullSimplify,i5];
SetSharedVariable[i6];
i7=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{y,x}]&,If[Head[i6[[i]]]===Plus,List@@(i6[[i]]),{i6[[i]]}]]],{i,1,Length[i6]}];
SetSharedVariable[i7];
i8 = ParallelMap[PolyLogTools`GIntegrate[#,y]&,i7];
SetSharedVariable[i8];
i9= ParallelMap[PolyLogTools`GCoefficientFullSimplify,(i4+i8)];
SetSharedVariable[i9];
i9=i9+(SeriesCoefficient[#,lowestorder+i]&/@BC);
result1 = AppendTo[result1,i9];,
{i,1,p0}];
Return[{lowestorder,result1/.SeriesCoefficient[a_,n_Integer]/;FreeQ[a,e]&&n>=1-> 0}];
];


(*Modified for Fd3*)
(*"shiftparameter" is changed a bit*)
parameterlist[exp_,eps_Symbol]:={Coefficient[exp,eps,0],Coefficient[exp,eps,1]};
shiftparameter[exp_List,eps_Symbol]:=Module[{},Return[Total[{1,eps}Which[exp[[2]]===0,exp,exp[[1]]<=0,{1,exp[[2]]},True,exp]]]];
shiftparameter2[exp_List,eps_Symbol]:=Module[{},Return[Total[{1,eps}Which[exp[[2]]===0,exp,exp[[1]]>0,{0,exp[[2]]},True,exp]]]];
ToHYPERDIREInput[FunctionName_,ParameterList_List,variables_List,eps_Symbol]:=Module[{shiftlist, noshiftlist,newpara1,newpara2,newpara3,HYPintegerlist},

(*Command works for F1,F2,F3,Fd3, H2*)
Which[ToString[FunctionName]==="F1", If[Length[ParameterList]===4,Nothing[];,Print["wrong parameterlist"];Abort[];];,
Or[ToString[FunctionName]==="F2",ToString[FunctionName]==="F3",ToString[FunctionName]==="H2"], If[Length[ParameterList]===5,Nothing[];,Print["wrong parameterlist"];Abort[];];,
ToString[FunctionName]==="Fd3", If[And[Length[ParameterList]===5,Length[variables]===3],Nothing[];,Print["wrong parameterlist"];Abort[];];,
True,Abort[];];

(*shifting the necessary parameters for the above functions*)
If[ToString[FunctionName]==="H2",
{noshiftlist,shiftlist}={ParameterList[[2;;4]],ParameterList[[{1,-1}]]};
newpara3=Join[ {shiftparameter2[ parameterlist[shiftlist[[1]],eps],eps]}, noshiftlist,{shiftparameter[ parameterlist[shiftlist[[-1]],eps],eps]}];
,

{noshiftlist,shiftlist}=Which[ToString[FunctionName]==="F2",{ParameterList[[1;;-3]],ParameterList[[-2;;-1]]},
Or[ToString[FunctionName]==="F1",ToString[FunctionName]==="F3",ToString[FunctionName]==="Fd3"],{ParameterList[[1;;-2]],ParameterList[[-1;;-1]]}
];
newpara1 = parameterlist[#,eps]&/@shiftlist;
newpara2 = shiftparameter[#,eps]&/@newpara1;
newpara3=Join[noshiftlist,newpara2] ;

];


(*Print[newpara3];*)
(*Return[];*)


Return[
Which[
(*For F1,F2,F3*)
Or[ToString[FunctionName]==="F2",ToString[FunctionName]==="F3",ToString[FunctionName]==="F1"],
{ToExpression["AppellF1F4`"<>ToString[FunctionName]<>"IndexChange"],
newpara3-ParameterList,
Join[ParameterList,variables]},
(*For Fd3*)
ToString[FunctionName]==="Fd3",
{ToExpression["Fd`FdIndexChange"],Join[#[[1;;1]],{#[[2;;-2]]},#[[-1;;-1]]]&@(newpara3-ParameterList),Join[Join[#[[1;;1]],{#[[2;;-2]]},#[[-1;;-1]]]&@ParameterList,{variables}]},
ToString[FunctionName]==="H2",
{ToExpression["HornFunctions`"<>ToString[FunctionName]<>"IndexChange"],
newpara3-ParameterList,
Join[ParameterList,variables]}

](*End which*)



](*end Return*)
]


(*When the IndexChangeVector is {0,0,...,0,0}, then the series is of Taylor Type and no need of the operators. Any vector is fine*)
FromHYPERDIREOutput[exp_]:=Module[{},
Return[If[DeleteDuplicates[Flatten[exp[[2]]]]==={0},{Join[{1},Table[0,{i,1,Length[Flatten[exp[[2]]]]-1}]],exp[[-1]]},Together[exp[[1]]@@(exp[[2;;3]])]]]
]


(*FromHYPERDIREOutput[exp_]:=Module[{},
Return[Together[exp[[1]]@@(exp[[2;;3]])]]
]*)


(* ::Section:: *)
(*Commands for two variables*)


ClearAll[Global`w,Global`z];


G3G1soln[x1_,y1_]={-(1/(12 x1 (-1+x1 y1)))(4+8 x1-12 x1^2 y1+(2 I 2^(1/3) (I+Sqrt[3]) (-1-4 x1+x1^2 (-1+6 y1)))/(-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)+2^(2/3) (1+I Sqrt[3]) (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)),(-18 I 2^(2/3) (-I+Sqrt[3]) x1^5 y1^2+2 2^(2/3) (3 I+Sqrt[3]) x1 Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))]+2^(1/3) Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))] (3 I 2^(1/3)+2^(1/3) Sqrt[3]-3 I (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)+Sqrt[3] (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3))+x1^2 (1+2 y1) (-2^(2/3)-I 2^(2/3) Sqrt[3]-4 (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(2/3)+(-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)-I Sqrt[3] (-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3))+2 x1^3 (-2 2^(2/3)-2 I 2^(2/3) Sqrt[3]-4 2^(2/3) y1-4 I 2^(2/3) Sqrt[3] y1+6 y1^2 (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(2/3)+(-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)-I Sqrt[3] (-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)+2 y1 (-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)-2 I Sqrt[3] y1 (-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3))+x1^4 (2 2^(2/3) (1+I Sqrt[3])+16 2^(2/3) (1+I Sqrt[3]) y1+3 2^(1/3) y1^2 (5 2^(1/3)+5 I 2^(1/3) Sqrt[3]-3 (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)+3 I Sqrt[3] (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3))))/(12 x1^2 y1 (-1+x1 y1) (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(2/3))};
G3G1factor[x1_,y1_,g_,g1_]=(1+1/(12 x1 (-1+x1 y1)) (4+8 x1-12 x1^2 y1+(2 I 2^(1/3) (I+Sqrt[3]) (-1-4 x1+x1^2 (-1+6 y1)))/(-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)+2^(2/3) (1+I Sqrt[3]) (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)))^g1 (1-(-18 I 2^(2/3) (-I+Sqrt[3]) x1^5 y1^2+2 2^(2/3) (3 I+Sqrt[3]) x1 Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))]+2^(1/3) Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))] (3 I 2^(1/3)+2^(1/3) Sqrt[3]-3 I (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)+Sqrt[3] (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3))+x1^2 (1+2 y1) (-2^(2/3)-I 2^(2/3) Sqrt[3]-4 (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(2/3)+(-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)-I Sqrt[3] (-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3))+2 x1^3 (-2 2^(2/3)-2 I 2^(2/3) Sqrt[3]-4 2^(2/3) y1-4 I 2^(2/3) Sqrt[3] y1+6 y1^2 (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(2/3)+(-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)-I Sqrt[3] (-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)+2 y1 (-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)-2 I Sqrt[3] y1 (-4-24 x1-54 x1^4 y1^2+6 x1^2 (-5+6 y1)+x1^3 (4+72 y1)+6 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3))+x1^4 (2 2^(2/3) (1+I Sqrt[3])+16 2^(2/3) (1+I Sqrt[3]) y1+3 2^(1/3) y1^2 (5 2^(1/3)+5 I 2^(1/3) Sqrt[3]-3 (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3)+3 I Sqrt[3] (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(1/3))))/(12 x1^2 y1 (-1+x1 y1) (-2-12 x1-27 x1^4 y1^2+3 x1^2 (-5+6 y1)+x1^3 (2+36 y1)+3 Sqrt[3] Sqrt[x1^4 (-1+x1 y1)^2 (-1-4 y1+27 x1^2 y1^2-2 x1 (2+9 y1))])^(2/3)))^g;
G3replacements[x_,y_]=MapThread[Rule,{{Global`w,Global`z},G3G1soln[x,y]}];


(*For 2 variables*)
SeriesExpand[FunctionName_,ParameterList_List,Variables_List,eps_Symbol,p0_Integer]/;And[Length[Variables]===2,p0>=1,p0<=6]:=Module[
{m0,uu,vv,operators,step1,step2,step3,step4,step5,step6,step7,step8,step9,step10,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,
Rtof2,Rtof21,Rtof22,Rtof2pre,temp,temp1,a0,b0,c0,d0,a10,a20,b10,b20,x0,y0,xx,yy,ind,N0,NN0,mm},

(*m0==If[p0>6,6,p0];*)
m0=p0;
If[MemberQ[ParameterList],eps===False,Return[]];

			
	
Which[
ToString[FunctionName]==="F4",
If[Length[ParameterList]=!=4,Print["wrong parameterlist"];Abort[];];
{a0,b0,c0,d0}=ParameterList;
{x0,y0}=Variables;
Which[
(*case1 F4(a0,b0,c0,b0) = ... F1(...)*)
d0===b0,
Rtof2pre = ((1 - (-1 + x0 + y0 + Sqrt[-4 x0 y0 + (-1 + x0 + y0)^2])/(
   2 x0))^a0) ((1 - (-1 + x0 + y0 + 
    Sqrt[-4 x0 y0 + (-1 + x0 + y0)^2])/(2 y0))^a0) ;
Rtof2 = SeriesExpand["F1",Simplify/@{a0, -b0 + c0, 1 + a0 - c0, c0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{(-1 + x0 + y0 + Sqrt[-4 x0 y0 + (-1 + x0 + y0)^2])/(2 y0), (-1 + x0 + 
  y0 + Sqrt[-4 x0 y0 + (-1 + x0 + y0)^2])^2/(4 x0 y0)}}];,
  
(* (*case4 F4(a0,b0,c0,b0) = ... H3(...)*)
c0===a0,
Rtof2pre = (1-x0-y0)^-b0 ;
Rtof2 = SeriesExpand["H3",Simplify/@{b0,d0-a0,d0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{(x0 y0)/(x0+y0-1)^2,x0/(x0+y0-1)}}];,*)

(*case2 F4(a0,b0,c0,a-c0+1) = ... F2(...)*)
d0===a0-c0+1,
Rtof2pre = 1;
Rtof2 = SeriesExpand["F2",Simplify/@{a0, b0, b0, c0, 1 + a0 - c0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{1/2 (1 + x0 - y0 - Sqrt[-4 x0 + (-1 - x0 + y0)^2]), 1/2 (1 - x0 + y0 -
    Sqrt[-4 x0 + (-1 - x0 + y0)^2])}}];,
 
 (*case3 F4(a0,b0,c0,a0+b0-c0+1) = 2F1(...) 2F1(...)*)
d0===a0+b0-c0+1,

(*Rtof21 = SeriesExpand[{mm},(Pochhammer[a0,mm]Pochhammer[b0,mm]xx^mm)/(Pochhammer[c0,mm]mm!),{xx},eps,p0];*)
Rtof21 = SeriesExpand[{{a0,b0},{c0}},{xx},eps,p0];

(*Rtof22 = SeriesExpand[{mm},(Pochhammer[a0,mm]Pochhammer[b0,mm]yy^mm)/(Pochhammer[1 + a0 + b0 - c0,mm]mm!),{yy},eps,p0];*)
Rtof22 = SeriesExpand[{{a0,b0},{1 + a0 + b0 - c0}},{yy},eps,p0];
temp= Series[Rtof21 Rtof22,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{1/2 (1 + x0 - y0 - Sqrt[-4 x0 + (-1 - x0 + y0)^2]),
1/2 (1 - x0 + y0 - Sqrt[-4 x0 + (-1 - x0 + y0)^2])}}];
,
(*symmetric relations*)
Or[c0===a0,c0===b0-d0+1,c0===a0+b0-d0+1],
temp=SeriesExpand["F4",Simplify/@{b0,a0,d0,c0},{y0,x0},eps,m0];
Return[temp];,
 
True,
Print["No connection formula of F_4 found"];Abort[];
];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];

Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]]];,


ToString[FunctionName]==="G1",
If[Length[ParameterList]=!=3,Print["wrong parameterlist"];Abort[];];
{a0,b0,c0}=ParameterList;
{x0,y0}=Variables;
Rtof2pre = (1/(1+x0+y0))^a0 ;
Rtof2 = SeriesExpand["F2",{1-b0-c0,a0,a0,1-b0,1-c0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{(1+2 x0-Sqrt[1-4 x0 y0])/(2 (1+x0+y0)),(1+2 y0-Sqrt[1-4 x0 y0])/(2 (1+x0+y0))}}];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];
Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]]];,

ToString[FunctionName]==="G2",
If[Length[ParameterList]=!=4,Print["wrong parameterlist"];Abort[];];
{a10,a20,b10,b20}=ParameterList;
{x0,y0}=Variables;
Rtof2pre = (1+x0)^-a10 (1+y0)^-a20 ;
Rtof2 = SeriesExpand["F2",{1-b10-b20,a10,a20,1-b10,1-b20},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{x0/(1+x0),y0/(1+y0)}}];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];
Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]]];,


ToString[FunctionName]==="G3",
If[Length[ParameterList]=!=2,Print["wrong parameterlist"];Abort[];];
{a10,a20}=ParameterList;
{x0,y0}=Variables;
Rtof2pre = G3G1factor[x0,y0,a10,a20] ;
Rtof2 = SeriesExpand["G1",{a10+a20,a10,a20},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{Global`w,Global`z}}];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];
Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]//Quiet]];,


ToString[FunctionName]==="H1",
If[Length[ParameterList]=!=4,Print["wrong parameterlist"];Abort[];];
{a0,b0,c0,d0}=ParameterList;
{x0,y0}=Variables;
Which[
(*case1 H1(d0-c0,b0,c0,d0) = ... H2(...)*)
d0-c0===a0,
Rtof2pre = (2^-b0) ((-((-1 + y0 + Sqrt[ 1 + 2 y0 - 4 x0 y0 + y0^2])/((-1 + x0) y0)))^b0) ;
Rtof2 = SeriesExpand["H2",Simplify/@{-c0 + d0, b0, b0, c0, d0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{(1 + y0 - Sqrt[1 + 2 y0 - 4 x0 y0 + y0^2])/(2 y0), (1 - y0 - Sqrt[
 1 + 2 y0 - 4 x0 y0 + y0^2])/(2 (-1 + x0))}}];,

(*case2 H1(a0,b0,c0,1/2(a0+b0+1)) = ... F2(...)*)
d0===Simplify[1/2 (a0+b0+1)]&&EvenQ[a0+b0/.eps-> 0],
Rtof2pre = (1 - 2 x0 - 2 Sqrt[(-1 + x0) x0])^( 1/2 (-a0 + b0)) (1 - 2 x0 + 2 Sqrt[(-1 + x0) x0])^(1/2 (-a0 - b0)) ;
Rtof2 = SeriesExpand["F2",Simplify/@{b0, (a0 + b0)/2, c0, a0 + b0, 1 - a0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{(4 Sqrt[(-1 + x0) x0])/(Sqrt[(1 - 2 x0)^2] + 
 2 Sqrt[(-1 + x0) x0]), -((1 - 2 x0 - 2 Sqrt[(-1 + x0) x0]) y0)}}];,
 
 
True,
Print["No connection formula of H_1 found"];Abort[];
];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];
Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]]];,



ToString[FunctionName]==="H3",
If[Length[ParameterList]=!=3,Print["wrong parameterlist"];Abort[];];
{a0,b0,c0}=ParameterList;
{x0,y0}=Variables;
Rtof2pre = (1-4 x0)^(-a0/2) ((1+Sqrt[1-4 x0]-4 x0)/(2-8 x0))^(1-c0) ;
Rtof2 = SeriesExpand["F3",{-a0+c0,a0,1+a0-c0,b0,c0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{(-1+Sqrt[1-4 x0]+4 x0)/(-2+8 x0),(y0-Sqrt[1-4 x0] y0)/(2 x0)}}];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];
Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]]];,


ToString[FunctionName]==="H4",
If[Length[ParameterList]=!=4,Print["wrong parameterlist"];Abort[];];
{a0,b0,c0,d0}=ParameterList;
If[Denominator[c0/.eps-> 0]=!=2,Print["The parameter c in H4[a,b,c,d,x,y] must have the form p*eps + q/2\np,q are nonzero-integers "];Abort[]];
{x0,y0}=Variables;
Rtof2pre =(1-(2 Sqrt[x0])/(1+2 Sqrt[x0]))^a0  ;
Rtof2 = SeriesExpand["F2",Simplify/@{a0,-(1/2)+c0,b0,-1+2 c0,d0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{(4 Sqrt[x0])/(1+2 Sqrt[x0]),y0/(1+2 Sqrt[x0])}}];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];
Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]]];,


ToString[FunctionName]==="H6",
If[Length[ParameterList]=!=3,Print["wrong parameterlist"];Abort[];];
{a0,b0,c0}=ParameterList;
{x0,y0}=Variables;
Rtof2pre =(1 + 4 x0)^(-a0/2) ((1 + 4 x0 + Sqrt[1 + 4 x0])/(2 + 8 x0))^b0  ;
Rtof2 = SeriesExpand["H2",Simplify/@{b0, c0, a0 + b0, 1 - a0 - b0, 1 - a0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{-(((1 + 4 x0 + Sqrt[1 + 4 x0]) y0)/(2 Sqrt[1 + 4 x0])), (-1 - 4 x0 + Sqrt[
 1 + 4 x0])/(2 + 8 x0)}}];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];
Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]]];,


ToString[FunctionName]==="H7",
If[Length[ParameterList]=!=4,Print["wrong parameterlist"];Abort[];];
{a0,b0,c0,d0}=ParameterList;
If[Denominator[d0/.eps-> 0]=!=2,Print["The parameter d in H7[a,b,c,d,x,y] must have the form p*eps + q/2\np,q are nonzero-integers "];Abort[]];
{x0,y0}=Variables;
Rtof2pre =(1 + 2 Sqrt[x0])^-a0 ;
Rtof2 = SeriesExpand["H2",Simplify/@{a0, -(1/2) + d0, b0, c0, -1 + 2 d0},{xx,yy},eps,m0];
temp=Series[Rtof2pre Rtof2,{eps,0,p0}]/.MapThread[Rule,{{xx,yy},{(4 Sqrt[x0])/(1 + 2 Sqrt[x0]), y0 + 2 Sqrt[x0] y0}}];
temp1 = ReplacePart[temp,3-> ParallelMap[PolyLogTools`LiToG,temp[[3]]]];
Return[ReplacePart[temp1,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,temp1[[3]]]]];


];





(* step1 - step2 : check if the expansion is of Taylor or Laurent type. If it is of Laurent type, use HYPERDIRE to find the appropriate differential
operator that relates the primary function with the secondary function*)


t0=SessionTime[];
step1 = ToHYPERDIREInput[FunctionName,Simplify/@ParameterList,{uu,vv},eps];
t1=SessionTime[];
PrintTemporary["Time : step 1 : ", t1-t0];

(*Print[step1/.MapThread[Rule,{{uu,vv},Variables}]];*)

step2 = FromHYPERDIREOutput[Simplify/@step1];
t2=SessionTime[];
PrintTemporary["Time : step 2 : ", t2-t1];
Print["HYPERDIRE o/p received ", SessionTime[]-t0];

(*Print[step2/.MapThread[Rule,{{uu,vv},Variables}]];*)

ind=ToSeriesData[#,eps,1]&/@(Series[#,{eps,0,0}]&/@(step2[[1]]));
(*Print[ind[[All,4]]];*)

(*N0 is same as lowest order*)
N0=Abs[Min[ind[[All,4]]]];
m0=N0+p0-1;
(*Print[{N0,p0,m0}];*)
(*Print[N0];
Print[IntegerQ[N0]];*)
	
Print["Differential operator has pole in e of order : ", N0];


operators = Which[Or[ToString[FunctionName]==="F2",ToString[FunctionName]==="F3",ToString[FunctionName]==="H2"],
{#&,uu PolyLogTools`DG[#,uu]&,vv PolyLogTools`DG[#,vv]&, uu vv PolyLogTools`DG[PolyLogTools`DG[#,uu],vv]&},
ToString[FunctionName]==="F1",
{#&,uu PolyLogTools`DG[#,uu]&,vv PolyLogTools`DG[#,vv]&}
];

(*Print[operators/.MapThread[Rule,{{uu,vv},Variables}]];*)

(* step 3 : find the canonical form of the Pfaffian form of the PDE of the function using CANONICA*)
step3= ToCANONICA[FunctionName,step2[[2,1;;-3]],{uu,vv},eps];
If[step3===False,Print["CANONICA fails to bring the system in canonical form"];Return[];];

t3=SessionTime[];
PrintTemporary["Time : step 3 : ", t3-t2];
Print["CANONICA o/p received ", SessionTime[]-t0];


(* Step4 - step5 : solve the PDE using PolyLogTools to find the expansion of the secondary function*)
step4 = SolutionVector[FunctionName,step3,{uu,vv},eps,m0+1];
t4=SessionTime[];
PrintTemporary["Time : step 4 : ", t4-t3];


(*step5 = Solution[step3,step4,eps,m0+1,m0+1];*)
step5 = Solution[step3,step4,eps,m0,m0];
PrintTemporary["Taylor expansion is done ", SessionTime[]-t0];
t5=SessionTime[];
PrintTemporary["Time : step 5 : ", t5-t4];
(*Print[step5/.MapThread[Rule,{{uu,vv},Variables}]];*)

If[DeleteDuplicates[step1[[2]]]==={0},Return[step5=ReplacePart[step5,-2-> step5[[-3]]+p0];step5/.MapThread[Rule,{{uu,vv},Variables}]]];

SetSharedVariable[operators,step2,step5];
(* step 6 : apply the differential operator to the expansion of the secondary function*)
step6= ParallelTable[If[step2[[1,i]]===0,0,step2[[1,i]](operators[[i]]@(Table[eps^(i-1),{i,1,Length[step5[[3]]]}]step5[[3]]))],{i,1,Length[operators]}];(*check if this is dependent on FunctionName*)
t6=SessionTime[];
PrintTemporary["Time : step 6 : ", t6-t5];
SetSharedVariable[step6];

(*Print[step6/.MapThread[Rule,{{uu,vv},Variables}]];*)

(* step 7 - step8 : Simplify the expression and expand in series*)
(*Print[step6/.MapThread[Rule,{{uu,vv},Variables}]];*)
step7 = ParallelTable[PolyLogTools`GCoefficientSimplify[Plus@@step6[[i]]],{i,1,Length[operators]}];
step7=Flatten[List@@@step7];
t7=SessionTime[];
PrintTemporary["Time : step 7 : ", t7-t6];
SetSharedVariable[step7];

(*Print[step7/.MapThread[Rule,{{uu,vv},Variables}]];*)

NN0=p0-N0;
SetSharedVariable[NN0];
step8=ParallelMap[Series[#,{eps,0,p0}]&,step7];
(*Print[step8/.MapThread[Rule,{{uu,vv},Variables}]];*)

t8=SessionTime[];
PrintTemporary["Time : step 8 : ", t8-t7];

(*Print[step8/.MapThread[Rule,{{uu,vv},Variables}]];*)
	
SetSharedVariable[step8];	
step9=Parallelize[Plus@@Flatten[step8]];
(*Print[step9/.MapThread[Rule,{{uu,vv},Variables}]];*)
step9=ReplacePart[step9,-2->step9[[-3]]+ p0];
t9=SessionTime[];
PrintTemporary["Time : step 9 : ", t9-t8];

(*Print[step9/.MapThread[Rule,{{uu,vv},Variables}]];*)


SetSharedVariable[step9];
step10=ReplacePart[step9,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,step9[[3]]]];
t10=SessionTime[];
PrintTemporary["Time : step 10 : ", t10-t9];


Return[step10/.MapThread[Rule,{{uu,vv},Variables}]];
]


(* ::Section:: *)
(*Commands for three variables (not using)*)


(*(*For three variables*)
SolutionVector[FunctionName_, exp_List(*ToCANONICA output*), variables_List, e_Symbol, p0_Integer]/;Length[variables]===3:=Module[{de,x,y,z,Omega,result1,i1,i2,i3,i4,i5,i6,
i7,i8,i9,i10,i11,i12,i13,i14,i15,i16,i17,i18,i19,BC,lowestorder},
{x,y,z}=variables;
SetSharedVariable[x,y,z];
Omega = exp[[2]];

BC = ToSeriesData[#,e,p0+1]&/@(Series[#,{e,0,p0}]&/@Which[
ToString[FunctionName]==="Fd3",Normal[Series[#,{x,0,0}]]&/@(Normal[Series[#,{y,0,0}]]&/@(Normal[Series[#,{z,0,0}]]&/@(Simplify[Inverse[exp[[1]]] . {1,0,0,0}])))]);			
(*Print[Normal[Series[#,{x,0,0}]]&/@(Normal[Series[#,{y,0,0}]]&/@(Normal[Series[#,{z,0,0}]]&/@(Simplify[Inverse[exp[[1]]] . {1,0,0,0}])))];	*)
(*Print[BC];*)
lowestorder=Min[BC[[All,4]]];
	
(*Print[lowestorder]*);
(*setting the boundary condition*)
i16=SeriesCoefficient[#,lowestorder]&/@BC;
(*i9 =If[ToString[FunctionName]==="F2",SeriesCoefficient[#,lowestorder]&/@BC];*)
(*i9 =If[ToString[FunctionName]==="F2",If[#==={},0,#[[1]]]&/@((ToSeriesData[#,e,p0]&/@BC)[[All,3]])(*SeriesCoefficient[#,0]&/@BC*)];*)
(*result1 =If[ToString[FunctionName]==="F2",{i9}];*)
result1={i16};
(*Print[i16];*)



(* i1 to i9 are steps for solving the DE*)
Table[
de =(# . i16)&/@Omega/.e-> 1;
SetSharedVariable[de];
i1 = ParallelMap[PolyLogTools`GCoefficientSimplify,de[[1]]];
SetSharedVariable[i1];
i2=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{x,y,z}]&,If[Head[i1[[i]]]===Plus,List@@(i1[[i]]),{i1[[i]]}]]],{i,1,Length[i1]}];
SetSharedVariable[i2];
(*Print[i2];*)
i3 = ParallelMap[PolyLogTools`GIntegrate[#,x]&,i2];
SetSharedVariable[i3];
i4=  ParallelMap[PolyLogTools`GCoefficientSimplify,i3];
SetSharedVariable[i4];
(*Print[i3];*)
i5 =de[[2]]- ParallelMap[PolyLogTools`DG[#,y]&,i4];
SetSharedVariable[i5];
i6 = ParallelMap[PolyLogTools`GCoefficientFullSimplify,i5];
SetSharedVariable[i6];
i7=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{y,z,x}]&,If[Head[i6[[i]]]===Plus,List@@(i6[[i]]),{i6[[i]]}]]],{i,1,Length[i6]}];
SetSharedVariable[i7];
i8 = ParallelMap[PolyLogTools`GIntegrate[#,y]&,i7];
SetSharedVariable[i8];
i9= ParallelMap[PolyLogTools`GCoefficientFullSimplify,i8];
SetSharedVariable[i9];
i10 =  ParallelMap[PolyLogTools`DG[#,z]&,i4];
SetSharedVariable[i10];
i11 =  ParallelMap[PolyLogTools`DG[#,z]&,i9];
SetSharedVariable[i11];
i12 =de[[3]]-i10-i11;
SetSharedVariable[i12];
i13 = ParallelMap[PolyLogTools`GCoefficientFullSimplify,i12];
SetSharedVariable[i13];
i14=Table[Total[ParallelMap[PolyLogTools`ToFibrationBasis[#,{z,x,y}]&,If[Head[i13[[i]]]===Plus,List@@(i13[[i]]),{i13[[i]]}]]],{i,1,Length[i13]}];
SetSharedVariable[i14];
i15 = ParallelMap[PolyLogTools`GIntegrate[#,z]&,i14];
SetSharedVariable[i15];
i16=  ParallelMap[PolyLogTools`GCoefficientSimplify,i15+i4+i9];
i16=i16+(SeriesCoefficient[#,lowestorder+i]&/@BC);
result1 = AppendTo[result1,i16];,
{i,1,p0}];
Return[{lowestorder,result1/.SeriesCoefficient[a_,n_Integer]/;FreeQ[a,e]&&n>=1-> 0}];
];*)


(*(*Not using it*)

(*For 3 variables*)
SeriesExpand[FunctionName_,ParameterList_List,Variables_List,eps_Symbol,p0_Integer]/;Length[Variables]===3&&FunctionName=!=FD3:=Module[
{m0=If[p0>6,6,p0],uu,vv,ww,operators,step1,step2,step3,step4,step5,step6,step7,step8,step9,step10,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,ind,N0,NN0,
tempser,mm,nn,pp},

If[MemberQ[ParameterList],eps===False,Return[]];

(* step1 - step2 : check if the expansion is of Taylor or Laurent type. If it is of Laurent type, use HYPERDIRE to find the appropriate differential
operator that relates the primary function with the secondary function*)



t0=SessionTime[];
step1 = ToHYPERDIREInput[FunctionName,ParameterList,{uu,vv,ww},eps];
t1=SessionTime[];
Print["Time : step 1 : ", t1-t0];
	

		
step2 = FromHYPERDIREOutput[step1];
t2=SessionTime[];
Print["Time : step 2 : ", t2-t1];
Print["HYPERDIRE o/p received ", SessionTime[]-t0];

(*Print[step2/.MapThread[Rule,{{uu,vv,ww},Variables}]];*)


ind=ToSeriesData[#,eps,1]&/@(Series[#,{eps,0,0}]&/@(step2[[1]]));
(*Print[ind[[All,4]]];*)
N0=Abs[Min[ind[[All,4]]]];
(*Print[N0];
Print[IntegerQ[N0]];*)
	
Print["Differential operator has pole in e of order : ", N0];

(*Print[ind/.MapThread[Rule,{{uu,vv,ww},Variables}]];*)



operators = Which[
ToString[FunctionName]==="FD3",
{#&,uu PolyLogTools`DG[#,uu]&,vv PolyLogTools`DG[#,vv]&,ww PolyLogTools`DG[#,ww]&}
];
		
(*Print[operators/.MapThread[Rule,{{uu,vv},Variables}]];*)
Print["CANONICA running"];

(* step 3 : find the canonical form of the Pfaffian form of the PDE of the function using CANONICA*)
step3= ToCANONICA[FunctionName,Flatten[step2[[2]]][[1;;-4]],{uu,vv,ww},eps];
If[step3===False,Print["CANONICA fails to bring the system in canonical form"];Return[];];
(*Print[step3/.MapThread[Rule,{{uu,vv,ww},Variables}]];*)
t3=SessionTime[];
Print["Time : step 3 : ", t3-t2];
Print["CANONICA o/p received ", SessionTime[]-t0];


(* Step4 - step5 : solve the PDE using PolyLogTools to find the expansion of the secondary function*)
step4 = SolutionVector[FunctionName,step3,{uu,vv,ww},eps,m0];
t4=SessionTime[];
Print["Time : step 4 : ", t4-t3];

	
step5 = Solution[step3,step4,eps,m0,m0];
Print["Taylor expansion is done ", SessionTime[]-t0];
t5=SessionTime[];
Print["Time : step 5 : ", t5-t4];

(*Print[step5/.MapThread[Rule,{{uu,vv,ww},Variables}]];*)

If[DeleteDuplicates[Flatten[step1[[2]]]]==={0},Return[step5/.MapThread[Rule,{{uu,vv,ww},Variables}]]];


SetSharedVariable[operators,step2,step5];
(* step 6 : apply the differential operator to the expansion of the secondary function*)
step6= ParallelTable[If[step2[[1,i]]===0,0,step2[[1,i]](operators[[i]]@(List@@Normal[step5]))],{i,1,Length[operators]}];(*check if this is dependent on FunctionName*)
t6=SessionTime[];
Print["Time : step 6 : ", t6-t5];
SetSharedVariable[step6];
(*Print[step6/.MapThread[Rule,{{uu,vv,ww},Variables}]];*)

(* step 7 - step8 : Simplify the expression and expand in series*)
(*Print[step6/.MapThread[Rule,{{uu,vv},Variables}]];*)
step7 = ParallelTable[PolyLogTools`GCoefficientSimplify[Plus@@step6[[i]]],{i,1,Length[operators]}];
(*step7=Flatten[List@@@step7];*)
t7=SessionTime[];
Print["Time : step 7 : ", t7-t6];
(*Print[step7/.MapThread[Rule,{{uu,vv,ww},Variables}]];*)

SetSharedVariable[step7];

NN0=p0-N0;
SetSharedVariable[NN0];
step8=ParallelMap[Series[#,{eps,0,NN0}]&,step7];
t8=SessionTime[];
Print["Time : step 8 : ", t8-t7];

(*Print[step8/.MapThread[Rule,{{uu,vv,ww},Variables}]];*)
	
SetSharedVariable[step8];	
step9=Parallelize[Plus@@Flatten[step8]];
t9=SessionTime[];
Print["Time : step 9 : ", t9-t8];




SetSharedVariable[step9];
step10=ReplacePart[step9,3-> ParallelMap[PolyLogTools`GCoefficientFullSimplify,step9[[3]]]];
t10=SessionTime[];
Print["Time : step 10 : ", t10-t9];


Return[step10/.MapThread[Rule,{{uu,vv,ww},Variables}]];
]*)


(* ::Section:: *)
(*Old Codes*)


(*SeriesExpand[FunctionName_,ParameterList_List,Variables_List,eps_Symbol,p0_Integer]/;Or[Length[Variables]===2,ToString[FunctionName]==="G1"]:=Module[
{g1f2,g1f2pre,a0,b0,c0,x0,y0,xx,yy},

Which[ToString[FunctionName]==="G1",
{a0,b0,c0}=ParameterList;
{x0,y0}=Variables;
g1f2pre = (1/(1+x0+y0))^a0 ;
Print[MapThread[Rule,{{Global`w,Global`z},{(1+2 x0-Sqrt[1-4 x0 y0])/(2 (1+x0+y0)),(1+2 y0-Sqrt[1-4 x0 y0])/(2 (1+x0+y0))}}]];
Print[{{1-b0-c0,a0,a0,1-b0,1-c0},{xx,yy}}];
g1f2 = SeriesExpand["F2",{1-b0-c0,a0,a0,1-b0,1-c0},{xx,yy},eps,5];
Return[Series[g1f2pre g1f2,{eps,0,p0}]];

];

];*)


(*OLD CODE v1.3*)

(*parameterlist[exp_,eps_Symbol]:={Coefficient[exp,eps,0],Coefficient[exp,eps,1]};
shiftparameter[exp_List,eps_Symbol]:=Module[{},Return[Total[{1,eps}If[exp[[2]]===0,exp,{1,exp[[2]]}]]]];
ToHYPERDIREInput[FunctionName_,ParameterList_List,{z1_,z2_},eps_Symbol]:=Module[{shiftlist, noshiftlist,newpara1,newpara2,newpara3,HYPintegerlist},
{noshiftlist,shiftlist}=Which[ToString[FunctionName]==="F1",{ParameterList[[1;;-3]],ParameterList[[-2;;-1]]},
Or[ToString[FunctionName]==="F2",ToString[FunctionName]==="F3"],{ParameterList[[1;;-2]],ParameterList[[-1;;-1]]}];
newpara1 = parameterlist[#,eps]&/@shiftlist;
newpara2 = shiftparameter[#,eps]&/@newpara1;
newpara3=Join[noshiftlist,newpara2] ;
Return[{ToExpression["AppellF1F4`"<>ToString[FunctionName]<>"IndexChange"],newpara3-ParameterList,
Join[ParameterList,{z1,z2}]}]
]*)


(*OLD CODE v1.3*)

(*FromHYPERDIREOutput[exp_]:=Module[{},
Return[Together[exp[[1]]@@(exp[[2;;3]])]]
]*)


(*SeriesExpand command for v1.0 *)


(*SeriesExpand[FunctionName_,ParameterList_List,variables_List,e_Symbol,p0_Integer]:=Module[{m0=6,step1,step2,step3},
step1= ToCANONICA[FunctionName,ParameterList,variables,e];
(*Print[step1];*)
step2 =  SolutionVector[FunctionName,step1,variables,e,m0];
(*Print[step2];*)
step3 = Solution[step1,step2,e,p0,m0];
Return[step3];
];*)


(* ::Section::Closed:: *)
(*Pfaff Systems*)


(* ::Subsection::Closed:: *)
(*F1	*)


mat1F1[{a_,b1_,b2_,c_},{x_,y_}]={{0,1/x,0},{(a b1)/(1-x),(-1+c+b2 (-1+x)-(a+b1) x)/((-1+x) x)+b2/(-x+y),(b1 (-1+y))/((-1+x) (x-y))},{0,(b2 y)/(x^2-x y),b1/(-x+y)}};


mat2F1[{a_,b1_,b2_,c_},{x_,y_}]={{0,0,1/y},{0,b2/(x-y),(b1 x)/(y (-x+y))},{(a b2)/(1-y),(b2 (-1+x))/((-1+y) (-x+y)),(1+a+b2-c)/(1-y)+b1/(x-y)+(1+b1-c)/y}};


(* ::Subsection::Closed:: *)
(*F2*)


mat1F2[{a_,b1_,b2_,c1_,c2_},{x_,y_}]  ={{0,1/x,0,0},{-((a b1)/(-1+x)),(-1-a-b1+c1)/(-1+x)+(1-c1)/x,-(b1/(-1+x)),-(1/(-1+x))},{0,0,0,1/x},{(a b1 b2)/(-1+x)-(a b1 b2)/(-1+x+y),(b2 (1+a+b1-c1))/(-1+x)-(b2 (1+a+b1-c1))/(-1+x+y),(b1 b2)/(-1+x)-(b1 (1+a+b2-c2))/(-1+x+y),-((1-c1-x-b2 x+c1 x)/((-1+x) x))+(-2-a-b1-b2+c1+c2)/(-1+x+y)}};


mat2F2[{a_,b1_,b2_,c1_,c2_},{x_,y_}] = {{0,0,1/y,0},{0,0,0,1/y},{-((a b2)/(-1+y)),-(b2/(-1+y)),(-1-a-b2+c2)/(-1+y)+(1-c2)/y,-(1/(-1+y))},{(a b1 b2)/(-1+y)-(a b1 b2)/(-1+x+y),(b1 b2)/(-1+y)+(-b2-a b2-b1 b2+b2 c1)/(-1+x+y),(b1 (1+a+b2-c2))/(-1+y)-(b1 (1+a+b2-c2))/(-1+x+y),b1/(-1+y)+(1-c2)/y+(-2-a-b1-b2+c1+c2)/(-1+x+y)}};


(* ::Subsection::Closed:: *)
(*F3*)


mat1F3[{a1_,a2_,b1_,b2_,c_},{x_,y_}]={{0,1/x,0,0},{-((a1 b1)/(-1+x)),(-1-a1-b1+c)/(-1+x)+(1-c)/x,0,1/(-1+x)-1/x},{0,0,0,1/x},{0,-((a2 b2)/((-1+x) x))-(a2 b2)/((-1+x) (-x-y+x y)),-((a1 b1)/(-1+x))-(a1 b1)/((-1+x) (-x-y+x y)),-((1+a2+b2-c+a1 x+b1 x)/((-1+x) x))+(-1-a1-a2-b1-b2+c)/((-1+x) (-x-y+x y))}};


mat2F3[{a1_,a2_,b1_,b2_,c_},{x_,y_}]={{0,0,1/y,0},{0,0,0,1/y},{-((a2 b2)/(-1+y)),0,(-1-a2-b2+c)/(-1+y)+(1-c)/y,1/(-1+y)-1/y},{0,-((a2 b2 (-1+x))/(-x-y+x y)),(a1 b1)/y-(a1 b1 (-1+x))/(-x+(-1+x) y),(1+a1+b1-c)/y-((1+a1+b1-c) (-1+x))/(-x+(-1+x) y)-((a2+b2) (-1+x))/(-x-y+x y)}};


(* ::Subsection::Closed:: *)
(*H2*)


mat1H2[{a_,b_,g1_,g2_,e_},{x_,y_}]={{0,1/x,0,0},{-((a b)/(-1+x)),(-1-a-b+e)/(-1+x)+(1-e)/x,b/(-1+x),1/(-1+x)},{0,0,0,1/x},{-((b g1 g2)/(-1+x))-(b g1 g2)/((-1+x) (-1-y+x y)),-((g1 g2)/(-1+x))-(g1 g2)/((-1+x) (-1-y+x y)),-((b (a+g1+g2))/(-1+x))-(b (a+g1+g2))/((-1+x) (-1-y+x y)),-((1-e+a x+b x+g1 x+g2 x)/((-1+x) x))+(-1-a-b+e-g1-g2)/((-1+x) (-1-y+x y))}};


mat2H2[{a_,b_,g1_,g2_,e_},{x_,y_}]={{0,0,1/y,0},{0,0,0,1/y},{-((g1 g2)/(1+y)),0,a/y+(-a-g1-g2)/(1+y),1/y-1/(1+y)},{(b g1 g2)/(1+y)-(b g1 g2 (-1+x))/(-1-y+x y),-((g1 g2 (-1+x))/(-1-y+x y)),(b (a+g1+g2))/(1+y)-(b (a+g1+g2) (-1+x))/(-1-y+x y),(1+a-e)/y+b/(1+y)-((1+a+b-e+g1+g2) (-1+x))/(-1-y+x y)}};


(* ::Subsection::Closed:: *)
(*FA3*)


FA3ser[{a_,b1_,b2_,b3_,c1_,c2_,c3_},{x_,y_,z_},{m_,n_,p_}]=(Pochhammer[a,m+n+p]Pochhammer[b1,m]Pochhammer[b2,n]Pochhammer[b3,p]x^m y^n z^p)/(Pochhammer[c1,m]Pochhammer[c2,n]Pochhammer[c3,p]m!n!p!);


mat1FA3[{a_,b1_,b2_,b3_,c1_,c2_,c3_},{x_,y_,z_}]={{0,1/x,0,0,0,0,0,0},{-((a b1)/(-1+x)),(-1-a-b1+c1)/(-1+x)+(1-c1)/x,-(b1/(-1+x)),-(b1/(-1+x)),-(1/(-1+x)),0,-(1/(-1+x)),0},{0,0,0,0,1/x,0,0,0},{0,0,0,0,0,0,1/x,0},{(a b1 b2)/(-1+x)-(a b1 b2)/(-1+x+y),(b2 (1+a+b1-c1))/(-1+x)-(b2 (1+a+b1-c1))/(-1+x+y),(b1 b2)/(-1+x)-(b1 (1+a+b2-c2))/(-1+x+y),(b1 b2)/(-1+x)-(b1 b2)/(-1+x+y),-((1-c1-x-b2 x+c1 x)/((-1+x) x))+(-2-a-b1-b2+c1+c2)/(-1+x+y),-(b1/(-1+x+y)),b2/(-1+x)-b2/(-1+x+y),-(1/(-1+x+y))},{0,0,0,0,0,0,0,1/x},{(a b1 b3)/(-1+x)-(a b1 b3)/(-1+x+z),(b3 (1+a+b1-c1))/(-1+x)-(b3 (1+a+b1-c1))/(-1+x+z),(b1 b3)/(-1+x)-(b1 b3)/(-1+x+z),(b1 b3)/(-1+x)-(b1 (1+a+b3-c3))/(-1+x+z),b3/(-1+x)-b3/(-1+x+z),-(b1/(-1+x+z)),-((1-c1-x-b3 x+c1 x)/((-1+x) x))+(-2-a-b1-b3+c1+c3)/(-1+x+z),-(1/(-1+x+z))},{-((a b1 b2 b3 y)/((-1+x) (-1+x+y)))+(a b1 b2 b3)/(-1+x+z)-(a b1 b2 b3)/(-1+x+y+z),-((b2 b3 (1+a+b1-c1) y)/((-1+x) (-1+x+y)))+(b2 b3 (1+a+b1-c1))/(-1+x+z)-(b2 b3 (1+a+b1-c1))/(-1+x+y+z),(b1 b3 (-1-a+c2+x+a x-c2 x-b2 y))/((-1+x) (-1+x+y))+(b1 b2 b3)/(-1+x+z)-(b1 b3 (1+a+b2-c2))/(-1+x+y+z),-((b1 b2 b3 y)/((-1+x) (-1+x+y)))+(b1 b2 (1+a+b3-c3))/(-1+x+z)-(b1 b2 (1+a+b3-c3))/(-1+x+y+z),(b3 (-2-a-b1+c1+c2+2 x+a x+b1 x-c1 x-c2 x-b2 y))/((-1+x) (-1+x+y))+(b2 b3)/(-1+x+z)-(b3 (2+a+b1+b2-c1-c2))/(-1+x+y+z),(b1 b3)/(-1+x+y)+(b1 b2)/(-1+x+z)-(b1 (2+a+b2+b3-c2-c3))/(-1+x+y+z),-((b2 b3 y)/((-1+x) (-1+x+y)))+(b2 (2+a+b1+b3-c1-c3))/(-1+x+z)-(b2 (2+a+b1+b3-c1-c3))/(-1+x+y+z),-((1-c1-x-b3 x+c1 x-y+c1 y)/(x (-1+x+y)))+b2/(-1+x+z)+(-3-a-b1-b2-b3+c1+c2+c3)/(-1+x+y+z)}};


mat2FA3[{a_,b1_,b2_,b3_,c1_,c2_,c3_},{x_,y_,z_}]={{0,0,1/y,0,0,0,0,0},{0,0,0,0,1/y,0,0,0},{-((a b2)/(-1+y)),-(b2/(-1+y)),(-1-a-b2+c2)/(-1+y)+(1-c2)/y,-(b2/(-1+y)),-(1/(-1+y)),-(1/(-1+y)),0,0},{0,0,0,0,0,1/y,0,0},{(a b1 b2)/(-1+y)-(a b1 b2)/(-1+x+y),(b1 b2)/(-1+y)-(b2 (1+a+b1-c1))/(-1+x+y),(b1 (1+a+b2-c2))/(-1+y)-(b1 (1+a+b2-c2))/(-1+x+y),(b1 b2)/(-1+y)-(b1 b2)/(-1+x+y),b1/(-1+y)+(1-c2)/y+(-2-a-b1-b2+c1+c2)/(-1+x+y),b1/(-1+y)-b1/(-1+x+y),-(b2/(-1+x+y)),-(1/(-1+x+y))},{(a b2 b3)/(-1+y)-(a b2 b3)/(-1+y+z),(b2 b3)/(-1+y)-(b2 b3)/(-1+y+z),(b3 (1+a+b2-c2))/(-1+y)-(b3 (1+a+b2-c2))/(-1+y+z),(b2 b3)/(-1+y)-(b2 (1+a+b3-c3))/(-1+y+z),b3/(-1+y)-b3/(-1+y+z),-((1-c2-y-b3 y+c2 y)/((-1+y) y))+(-2-a-b2-b3+c2+c3)/(-1+y+z),-(b2/(-1+y+z)),-(1/(-1+y+z))},{0,0,0,0,0,0,0,1/y},{-((a b1 b2 b3 x)/((-1+y) (-1+x+y)))+(a b1 b2 b3)/(-1+y+z)-(a b1 b2 b3)/(-1+x+y+z),(b2 b3 (-1-a+c1-b1 x+y+a y-c1 y))/((-1+y) (-1+x+y))+(b1 b2 b3)/(-1+y+z)-(b2 b3 (1+a+b1-c1))/(-1+x+y+z),-((b1 b3 (1+a+b2-c2) x)/((-1+y) (-1+x+y)))+(b1 b3 (1+a+b2-c2))/(-1+y+z)-(b1 b3 (1+a+b2-c2))/(-1+x+y+z),-((b1 b2 b3 x)/((-1+y) (-1+x+y)))+(b1 b2 (1+a+b3-c3))/(-1+y+z)-(b1 b2 (1+a+b3-c3))/(-1+x+y+z),(b3 (-2-a-b2+c1+c2-b1 x+2 y+a y+b2 y-c1 y-c2 y))/((-1+y) (-1+x+y))+(b1 b3)/(-1+y+z)-(b3 (2+a+b1+b2-c1-c2))/(-1+x+y+z),-((b1 b3 x)/((-1+y) (-1+x+y)))+(b1 (2+a+b2+b3-c2-c3))/(-1+y+z)-(b1 (2+a+b2+b3-c2-c3))/(-1+x+y+z),(b2 b3)/(-1+x+y)+(b1 b2)/(-1+y+z)-(b2 (2+a+b1+b3-c1-c3))/(-1+x+y+z),-((1-c2-x+c2 x-y-b3 y+c2 y)/(y (-1+x+y)))+b1/(-1+y+z)+(-3-a-b1-b2-b3+c1+c2+c3)/(-1+x+y+z)}};


mat3FA3[{a_,b1_,b2_,b3_,c1_,c2_,c3_},{x_,y_,z_}]={{0,0,0,1/z,0,0,0,0},{0,0,0,0,0,0,1/z,0},{0,0,0,0,0,1/z,0,0},{-((a b3)/(-1+z)),-(b3/(-1+z)),-(b3/(-1+z)),(-1-a-b3+c3)/(-1+z)+(1-c3)/z,0,-(1/(-1+z)),-(1/(-1+z)),0},{0,0,0,0,0,0,0,1/z},{(a b2 b3)/(-1+z)-(a b2 b3)/(-1+y+z),(b2 b3)/(-1+z)-(b2 b3)/(-1+y+z),(b2 b3)/(-1+z)-(b3 (1+a+b2-c2))/(-1+y+z),(b2 (1+a+b3-c3))/(-1+z)-(b2 (1+a+b3-c3))/(-1+y+z),-(b3/(-1+y+z)),b2/(-1+z)+(1-c3)/z+(-2-a-b2-b3+c2+c3)/(-1+y+z),b2/(-1+z)-b2/(-1+y+z),-(1/(-1+y+z))},{(a b1 b3)/(-1+z)-(a b1 b3)/(-1+x+z),(b1 b3)/(-1+z)-(b3 (1+a+b1-c1))/(-1+x+z),(b1 b3)/(-1+z)-(b1 b3)/(-1+x+z),(b1 (1+a+b3-c3))/(-1+z)-(b1 (1+a+b3-c3))/(-1+x+z),-(b3/(-1+x+z)),b1/(-1+z)-b1/(-1+x+z),b1/(-1+z)+(1-c3)/z+(-2-a-b1-b3+c1+c3)/(-1+x+z),-(1/(-1+x+z))},{-((a b1 b2 b3)/(-1+z))+(a b1 b2 b3)/(-1+x+z)+(a b1 b2 b3)/(-1+y+z)-(a b1 b2 b3)/(-1+x+y+z),-((b1 b2 b3)/(-1+z))+(b2 b3 (1+a+b1-c1))/(-1+x+z)+(b1 b2 b3)/(-1+y+z)-(b2 b3 (1+a+b1-c1))/(-1+x+y+z),-((b1 b2 b3)/(-1+z))+(b1 b2 b3)/(-1+x+z)+(b1 b3 (1+a+b2-c2))/(-1+y+z)-(b1 b3 (1+a+b2-c2))/(-1+x+y+z),-((b1 b2 (1+a+b3-c3))/(-1+z))+(b1 b2 (1+a+b3-c3))/(-1+x+z)+(b1 b2 (1+a+b3-c3))/(-1+y+z)-(b1 b2 (1+a+b3-c3))/(-1+x+y+z),(b2 b3)/(-1+x+z)+(b1 b3)/(-1+y+z)-(b3 (2+a+b1+b2-c1-c2))/(-1+x+y+z),-((b1 b2)/(-1+z))+(b1 b2)/(-1+x+z)+(b1 (2+a+b2+b3-c2-c3))/(-1+y+z)-(b1 (2+a+b2+b3-c2-c3))/(-1+x+y+z),-((b1 b2)/(-1+z))+(b2 (2+a+b1+b3-c1-c3))/(-1+x+z)+(b1 b2)/(-1+y+z)-(b2 (2+a+b1+b3-c1-c3))/(-1+x+y+z),(1-c3)/z+b2/(-1+x+z)+b1/(-1+y+z)+(-3-a-b1-b2-b3+c1+c2+c3)/(-1+x+y+z)}};


(* ::Subsection::Closed:: *)
(*FB3*)


FB3ser[{a1_,a2_,a3_,b1_,b2_,b3_,c_},{x_,y_,z_},{m_,n_,p_}]=(Pochhammer[a1,m]Pochhammer[a2,n]Pochhammer[a3,p]Pochhammer[b1,m]Pochhammer[b2,n]Pochhammer[b3,p]x^m y^n z^p)/(Pochhammer[c,m+n+p]m!n!p!);


mat1FB3[{a1_,a2_,a3_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,1/x,0,0,0,0,0,0},{-((a1 b1)/(-1+x)),(-1-a1-b1+c)/(-1+x)+(1-c)/x,0,0,1/(-1+x)-1/x,0,1/(-1+x)-1/x,0},{0,0,0,0,1/x,0,0,0},{0,0,0,0,0,0,1/x,0},{0,-((a2 b2)/((-1+x) x))-(a2 b2)/((-1+x) (-x-y+x y)),-((a1 b1)/(-1+x))-(a1 b1)/((-1+x) (-x-y+x y)),0,-((1+a2+b2-c+a1 x+b1 x)/((-1+x) x))+(-1-a1-a2-b1-b2+c)/((-1+x) (-x-y+x y)),0,0,1/((-1+x) x)+1/((-1+x) (-x-y+x y))},{0,0,0,0,0,0,0,1/x},{0,-((a3 b3)/((-1+x) x))-(a3 b3)/((-1+x) (-x-z+x z)),0,-((a1 b1)/(-1+x))-(a1 b1)/((-1+x) (-x-z+x z)),0,0,-((1+a3+b3-c+a1 x+b1 x)/((-1+x) x))+(-1-a1-a3-b1-b3+c)/((-1+x) (-x-z+x z)),1/((-1+x) x)+1/((-1+x) (-x-z+x z))},{0,0,0,0,-((a3 b3 y)/(x (-x-y+x y)))-(a3 b3 y^2)/((-x-y+x y) (-x y-x z-y z+x y z)),-((a1 b1 (-1+y))/(-x-y+x y))-(a1 b1 y^2)/((-x-y+x y) (-x y-x z-y z+x y z)),-((a2 b2 y)/(x (-x-y+x y)))-(a2 b2 y^2)/((-x-y+x y) (-x y-x z-y z+x y z)),-((-a1 x-b1 x+y+a2 y+a3 y+b2 y+b3 y-c y+a1 x y+b1 x y)/(x (-x-y+x y)))-((1+a1+a2+a3+b1+b2+b3-c) y^2)/((-x-y+x y) (-x y-x z-y z+x y z))}};


mat2FB3[{a1_,a2_,a3_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,0,1/y,0,0,0,0,0},{0,0,0,0,1/y,0,0,0},{-((a2 b2)/(-1+y)),0,(-1-a2-b2+c)/(-1+y)+(1-c)/y,0,1/(-1+y)-1/y,1/(-1+y)-1/y,0,0},{0,0,0,0,0,1/y,0,0},{0,-((a2 b2 (-1+x))/(-x-y+x y)),(a1 b1)/y-(a1 b1 (-1+x))/(-x-y+x y),0,(1+a1+b1-c)/y+(1+a1+a2+b1+b2-c-x-a1 x-a2 x-b1 x-b2 x+c x)/(-x-y+x y),0,0,-(1/y)+(-1+x)/(-x-y+x y)},{0,0,-((a3 b3)/((-1+y) y))-(a3 b3)/((-1+y) (-y-z+y z)),-((a2 b2)/(-1+y))-(a2 b2)/((-1+y) (-y-z+y z)),0,-((1+a3+b3-c+a2 y+b2 y)/((-1+y) y))+(-1-a2-a3-b2-b3+c)/((-1+y) (-y-z+y z)),0,1/((-1+y) y)+1/((-1+y) (-y-z+y z))},{0,0,0,0,0,0,0,1/y},{0,0,0,0,-((a3 b3 x)/(y (-x-y+x y)))-(a3 b3 x^2)/((-x-y+x y) (-x y-x z-y z+x y z)),-((a1 b1 x)/(y (-x-y+x y)))-(a1 b1 x^2)/((-x-y+x y) (-x y-x z-y z+x y z)),-((a2 b2 (-1+x))/(-x-y+x y))-(a2 b2 x^2)/((-x-y+x y) (-x y-x z-y z+x y z)),-((x+a1 x+a3 x+b1 x+b3 x-c x-a2 y-b2 y+a2 x y+b2 x y)/(y (-x-y+x y)))-((1+a1+a2+a3+b1+b2+b3-c) x^2)/((-x-y+x y) (-x y-x z-y z+x y z))}};


mat3FB3[{a1_,a2_,a3_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,0,0,1/z,0,0,0,0},{0,0,0,0,0,0,1/z,0},{0,0,0,0,0,1/z,0,0},{-((a3 b3)/(-1+z)),0,0,(-1-a3-b3+c)/(-1+z)+(1-c)/z,0,1/(-1+z)-1/z,1/(-1+z)-1/z,0},{0,0,0,0,0,0,0,1/z},{0,0,-((a3 b3 (-1+y))/(-y-z+y z)),(a2 b2)/z-(a2 b2 (-1+y))/(-y-z+y z),0,(1+a2+b2-c)/z+(1+a2+a3+b2+b3-c-y-a2 y-a3 y-b2 y-b3 y+c y)/(-y-z+y z),0,-(1/z)+(-1+y)/(-y-z+y z)},{0,-((a3 b3 (-1+x))/(-x-z+x z)),0,(a1 b1)/z-(a1 b1 (-1+x))/(-x-z+x z),0,0,(1+a1+b1-c)/z+(1+a1+a3+b1+b3-c-x-a1 x-a3 x-b1 x-b3 x+c x)/(-x-z+x z),-(1/z)+(-1+x)/(-x-z+x z)},{0,0,0,0,-((a3 b3 (-x-y+x y))/(-x y-x z-y z+x y z)),(a1 b1)/z-(a1 b1 (-x-y+x y))/(-x y-x z-y z+x y z),(a2 b2)/z-(a2 b2 (-x-y+x y))/(-x y-x z-y z+x y z),(1+a1+a2+b1+b2-c)/z+1/(-x y-x z-y z+x y z) (x+a1 x+a2 x+a3 x+b1 x+b2 x+b3 x-c x+y+a1 y+a2 y+a3 y+b1 y+b2 y+b3 y-c y-x y-a1 x y-a2 x y-a3 x y-b1 x y-b2 x y-b3 x y+c x y)}};


(* ::Subsection::Closed:: *)
(*FD3*)


FD3ser[{a_,b1_,b2_,b3_,c_},{x_,y_,z_},{m_,n_,p_}]=(Pochhammer[a,m+n+p]Pochhammer[b1,m]Pochhammer[b2,n]Pochhammer[b3,p]x^m y^n z^p)/(Pochhammer[c,m+n+p]m!n!p!);


mat1FD3[{a_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,1/x,0,0},{-((a b1)/(-1+x)),(-x-b3 x+c x-a x^2-b1 x^2+b3 x^2+y+b2 y+b3 y-c y+a x y+b1 x y-b2 x y-b3 x y)/((-1+x) x (x-y))+b3/(-x+z),-(b1/(-1+x))-b1/(-x+y),-(b1/(-1+x))-b1/(-x+z)},{0,-(b2/x)-b2/(-x+y),b1/(-x+y),0},{0,-(b3/x)-b3/(-x+z),0,b1/(-x+z)}};


mat2FD3[{a_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,0,1/y,0},{0,-(b2/(-x+y)),-(b1/y)+b1/(-x+y),0},{-((a b2)/(-1+y)),-(b2/(-1+y))+b2/(-x+y),(x+b1 x+b3 x-c x-y-b3 y+c y+a x y-b1 x y+b2 x y-b3 x y-a y^2-b2 y^2+b3 y^2)/((-1+y) y (-x+y))+b3/(-y+z),-(b2/(-1+y))-b2/(-y+z)},{0,0,-(b3/y)-b3/(-y+z),b2/(-y+z)}};


mat3FD3[{a_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,0,0,1/z},{0,-(b3/(-x+z)),0,-(b1/z)+b1/(-x+z)},{0,0,-(b3/(-y+z)),-(b2/z)+b2/(-y+z)},{-((a b3)/(-1+z)),-(b3/(-1+z))+b3/(-x+z),-(b3/(-1+z))+b3/(-y+z),(-1-a-b3+c)/(-1+z)+(1+b1+b2-c)/z-b1/(-x+z)-b2/(-y+z)}};


(* ::Subsection::Closed:: *)
(*FN3*)


FN3ser[{a1_,a2_,a3_,b1_,b2_,c1_,c2_},{x_,y_,z_},{m_,n_,p_}]=(Pochhammer[a1,m]Pochhammer[a2,n]Pochhammer[a3,p]Pochhammer[b1,m+p]Pochhammer[b2,n]x^m y^n z^p)/(Pochhammer[c1,m]Pochhammer[c2,n+p]m!n!p!);


mat1FN3[{a1_,a2_,a3_,b1_,b2_,c1_,c2_},{x_,y_,z_}]={{0,1/x,0,0,0,0,0,0},{-((a1 b1)/(-1+x)),(1-c1+a1 x+b1 x)/(x-x^2),0,a1/(1-x),0,0,1/(1-x),0},{0,0,0,0,1/x,0,0,0},{0,0,0,0,0,0,1/x,0},{0,0,-((a1 b1)/(-1+x)),0,(1-c1+a1 x+b1 x)/(x-x^2),a1/(1-x),0,1/(1-x)},{0,0,0,0,0,0,0,1/x},{(a1 a3 b1 z)/((-1+x) (-1+x+z)),(a3 (1+a1+b1-c1) z)/((-1+x) (-1+x+z)),0,(a1 (1+b1+c2 (-1+x)-x-b1 x+a3 z))/((-1+x) (-1+x+z)),0,a1/(-1+x+z),(1+b1 x-c2 x-a1 (-1+x) x-x^2-b1 x^2+c2 x^2-z+x z+a3 x z+c1 (-1+x+z-x z))/((-1+x) x (-1+x+z)),1/(-1+x+z)},{0,0,(a1 a3 b1 (-1+y) z)/((-1+x) (-z+y (-1+x+z))),(a1 a2 b2 y)/(z-y (-1+x+z)),(a3 (1+a1+b1-c1) (-1+y) z)/((-1+x) (-z+y (-1+x+z))),(a1 (-a3 z+y (1+a2+b1+b2-c2-x-a2 x-b1 x-b2 x+c2 x+a3 z)))/((-1+x) (-z+y (-1+x+z))),(a2 b2 y)/(z-y (-1+x+z)),((1+c1 (-1+x)-(1+a3) x) z-y (-1-a2 x-b1 x-b2 x+c2 x+a1 (-1+x) x+x^2+a2 x^2+b1 x^2+b2 x^2-c2 x^2+c1 (-1+x) (-1+z)+z-x z-a3 x z))/((-1+x) x (-z+y (-1+x+z)))}};


mat2FN3[{a1_,a2_,a3_,b1_,b2_,c1_,c2_},{x_,y_,z_}]={{0,0,1/y,0,0,0,0,0},{0,0,0,0,1/y,0,0,0},{(a2 b2)/(1-y),0,(1-c2+a2 y+b2 y)/(y-y^2),0,0,1/((-1+y) y),0,0},{0,0,0,0,0,1/y,0,0},{0,(a2 b2)/(1-y),0,0,(1-c2+a2 y+b2 y)/(y-y^2),0,0,1/((-1+y) y)},{0,0,(a3 b1 z)/(y (y+z-y z)),-((a2 b2 (-1+z))/(y (-1+z)-z)),(a3 z)/(y (y+z-y z)),(-(1+a3+b1-c2) z+a2 (y-y z)+b2 (y-y z))/(y (y (-1+z)-z)),0,z/(y (y+z-y z))},{0,0,0,0,0,0,0,1/y},{0,0,(a1 a3 b1 x z)/((y (-1+z)-z) (-z+y (-1+x+z))),(a1 a2 b2 x z)/((y (-1+z)-z) (-z+y (-1+x+z))),-((a3 z ((-1-b1+c1) z+y (-1+c1-a1 x+b1 (-1+z)+z-c1 z)))/(y (y (-1+z)-z) (-z+y (-1+x+z)))),(a1 (1+a2+a3+b1+b2-c2) x z)/((y (-1+z)-z) (-z+y (-1+x+z))),-((a2 b2 (-1+x+z))/(-z+y (-1+x+z))),(-a2 y (y (-1+z)-z) (-1+x+z)-b2 y (y (-1+z)-z) (-1+x+z)+z ((2+a3+b1-c1-c2) z+y (2+a3+b1-c1-c2+a1 x-2 z-a3 z-b1 z+c1 z+c2 z)))/(y (y (-1+z)-z) (-z+y (-1+x+z)))}};


mat3FN3[{a1_,a2_,a3_,b1_,b2_,c1_,c2_},{x_,y_,z_}]={{0,0,0,1/z,0,0,0,0},{0,0,0,0,0,0,1/z,0},{0,0,0,0,0,1/z,0,0},{(a3 b1)/(1-z),a3/(1-z),0,(1-c2+a3 z+b1 z)/(z-z^2),0,1/((-1+z) z),1/(1-z),0},{0,0,0,0,0,0,0,1/z},{0,0,-((a3 b1 (-1+y))/(y (-1+z)-z)),(a2 b2 y)/(z (y+z-y z)),-((a3-a3 y)/(y+z-y z)),((a3+b1) z-y (1+a2+b2-c2+a3 z+b1 z))/((y (-1+z)-z) z),0,-((1-y)/(y+z-y z))},{(a1 a3 b1 x)/((-1+z) (-1+x+z)),(a3 (1+b1+a1 x+c1 (-1+z)-z-b1 z))/((-1+z) (-1+x+z)),0,(a1 (1+a3+b1-c2) x)/((-1+z) (-1+x+z)),0,-((a1 x)/((-1+z) (-1+x+z))),(x (-1+z+a1 z)-(-1+z) (1+(1+a3+b1-c1) z)+c2 (-1+x+z-x z))/((-1+z) z (-1+x+z)),(1-x)/(z (-1+x+z))},{0,0,(a1 a3 b1 x (-1+y) y)/((y (-1+z)-z) (-z+y (-1+x+z))),(a1 a2 b2 x (-1+y) y)/((y (-1+z)-z) (-z+y (-1+x+z))),-((a3 (-1+y) ((-1-b1+c1) z+y (-1+c1-a1 x+b1 (-1+z)+z-c1 z)))/((y (-1+z)-z) (-z+y (-1+x+z)))),(a1 (1+a2+a3+b1+b2-c2) x (-1+y) y)/((y (-1+z)-z) (-z+y (-1+x+z))),(a2 b2 (-1+x) y)/(z (-z+y (-1+x+z))),(-(1+a3+b1-c1) z^2-y z (a3+b1-b2-c1+c2+a2 (-1+x)+x+a1 x+b2 x-c2 x-2 z-2 a3 z-2 b1 z+2 c1 z)+y^2 (1-c2-x+c2 x+a2 (-1+x) (-1+z)+b2 (-1+x) (-1+z)+a3 z+b1 z-c1 z+c2 z+x z+a1 x z-c2 x z-z^2-a3 z^2-b1 z^2+c1 z^2))/((y (-1+z)-z) z (-z+y (-1+x+z)))}};


(* ::Subsection::Closed:: *)
(*FS3*)


FS3ser[{a1_,a2_,b1_,b2_,b3_,c_},{x_,y_,z_},{m_,n_,p_}]=(Pochhammer[a1,m]Pochhammer[a2,n+p]Pochhammer[b1,m]Pochhammer[b2,n]Pochhammer[b3,p]x^m y^n z^p)/(Pochhammer[c,m+n+p]m!n!p!);


mat1FS3[{a1_,a2_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,1/x,0,0,0,0},{-((a1 b1)/(-1+x)),(-1-a1-b1+c)/(-1+x)+(1-c)/x,0,0,1/(-1+x)-1/x,1/(-1+x)-1/x},{0,0,0,0,1/x,0},{0,0,0,0,0,1/x},{0,-((a2 b2)/((-1+x) x))-(a2 b2)/((-1+x) (-x-y+x y)),-((a1 b1)/(-1+x))-(a1 b1)/((-1+x) (-x-y+x y)),0,-((1+a2+b2-c+a1 x+b1 x)/((-1+x) x))+(-1-a1-a2-b1-b2+c)/((-1+x) (-x-y+x y)),-(b2/((-1+x) x))-b2/((-1+x) (-x-y+x y))},{0,-((a2 b3)/((-1+x) x))-(a2 b3)/((-1+x) (-x-z+x z)),0,-((a1 b1)/(-1+x))-(a1 b1)/((-1+x) (-x-z+x z)),-(b3/((-1+x) x))-b3/((-1+x) (-x-z+x z)),-((1+a2+b3-c+a1 x+b1 x)/((-1+x) x))+(-1-a1-a2-b1-b3+c)/((-1+x) (-x-z+x z))}};


mat2FS3[{a1_,a2_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,0,1/y,0,0,0},{0,0,0,0,1/y,0},{-((a2 b2)/(-1+y)),0,(-1-b3+c-a2 y-b2 y+b3 y)/((-1+y) y)+b3/(-y+z),-(b2/(-1+y))-b2/(-y+z),1/(-1+y)-1/y,0},{0,0,-(b3/y)-b3/(-y+z),b2/(-y+z),0,0},{0,-((a2 b2 (-1+x))/(-x-y+x y)),(a1 b1)/y-(a1 b1 (-1+x))/(-x-y+x y),0,(-x-a1 x-b1 x-b3 x+c x+a2 y+b2 y-b3 y-a2 x y-b2 x y+b3 x y)/(y (-x-y+x y))+b3/(-y+z),-((b2 (-1+x))/(-x-y+x y))-b2/(-y+z)},{0,0,0,0,-(b3/y)-b3/(-y+z),b2/(-y+z)}};


mat3FS3[{a1_,a2_,b1_,b2_,b3_,c_},{x_,y_,z_}]={{0,0,0,1/z,0,0},{0,0,0,0,0,1/z},{0,0,-(b3/(-y+z)),-(b2/z)+b2/(-y+z),0,0},{-((a2 b3)/(-1+z)),0,-(b3/(-1+z))+b3/(-y+z),(-1-a2-b3+c)/(-1+z)+(1+b2-c)/z-b2/(-y+z),0,1/(-1+z)-1/z},{0,0,0,0,-(b3/(-y+z)),-(b2/z)+b2/(-y+z)},{0,-((a2 b3 (-1+x))/(-x-z+x z)),0,(a1 b1)/z-(a1 b1 (-1+x))/(-x-z+x z),b3/(-y+z)-(b3 (-1+x))/(-x-z+x z),(1+a1+b1+b2-c)/z-b2/(-y+z)-((1+a1+a2+b1+b3-c) (-1+x))/(-x-z+x z)}};


(* ::Subsection::Closed:: *)
(*FK3*)


FK3ser[{a1_,a2_,b1_,b2_,c1_,c2_,c3_},{x_,y_,z_},{m_,n_,p_}]=(Pochhammer[a1,m]Pochhammer[a2,n+p]Pochhammer[b1,m+p]Pochhammer[b2,n]x^m y^n z^p)/(Pochhammer[c1,m]Pochhammer[c2,n]Pochhammer[c3,p]m!n!p!);


mat1FK3[{a1_,a2_,b1_,b2_,c1_,c2_,c3_},{x_,y_,z_}]={{0,1/x,0,0,0,0,0,0},{-((a1 b1)/(-1+x)),(1-c1+a1 x+b1 x)/(x-x^2),0,a1/(1-x),0,0,1/(1-x),0},{0,0,0,0,1/x,0,0,0},{0,0,0,0,0,0,1/x,0},{0,0,-((a1 b1)/(-1+x)),0,(1-c1+a1 x+b1 x)/(x-x^2),a1/(1-x),0,1/(1-x)},{0,0,0,0,0,0,0,1/x},{(a1 a2 b1 z)/((-1+x) (-1+x+z)),(a2 (1+a1+b1-c1) z)/((-1+x) (-1+x+z)),(a1 b1 z)/((-1+x) (-1+x+z)),(a1 (1+b1+c3 (-1+x)-x-b1 x+a2 z))/((-1+x) (-1+x+z)),((1+a1+b1-c1) z)/((-1+x) (-1+x+z)),(a1 z)/((-1+x) (-1+x+z)),(1+b1 x-c3 x-a1 (-1+x) x-x^2-b1 x^2+c3 x^2-z+x z+a2 x z+c1 (-1+x+z-x z))/((-1+x) x (-1+x+z)),z/((-1+x) (-1+x+z))},{-((a1 a2 b1 b2 y z)/((1+x (-1+y)-y-z) (-1+x+z))),-((a2 b2 (1+a1+b1-c1) y z)/((1+x (-1+y)-y-z) (-1+x+z))),-((a1 b1 z (-1+x-b2 y+b2 x y+z+a2 (-1+x+z)-c2 (-1+x+z)))/((-1+x) (1+x (-1+y)-y-z) (-1+x+z))),-((a1 b2 (1+a2+b1-c3) y z)/((1+x (-1+y)-y-z) (-1+x+z))),-(((1+a1+b1-c1) z (-1+x-b2 y+b2 x y+z+a2 (-1+x+z)-c2 (-1+x+z)))/((-1+x) (1+x (-1+y)-y-z) (-1+x+z))),-((a1 (-1+2 x-x^2+y-2 x y+x^2 y-a2 z+c2 z+a2 x z-c2 x z-y z-b2 y z+x y z+b2 x y z+z^2+a2 z^2-c2 z^2+b1 (-1+x) (-1+y) (-1+x+z)-c3 (-1+x) (-1+y) (-1+x+z)))/((-1+x) (1+x (-1+y)-y-z) (-1+x+z))),-((b2 (2+a1+a2+b1-c1-c3) y z)/((1+x (-1+y)-y-z) (-1+x+z))),((-1-a1-b1+c3) x^3 (-1+y)+(-1+z) (-1+y+z)+c1 (-1+x) (-1+x+z) (-1+y+z)-x^2 (1-2 c3-y+2 c3 y+a1 (-1+y) (-2+z)+b1 (-1+y) (-2+z)+z+a2 z-c2 z+c3 z+y z+b2 y z-c3 y z)+x (-1-c3+y+c3 y+a1 (-1+y) (-1+z)+b1 (-1+y) (-1+z)+3 z+a2 z-c2 z+c3 z+b2 y z-c3 y z-2 z^2-a2 z^2+c2 z^2))/((-1+x) x (1+x (-1+y)-y-z) (-1+x+z))}};


mat2FK3[{a1_,a2_,b1_,b2_,c1_,c2_,c3_},{x_,y_,z_}]={{0,0,1/y,0,0,0,0,0},{0,0,0,0,1/y,0,0,0},{-((a2 b2)/(-1+y)),0,(1-c2+a2 y+b2 y)/(y-y^2),b2/(1-y),0,1/(1-y),0,0},{0,0,0,0,0,1/y,0,0},{0,-((a2 b2)/(-1+y)),0,0,(1-c2+a2 y+b2 y)/(y-y^2),0,b2/(1-y),1/(1-y)},{(a2 b1 b2 z)/((-1+y) (-1+y+z)),(a2 b2 z)/((-1+y) (-1+y+z)),(b1 (1+a2+b2-c2) z)/((-1+y) (-1+y+z)),(b2 (1+a2+c3 (-1+y)-y-a2 y+b1 z))/((-1+y) (-1+y+z)),((1+a2+b2-c2) z)/((-1+y) (-1+y+z)),(1+b2 y-c3 y-a2 (-1+y) y-y^2-b2 y^2+c3 y^2-z+y z+b1 y z+c2 (-1+y+z-y z))/((-1+y) y (-1+y+z)),(b2 z)/((-1+y) (-1+y+z)),z/((-1+y) (-1+y+z))},{0,0,0,0,0,0,0,1/y},{-((a1 a2 b1 b2 x z)/((1+x (-1+y)-y-z) (-1+y+z))),-((a2 b2 z (-1-a1 x+y+a1 x y+z+b1 (-1+y+z)-c1 (-1+y+z)))/((-1+y) (1+x (-1+y)-y-z) (-1+y+z))),-((a1 b1 (1+a2+b2-c2) x z)/((1+x (-1+y)-y-z) (-1+y+z))),-((a1 b2 (1+a2+b1-c3) x z)/((1+x (-1+y)-y-z) (-1+y+z))),-(((1+a2+b2-c2) z (-1-a1 x+y+a1 x y+z+b1 (-1+y+z)-c1 (-1+y+z)))/((-1+y) (1+x (-1+y)-y-z) (-1+y+z))),-((a1 (2+a2+b1+b2-c2-c3) x z)/((1+x (-1+y)-y-z) (-1+y+z))),-((b2 (-1+x+2 y-2 x y-y^2+x y^2-b1 z+c1 z-x z-a1 x z+b1 y z-c1 y z+x y z+a1 x y z+z^2+b1 z^2-c1 z^2+a2 (-1+x) (-1+y) (-1+y+z)-c3 (-1+x) (-1+y) (-1+y+z)))/((-1+y) (1+x (-1+y)-y-z) (-1+y+z))),(c2 (-1+y) (-1+x+z) (-1+y+z)+(-1+y+z) (-1+c3 y+a2 (-1+y) y+b2 (-1+y) y+y^2-c3 y^2+z-2 y z-b1 y z+c1 y z)-x (-1+y) (-1+c3 y+y^2-c3 y^2+z+y z+a1 y z-c3 y z+a2 y (-1+y+z)+b2 y (-1+y+z)))/((-1+y) y (1+x (-1+y)-y-z) (-1+y+z))}};


mat3FK3[{a1_,a2_,b1_,b2_,c1_,c2_,c3_},{x_,y_,z_}]={{0,0,0,1/z,0,0,0,0},{0,0,0,0,0,0,1/z,0},{0,0,0,0,0,1/z,0,0},{(a2 b1)/(1-z),a2/(1-z),b1/(1-z),(1-c3+a2 z+b1 z)/(z-z^2),1/(1-z),1/(1-z),1/(1-z),0},{0,0,0,0,0,0,0,1/z},{(a2 b1 b2 y)/((-1+z) (-1+y+z)),(a2 b2 y)/((-1+z) (-1+y+z)),(b1 (1+a2+b2 y+c2 (-1+z)-z-a2 z))/((-1+z) (-1+y+z)),(b2 (1+a2+b1-c3) y)/((-1+z) (-1+y+z)),(1+a2+b2 y+c2 (-1+z)-z-a2 z)/((-1+z) (-1+y+z)),(y (-1+z+b2 z)-(-1+z) (1+(1+a2+b1-c2) z)+c3 (-1+y+z-y z))/((-1+z) z (-1+y+z)),(b2 y)/((-1+z) (-1+y+z)),-(1/(-1+y+z))},{(a1 a2 b1 x)/((-1+z) (-1+x+z)),(a2 (1+b1+a1 x+c1 (-1+z)-z-b1 z))/((-1+z) (-1+x+z)),(a1 b1 x)/((-1+z) (-1+x+z)),(a1 (1+a2+b1-c3) x)/((-1+z) (-1+x+z)),(1+b1+a1 x+c1 (-1+z)-z-b1 z)/((-1+z) (-1+x+z)),(a1 x)/((-1+z) (-1+x+z)),(x (-1+z+a1 z)-(-1+z) (1+(1+a2+b1-c1) z)+c3 (-1+x+z-x z))/((-1+z) z (-1+x+z)),-(1/(-1+x+z))},{-((a1 a2 b1 b2 x y (1+x (-1+y)-y-z^2))/((1+x (-1+y)-y-z) (-1+z) (-1+x+z) (-1+y+z))),(a2 b2 y (-1+x-a1 x+a1 x^2+y-x y+a1 x y-a1 x^2 y+2 z-2 x z-y z+x y z-z^2+x z^2+a1 x z^2+b1 (-1+x) (-1+z) (-1+y+z)-c1 (-1+x) (-1+z) (-1+y+z)))/((1+x (-1+y)-y-z) (-1+z) (-1+x+z) (-1+y+z)),(a1 b1 x (-1+x+y-b2 y-x y+b2 x y+b2 y^2-b2 x y^2+2 z-x z-2 y z+x y z-z^2+y z^2+b2 y z^2+a2 (-1+y) (-1+z) (-1+x+z)-c2 (-1+y) (-1+z) (-1+x+z)))/((1+x (-1+y)-y-z) (-1+z) (-1+x+z) (-1+y+z)),-((a1 b2 (1+a2+b1-c3) x y (1+x (-1+y)-y-z^2))/((1+x (-1+y)-y-z) (-1+z) (-1+x+z) (-1+y+z))),(-1+c1+c2-c1 c2+x-a1 x-c1 x-c2 x+a1 c2 x+c1 c2 x+a1 x^2-a1 c2 x^2+y-b2 y-c1 y+b2 c1 y-c2 y+c1 c2 y-x y+a1 x y+b2 x y-a1 b2 x y+c1 x y-b2 c1 x y+c2 x y-a1 c2 x y-c1 c2 x y-a1 x^2 y+a1 b2 x^2 y+a1 c2 x^2 y+b2 y^2-b2 c1 y^2-b2 x y^2+a1 b2 x y^2+b2 c1 x y^2-a1 b2 x^2 y^2+3 z-3 c1 z-3 c2 z+3 c1 c2 z-2 x z+2 a1 x z+2 c1 x z+2 c2 x z-2 a1 c2 x z-2 c1 c2 x z-a1 x^2 z+a1 c2 x^2 z-2 y z+2 b2 y z+2 c1 y z-2 b2 c1 y z+2 c2 y z-2 c1 c2 y z+x y z-2 a1 x y z-2 b2 x y z-c1 x y z+2 b2 c1 x y z-c2 x y z+2 a1 c2 x y z+c1 c2 x y z+a1 x^2 y z-a1 c2 x^2 y z-b2 y^2 z+b2 c1 y^2 z+b2 x y^2 z-b2 c1 x y^2 z-3 z^2+3 c1 z^2+3 c2 z^2-3 c1 c2 z^2+x z^2-a1 x z^2-c1 x z^2-c2 x z^2+a1 c2 x z^2+c1 c2 x z^2+y z^2-b2 y z^2-c1 y z^2+b2 c1 y z^2-c2 y z^2+c1 c2 y z^2+a1 x y z^2+b2 x y z^2+a1 b2 x y z^2-b2 c1 x y z^2-a1 c2 x y z^2+z^3-c1 z^3-c2 z^3+c1 c2 z^3-b1 (-1+z) (-1+y+z) (1+b2 y-x (1+b2 y)-z+c2 (-1+x+z))+a2 (-1+z) (-1+x+z) (-1-a1 x+y+a1 x y+z+b1 (-1+y+z)-c1 (-1+y+z)))/((1+x (-1+y)-y-z) (-1+z) (-1+x+z) (-1+y+z)),(a1 x (-2+c2+c3+2 x-c2 x-c3 x+2 y-b2 y-c2 y-c3 y-2 x y+b2 x y+c2 x y+c3 x y+b2 y^2-b2 x y^2+4 z-2 c2 z-2 c3 z-2 x z+c2 x z+c3 x z-4 y z+2 c2 y z+2 c3 y z+2 x y z-c2 x y z-c3 x y z-2 z^2+c2 z^2+c3 z^2+2 y z^2+b2 y z^2-c2 y z^2-c3 y z^2+a2 (-1+y) (-1+z) (-1+x+z)+b1 (-1+y) (-1+z) (-1+x+z)))/((1+x (-1+y)-y-z) (-1+z) (-1+x+z) (-1+y+z)),(b2 y (-2+c1+c3+2 x-a1 x-c1 x-c3 x+a1 x^2+2 y-c1 y-c3 y-2 x y+a1 x y+c1 x y+c3 x y-a1 x^2 y+4 z-2 c1 z-2 c3 z-4 x z+2 c1 x z+2 c3 x z-2 y z+c1 y z+c3 y z+2 x y z-c1 x y z-c3 x y z-2 z^2+c1 z^2+c3 z^2+2 x z^2+a1 x z^2-c1 x z^2-c3 x z^2+a2 (-1+x) (-1+z) (-1+y+z)+b1 (-1+x) (-1+z) (-1+y+z)))/((1+x (-1+y)-y-z) (-1+z) (-1+x+z) (-1+y+z)),(-c3 (-1+x) (-1+y) (-1+x+z) (-1+y+z)+x^2 (-1+y) (-1+y+z+a1 z)-(-1+y+z) (y (-1+z+b2 z)-(-1+z) (1+(2+a2+b1-c1-c2) z))+x (y^2 (-2+z+b2 z)-(-1+z) (-2+(-1+a1-a2-b1+c1+c2) z)+y (4-(2+a1-a2-b1+b2+c1+c2) z+(1+a1+b2) z^2)))/((1+x (-1+y)-y-z) z (-1+x+z) (-1+y+z))}};


(* ::Subsection:: *)
(*FM3 *)


FM3ser[{a1_,a2_,b1_,b2_,c1_,c2_},{x_,y_,z_},{m_,n_,p_}]=(Pochhammer[a1,m]Pochhammer[a2,n+p]Pochhammer[b1,m+p]Pochhammer[b2,n]x^m y^n z^p)/(Pochhammer[c1,m]Pochhammer[c2,n+p]m!n!p!);


mat1FM3[{a1_,a2_,b1_,b2_,c1_,c2_},{x_,y_,z_}]={{0,1/x,0,0,0,0},{(a1 b1)/(1-x),(1-c1+a1 x+b1 x)/(x-x^2),0,a1/(1-x),0,1/(1-x)},{0,0,0,0,1/x,0},{0,0,0,0,0,1/x},{0,0,-((a1 b1 y)/((-1+x) y+z)),(a1 b2 y)/((-1+x) y+z),((-1+c1-a1 x-b1 x) y+z-c1 z)/(x ((-1+x) y+z)),(b2 y)/((-1+x) y+z)},{(a1 a2 b1 z)/((-1+x) (-1+x+z)),(a2 (1+a1+b1-c1) z)/((-1+x) (-1+x+z)),(a1 b1 (-1+y) z)/((-1+x+z) ((-1+x) y+z)),-((a1 (z (-1+c2+b1 (-1+x)+x-c2 x-a2 z)+(-1+x) y (-1+c2+b1 (-1+x)+x-c2 x-a2 z+b2 (-1+x+z))))/((-1+x) (-1+x+z) ((-1+x) y+z))),((1+a1+b1-c1) (-1+y) z)/((-1+x+z) ((-1+x) y+z)),-((z (-1-b1 x+c2 x+a1 (-1+x) x+x^2+b1 x^2-c2 x^2+c1 (-1+x) (-1+z)+z-x z-a2 x z)+(-1+x) y (-1-b1 x-b2 x+c2 x+a1 (-1+x) x+x^2+b1 x^2+b2 x^2-c2 x^2+c1 (-1+x) (-1+z)+z-x z-a2 x z+b2 x z))/((-1+x) x (-1+x+z) ((-1+x) y+z)))}};


mat2FM3[{a1_,a2_,b1_,b2_,c1_,c2_},{x_,y_,z_}]={{0,0,1/y,0,0,0},{0,0,0,0,1/y,0},{-((a2 b2)/(-1+y)),0,((-a2-b2) y^2+(1+b1-c2) z+y (-1+c2+a2 z-b1 z+b2 z))/((-1+y) y (y-z)),(b2 (-1+z))/((-1+y) (y-z)),-(z/(y^2-y z)),0},{0,0,(b1 z)/(y^2-y z),b2/(-y+z),z/(y^2-y z),0},{0,-((a2 b2)/(-1+y)),(a1 b1 x z)/((y-z) ((-1+x) y+z)),-((a1 b2 x z)/((y-z) ((-1+x) y+z))),((-a2-b2) (-1+x) y^3+(2+b1-c1-c2) z^2+y z (-3+c1+2 c2+x-a1 x-c2 x-z+a2 z+b2 z+c1 z-b1 (1+z))+y^2 (1+c2 (-1+x)+z-2 a2 z+b1 z-2 b2 z-c1 z+x (-1+a1 z+a2 z+b2 z)))/((-1+y) y (y-z) ((-1+x) y+z)),-((b2 (-1+x+z))/((-1+y) ((-1+x) y+z)))},{0,0,-((a1 b1 x z)/((y-z) ((-1+x) y+z))),(a1 b2 x z)/((y-z) ((-1+x) y+z)),(z ((-1-b1+c1-a1 x) y+(1+b1-c1) z))/(y (y-z) ((-1+x) y+z)),(b2-b2 x)/((-1+x) y+z)}};


mat3FM3[{a1_,a2_,b1_,b2_,c1_,c2_},{x_,y_,z_}]={{0,0,0,1/z,0,0},{0,0,0,0,0,1/z},{0,0,b1/(y-z),-((b2 y)/(y z-z^2)),1/(y-z),0},{(a2 b1)/(1-z),a2/(1-z),(b1-b1 y)/((y-z) (-1+z)),(-z (1-c2+a2 z+b1 z)+y (1+b2-c2+a2 z+b1 z-b2 z))/((-1+z) z (-y+z)),(1-y)/((y-z) (-1+z)),1/(1-z)},{0,0,-((a1 b1 x y)/((y-z) ((-1+x) y+z))),(a1 b2 x y)/((y-z) ((-1+x) y+z)),((-1-b1+c1-a1 x) y+(1+b1-c1) z)/((y-z) ((-1+x) y+z)),-((b2 (-1+x) y)/(z ((-1+x) y+z)))},{(a1 a2 b1 x)/((-1+z) (-1+x+z)),(a2 (1+b1+a1 x+c1 (-1+z)-z-b1 z))/((-1+z) (-1+x+z)),(a1 b1 x (-1+y) ((-1+x) y+z^2))/((y-z) (-1+z) (-1+x+z) ((-1+x) y+z)),(a1 x (-(((1+a2+b1-c2) (y-z))/((-1+z) (-1+x+z)))+(b2 y)/((-1+x) y+z)))/(-y+z),((-1+y) ((-1+x) y (-1+c1-a1 x+b1 (-1+z)+z-c1 z)+z (-1+x+c1 (-1+x) (-1+z)+z-x z-a1 x z+b1 (-1+x+z-x z))))/((-1+z) (-1+x+z) (-y+z) ((-1+x) y+z)),(-z (c2 (-1+x) (-1+z)-x (-1+z+a1 z)+(-1+z) (1+(1+a2+b1-c1) z))+(-1+x) y (1-x+a2 z+b1 z-c1 z+x z+a1 x z-z^2-a2 z^2-b1 z^2+c1 z^2+b2 (-1+z) (-1+x+z)+c2 (-1+x+z-x z)))/((-1+z) z (-1+x+z) ((-1+x) y+z))}};


(* ::Section:: *)
(*End Package*)


End[]


EndPackage[]
