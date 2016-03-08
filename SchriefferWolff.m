(* ::Package:: *)

BeginPackage["SchriefferWolff`"]

computePerturbativeTranformation::usage =
	"computePerturbativeTranformation[ inputH0, inputV, inputOrder] takes in a base Hamiltonian inputH0, a perturbation inputV, and calculates the SW Transformation to the inputOrder-th order."

Begin["Private`"]
computePerturbativeTranformation[inputH0_,inputV_,inputOrder_]:=Module[{H0=inputH0,V=inputV,order=inputOrder,HTot,dim,Hextra,Comm,Sadd,Hextra2,n,m,i,S2,H1,H2,Heff,Stot,A},(
HTot = H0+V;
dim=Dimensions[H0][[1]];
Comm[A_,B_]:=A.B-B.A;

Stot = Table [0,{n,1,dim},{m,1,dim}];
For[i=1,i<=order,i++,
Hextra = If[i==1,V,If[i==2,\[Minus]Comm[Stot,Comm[Stot,H0]]/2+I Comm[Stot,V],Throw[i]]];
Sadd=Simplify[Table [If[ n==m,0,-((I Hextra[[n,m]])/(H0[[n]][[n]]\[Minus]H0[[m]][[m]]))],{n,1,dim},{m,1,dim}]];
Stot = Stot + (A^i)*Sadd;
];
Normal[Series[MatrixExp[-I*Stot],{A,0,order}]]/.{A->1}
)]
End[]
EndPackage[]
