(* ::Package:: *)

BeginPackage["SchriefferWolff`"]

computePerturbativeTranformation::usage =
	"computePerturbativeTranformation[ inputH0, inputV, inputOrder] takes in a base Hamiltonian inputH0, a perturbation inputV, and calculates the SW Transformation to the inputOrder-th order."

Begin["Private`"]
computePerturbativeTranformation[inputH0_,inputV_,inputOrder_]:=Module[{H0=inputH0,V=inputV,order=inputOrder,HTot,dim,Hextra,Comm,Sadd,Hextra2,n,m,i,S2,H1,H2,Heff,Stot,A,toReturn},(
HTot = H0+V;
dim=Dimensions[H0][[1]];
Comm[A_,B_]:=A.B-B.A;

Stot = 0;
For[i = 1, i <= order, i++, Stot = Stot + S[[i]] \[Lambda]^i;];
Heff = H0 + \[Lambda] V;

Adjoint[x_] := Comm[I*Stot, x]

For[i = 1, i <= order, i++, Heff = Heff + Nest[Adjoint, H0 + \[Lambda] V, i]/i!];

Hx[i_] := Simplify[SeriesCoefficient[TensorExpand[Heff, Assumptions -> \[Lambda] \[Element] Reals], {\[Lambda], 0, i}] + I*Comm[H0, S[[i]]]];

Hextras = Array[Hx, order];

Stot = Table [0,{n,1,dim},{m,1,dim}];
For[i=1,i<=order,i++,
Hextra = Hextras[[i]];
Sadd=Simplify[Table [If[ n==m,0,-((I Hextra[[n,m]])/(H0[[n]][[n]]\[Minus]H0[[m]][[m]]))],{n,1,dim},{m,1,dim}]];
Stot = Stot + (A^i)*Sadd;
];
toReturn = IdentityMatrix[dim];
Do[toReturn = toReturn+MatrixPower[-I*Stot,i]/i!,{i,1,order,1}];
Normal[Series[toReturn,{A,0,order}]]/.{A->1}
)]
End[]
EndPackage[]
