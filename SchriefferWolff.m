(* ::Package:: *)

BeginPackage["SchriefferWolff`"]

computePerturbativeTranformation::usage =
	"computePerturbativeTranformation[ inputH0, inputV, inputOrder] takes in a base Hamiltonian inputH0, a perturbation inputV, and calculates the SW Transformation to the inputOrder-th order."

Begin["Private`"]
computePerturbativeTranformation[inputH0_,inputV_,inputOrder_]:=Module[{H0=inputH0,V=inputV,order=inputOrder,HTot,dim,Hextra,Comm,Sadd,Hextra2,n,m,i,S2,H1,H2,Heff,Stot,A,toReturn},(
HTot = H0+V;
dim=Dimensions[H0][[1]];
Comm[A_,B_]:=A.B-B.A;

Stot = Table [0,{n,1,dim},{m,1,dim}];
For[i=1,i<=order,i++,
Hextra = If[i==1,V,If[i==2,\[Minus]Comm[Stot,Comm[Stot,H0]]/2+I Comm[Stot,V],Throw[i]]];
Sadd=Simplify[Table [If[ n==m,0,-((I Hextra[[n,m]])/(H0[[n]][[n]]\[Minus]H0[[m]][[m]]))],{n,1,dim},{m,1,dim}]];
Stot = Stot + (A^i)*Sadd;
];
toReturn = IdentityMatrix[dim];
Do[toReturn = toReturn+MatrixPower[-I*Stot,i]/i!,{i,1,order,1}];
Normal[Series[toReturn,{A,0,order}]]/.{A->1}
)]
End[]

computeSubspaceTranformation::usage =
"computeSubspaceTranformation[ inputHam, inputOrder,lastLowEnergyIndex] takes in a base Hamiltonian inputH0, a perturbation inputV, and calculates the SW Transformation to the inputOrder-th order, separating the subspace that ends at index lastLowEnergyIndex."

Begin["Private`"]
computeSubspaceTranformation[inputHam_,inputOrder_,lastLowEnergyIndex_]:=Module[{HTot=inputHam,order=inputOrder,cutoff=lastLowEnergyIndex,H0,V,dim,Hextra,Comm,Sadd,Hextra2,n,m,i,S2,H1,H2,Heff,Stot,A,toReturn},(
dim=Dimensions[HTot][[1]];
Comm[A_,B_]:=A.B-B.A;

V = Table[If[((m <= cutoff && n <= cutoff) || (m > cutoff && n > cutoff)),0, HTot[[n, m]]], {n, 1, dim}, {m, 1, dim}];
H0 = Table[If[(!(m <= cutoff && n <= cutoff) && !(m > cutoff && n > cutoff)),0, HTot[[n, m]]], {n, 1, dim}, {m, 1, dim}];
Stot = Table [0,{n,1,dim},{m,1,dim}];
For[i=1,i<=order,i++,
Hextra = If[i==1,V,If[i==2,\[Minus]Comm[Stot,Comm[Stot,H0]]/2+I Comm[Stot,V],Throw[i]]];
Sadd=Simplify[Table [If[ ((m <= cutoff && n <= cutoff) || (m > cutoff && n > cutoff)),0,-((I Hextra[[n,m]])/(H0[[n]][[n]]\[Minus]H0[[m]][[m]]))],{n,1,dim},{m,1,dim}]];
Stot = Stot + (A^i)*Sadd;
];
toReturn = IdentityMatrix[dim];
Do[toReturn = toReturn+MatrixPower[-I*Stot,i]/i!,{i,1,order,1}];
Normal[Series[toReturn,{A,0,order}]]/.{A->1}
)]
End[]
EndPackage[]
