(* ::Package:: *)

BeginPackage["SchriefferWolff`"]
Off[Part::partd]
computePerturbativeTranformation::usage =
	"computePerturbativeTranformation[ inputH0, inputV, inputOrder] takes in a base Hamiltonian inputH0, a perturbation inputV, and calculates the SW Transformation to the inputOrder-th order."

Begin["Private`"]
computePerturbativeTranformation[inputH0_,inputV_,inputOrder_]:=Module[{H0Mat=inputH0,VMat=inputV,order=inputOrder,HTot,dim,Hextra,Comm,Sadd,Hextra2,n,m,i,S2,H1,H2,Heff,Stot,A,toReturn,H0,V},(

HTotMat = H0Mat+VMat;
dim=Dimensions[H0Mat][[1]];
Comm[A_,B_]:=A.B-B.A;

(******************Creates symbolic list of Hxs *********************)
Stot = 0;
For[i = 1, i <= order, i++, Stot = Stot + S[[i]] \[Lambda]^i;];
Heff = H0 + \[Lambda] V;

Adjoint[x_] := Comm[I*Stot, x];

For[i = 1, i <= order, i++, Heff = Heff + Nest[Adjoint, H0 + \[Lambda] V, i]/i!];

Hx[i_] := Simplify[SeriesCoefficient[TensorExpand[Heff, Assumptions -> \[Lambda] \[Element] Reals], {\[Lambda], 0, i}] + I*Comm[H0, S[[i]]]];

HextrasSYM = Array[Hx, order];
(***************************************)
Stot = Table [0,{n,1,dim},{m,1,dim}];
SMats = Table[Table [0,{n,1,dim},{m,1,dim}],order];
For[i=1,i<=order,i++,
Rules = {H0->H0Mat,V->VMat};
rule[i_] := S[[i]]->SMats[[i]];
Rules = Join[Rules,Array[rule,order]];
Hextra = HextrasSYM[[i]]/.Rules;
SMats[[i]]=Simplify[Table [If[ n==m,0,-((I Hextra[[n,m]])/(H0Mat[[n]][[n]]\[Minus]H0Mat[[m]][[m]]))],{n,1,dim},{m,1,dim}]];
Stot = Stot + (A^i)*SMats[[i]];
];
toReturn = IdentityMatrix[dim];
Do[toReturn = toReturn+MatrixPower[-I*Stot,i]/i!,{i,1,order,1}];
Normal[Series[toReturn,{A,0,order}]]/.{A->1}
)]
End[]

computeSubspaceTranformation::usage =
	"computeSubspaceTranformation[ inputHam, inputOrder, lastLowEnergyIndex] takes in a base Hamiltonian inputHam, and calculates the SW Transformation to the inputOrder-th order, separating the subspace that ends at index lastLowEnergyIndex."

Begin["Private`"]
computeSubspaceTranformation[inputHam_,inputOrder_,lastLowEnergyIndex_]:=Module[{HTotMat=inputHam,order=inputOrder,cutoff=lastLowEnergyIndex,VMat,H0Mat,dim,Hextra,Comm,Sadd,Hextra2,n,m,i,S2,H1,H2,Heff,Stot,A,toReturn,H0,V},(

dim=Dimensions[HTotMat][[1]];
VMat = Table[If[((m <= cutoff && n <= cutoff) || (m > cutoff && n > cutoff)),0, HTotMat[[n, m]]], {n, 1, dim}, {m, 1, dim}];
H0Mat = Table[If[(!(m <= cutoff && n <= cutoff) && !(m > cutoff && n > cutoff)),0, HTotMat[[n, m]]], {n, 1, dim}, {m, 1, dim}];

Comm[A_,B_]:=A.B-B.A;

(******************Creates symbolic list of Hxs *********************)
Stot = 0;
For[i = 1, i <= order, i++, Stot = Stot + S[[i]] \[Lambda]^i;];
Heff = H0 + \[Lambda] V;

Adjoint[x_] := Comm[I*Stot, x];

For[i = 1, i <= order, i++, Heff = Heff + Nest[Adjoint, H0 + \[Lambda] V, i]/i!];

Hx[i_] := Simplify[SeriesCoefficient[TensorExpand[Heff, Assumptions -> \[Lambda] \[Element] Reals], {\[Lambda], 0, i}] + I*Comm[H0, S[[i]]]];

HextrasSYM = Array[Hx, order];
(***************************************)
Stot = Table [0,{n,1,dim},{m,1,dim}];
SMats = Table[Table [0,{n,1,dim},{m,1,dim}],order];
For[i=1,i<=order,i++,
Rules = {H0->H0Mat,V->VMat};
rule[i_] := S[[i]]->SMats[[i]];
Rules = Join[Rules,Array[rule,order]];
Hextra = HextrasSYM[[i]]/.Rules;
SMats[[i]]=Simplify[Table [If[ ((m <= cutoff && n <= cutoff) || (m > cutoff && n > cutoff)),0,-((I Hextra[[n,m]])/(H0Mat[[n]][[n]]\[Minus]H0Mat[[m]][[m]]))],{n,1,dim},{m,1,dim}]];
Stot = Stot + (A^i)*SMats[[i]];
];
toReturn = IdentityMatrix[dim];
Do[toReturn = toReturn+MatrixPower[-I*Stot,i]/i!,{i,1,order,1}];
Normal[Series[toReturn,{A,0,order}]]/.{A->1}
)]
End[]
EndPackage[]
