(* ::Package:: *)

(* ::Input::Initialization:: *)
BeginPackage["InversePNUtils`"];


(* ::Section::Initialization:: *)
(*Utils functions definition*)


(* ::Input::Initialization:: *)
\[Eta][m1_,m2_]:=(m1 m2)/(m1+m2)^2


(* ::Input::Initialization:: *)
M[m1_,m2_]:=m1+m2


(* ::Input::Initialization:: *)
Mc[\[Eta]_,M_]:=M \[Eta]^0.6


(* ::Input::Initialization:: *)
(*M[\[Eta]_,Mc_]:=Mc \[Eta]^-0.6*)


(* ::Input::Initialization:: *)
\[Delta][x_List]:= Abs[x[[1]]-x[[2]]]/(x[[1]]+x[[2]]) /; 
Length[x]==2


(* ::Input::Initialization:: *)
\[Delta][\[Eta]_]:= ((1-4\[Eta])^0.5) 


(* ::Input:: *)
(*\[Delta][{27,98}]*)


(* ::Input:: *)
(*\[Delta][\[Eta][27,98]]*125/71*)


(* ::Input:: *)
(*\[Delta][{50,60}]*)


(* ::Input:: *)
(*\[Delta][\[Eta][50,60]]*11*)


(* ::Input::Initialization:: *)
(*\[Delta][\[Eta]_]:= (1-4\[Eta])^0.5 /; 
0<\[Eta]<=1/4*)


(* ::Input::Initialization:: *)
PiOneThirds[j_]:=(\[Pi])^(j/3)


(* ::Input::Initialization:: *)
PiMOneThirds[j_,M_]:=(\[Pi] M)^(j/3)


(* ::Input::Initialization:: *)
\[Chi]s[\[Chi]1L_,\[Chi]2L_]:=(\[Chi]1L+\[Chi]2L)/2


(* ::Input::Initialization:: *)
\[Chi]a[\[Chi]1L_,\[Chi]2L_]:=(\[Chi]1L-\[Chi]2L)/2


(* ::Input::Initialization:: *)
\[Chi]1L[\[Chi]s_,\[Chi]a_]:=\[Chi]s+\[Chi]a


(* ::Input::Initialization:: *)
\[Chi]2L[\[Chi]s_,\[Chi]a_]:=\[Chi]s-\[Chi]a


(* ::Input::Initialization:: *)
EG=N[EulerGamma,50]


(* ::Input::Initialization:: *)
m1s=Range[10,80,1/1000];
m2s=Range[10,80,1/1000];
etas=Range[1/1000,1/4,1/100];
chi1Ls=Range[-1,1-1/100,1/100];
chi2Ls=Range[-1,1-1/100,1/100];
chias=Range[-2,2-1/1000,1/1000];
chiss=Range[-2,2-1/1000,1/1000];

sparseChias=chias[[;; ;; 10]];
sparseChiss=chias[[;; ;; 10]];


(* ::Input::Initialization:: *)
normalizationFactor["eta"]=1/4;
normalizationFactor["chi1L"]=1;
normalizationFactor["chi2L"]=1;
normalizationFactor["chia"]=2;
normalizationFactor["chis"]=2;


(* ::Input::Initialization:: *)
G=6.6743015 10^-11 (* m^3kg^-1s^-2 *);
c=299792458.0(* m/s *);
solarMass=1.9884754153381438 10^30(* kg *);

(*timeMetric=(solarMass G)/c^3(* s *);*)

timeMetric=4.925490947641266978197229498498379006 10^-6;
TimeMetricOneThirds[j_]:=timeMetric^(j/3)(* s^(1/3) *);


(* ::Input:: *)
(*0.0026/timeMetric*)


(* ::Input:: *)
(*timeMetricOneThirds[3-5]*)


(* ::Section::Initialization:: *)
(*Common coefficient of psij*)


(* ::Input::Initialization:: *)
CE[\[Eta]_]:=(3/128) (1/\[Eta]) 


(* ::Section::Initialization:: *)
(*alpha-s*)
(*(IMRPhenomXAS, arXiv/2001.11412)*)
(*https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/_l_a_l_sim_i_m_r_phenom_x__internals_8c_source.html*)


(* ::Input::Initialization:: *)
Clear["alpha"]


(* ::Text::Initialization:: *)
(*Pay attention to the dimension of M*)


(* ::Input::Initialization:: *)
(*\:5b9a\:4e49\:540e\:725b\:987f\:76f8\:4f4d\:7cfb\:6570\:51fd\:6570*)
AlphaNS[j_,\[Eta]_,\[Chi]1L_,\[Chi]2L_]:=Switch[j,
(*0/2 PN*)0,1.0,
(*1/2 PN*)1,0.0,
(*2/2 PN*)2,3715/756+(55 \[Eta])/9,
(*3/2 PN*)3,-16\[Pi],
(*4/2 PN*)4,(15293365/508032+(27145 \[Eta])/504+(3085 \[Eta]^2)/72) ,
(*5/2 PN*)5,0.0,
(*6/2 PN*)6,(11583231236531/4694215680-(5 \[Eta] (3147553127+588 \[Eta] (-45633+102260 \[Eta])))/3048192-6848 EulerGamma/21-640 \[Pi]^2/3+(2255 \[Eta] \[Pi]^2)/12-13696 Log[2]/21-6848 Log[\[Pi]]/63) ,
(*7/2 PN*)7,5 (15419335+168 (75703-29618 \[Eta]) \[Eta] \[Pi])/254016 ,
(*8/2 PN*)8,0.0,
(*9/2 PN*)9,0.0,
(*DEFALUT*)_,0.0
]/.
\[Delta]->\[Delta][\[Eta]]



(* ::Input::Initialization:: *)
N[3715/756+(55 1/4)/9]


(* ::Input::Initialization:: *)
AlphaS[j_,\[Eta]_,\[Chi]1L_,\[Chi]2L_]:=Switch[j,
(*0/2 PN*)0, 0.0,
(*1/2 PN*)1, 0.0,
(*2/2 PN*)2, 0.0,
(*3/2 PN*)3,(*Spin-Orbit*)(113 (\[Chi]1L+\[Chi]2L+\[Chi]1L \[Delta]-\[Chi]2L \[Delta])-76 (\[Chi]1L+\[Chi]2L) \[Eta])/6,
(*4/2 PN*)4,(*Spin-Spin*)(-5 (81 \[Chi]1L^2 (1+\[Delta]-2 \[Eta])+316 \[Chi]1L \[Chi]2L \[Eta]-81 \[Chi]2L^2 (-1+\[Delta]+2 \[Eta])) )/16 ,
(*5/2 PN*)5,0.0,
(*6/2 PN*)6,(*Spin-Orbit*)5 (227 (\[Chi]1L+\[Chi]2L+\[Chi]1L \[Delta]-\[Chi]2L \[Delta])-156 (\[Chi]1L+\[Chi]2L) \[Eta]) \[Pi]/3
+(*Spin-Spin*)5 (20 \[Chi]1L \[Chi]2L \[Eta] (11763+12488 \[Eta])+7 \[Chi]2L^2 (-15103 (-1+\[Delta])+2 (-21683+6580 \[Delta]) \[Eta]-9808 \[Eta]^2)-7 \[Chi]1L^2 (-15103 (1+\[Delta])+2 (21683+6580 \[Delta]) \[Eta]+9808 \[Eta]^2)) /4032,
(*7/2 PN*)7,(*Spin-Orbit*)5 (-5030016755 (\[Chi]1L+\[Chi]2L+\[Chi]1L \[Delta]-\[Chi]2L \[Delta])+4 (2113331119 (\[Chi]1L+\[Chi]2L)+675484362 (\[Chi]1L-\[Chi]2L) \[Delta]) \[Eta]-1008 (208433 (\[Chi]1L+\[Chi]2L)+25011 (\[Chi]1L-\[Chi]2L) \[Delta]) \[Eta]^2+90514368 (\[Chi]1L+\[Chi]2L) \[Eta]^3) /6096384-(*Spin-Spin*)5 (57 \[Chi]1L^2 (1+\[Delta]-2 \[Eta])+220 \[Chi]1L \[Chi]2L \[Eta]-57 \[Chi]2L^2 (-1+\[Delta]+2 \[Eta])) +(*Cubic-in-Spin*)(14585 (\[Chi]1L^3 (1+\[Delta])-\[Chi]2L^3 (-1+\[Delta]))-5 (8819 \[Chi]2L^3+8439 \[Chi]1L \[Chi]2L^2 (-1+\[Delta])-8439 \[Chi]1L^2 \[Chi]2L (1+\[Delta])+\[Chi]1L^3 (8819+2985 \[Delta])) \[Eta]+40 (\[Chi]1L+\[Chi]2L)(17 \[Chi]1L^2-14 \[Chi]1L \[Chi]2L+17 \[Chi]2L^2) \[Eta]^2) /48,
(*8/2 PN*)8, 0.0,
(*9/2 PN*)9,0.0,
(*DEFAULT*)_,0.0
]/.
\[Delta]->\[Delta][\[Eta]]


(* ::Input:: *)
(*N[AlphaNS[3,\[Eta][30,40],-0.3,0.09]*\[Pi]]*)


(* ::Input:: *)
(*AlphaS[3,\[Eta][30,40],-0.3,0.09]*\[Pi]*)


(* ::Input:: *)
(*N[((113*(-0.3+0.09-0.3*\[Delta][\[Eta][30,40]]-0.09*\[Delta][\[Eta][30,40]])-76*(-0.3+0.09)*\[Eta][30,40])/6.)*\[Pi],16]*)


(* ::Input:: *)
(*\[Delta][\[Eta][30,40]]*)


(* ::Input::Initialization:: *)
AlphaLNS[j_,\[Eta]_,\[Chi]1L_,\[Chi]2L_]:=Switch[j,
(*0/2 PN*)0, 0.0,
(*1/2 PN*)1, 0.0,
(*2/2 PN*)2, 0.0,
(*3/2 PN*)3, 0.0,
(*4/2 PN*)4, 0.0,
(*5/2 PN*)5,5 (46374-6552 \[Eta]) \[Pi]/4536,
(*6/2 PN*)6,-6848 /63,
(*7/2 PN*)7,0.0,
(*8/2 PN*)8,0.0,
(*9/2 PN*)9,0.0,
(*DEFAULT*)_,0.0
]/.
\[Delta]->\[Delta][\[Eta]]


(* ::Input::Initialization:: *)
AlphaLS[j_,\[Eta]_,\[Chi]1L_,\[Chi]2L_]:=Switch[j,
(*0/2 PN*)0, 0.0,
(*1/2 PN*)1, 0.0,
(*2/2 PN*)2, 0.0,
(*3/2 PN*)3, 0.0,
(*4/2 PN*)4, 0.0,
(*5/2 PN*)5,(-732985*(\[Chi]1L+\[Chi]2L+\[Chi]1L*\[Delta]-\[Chi]2L*\[Delta])-560*(-1213*(\[Chi]1L+\[Chi]2L)+63*(\[Chi]1L-\[Chi]2L)*\[Delta])*\[Eta]+85680*(\[Chi]1L+\[Chi]2L)*\[Eta]^2)/4536,
(*6/2 PN*)6,0.0,
(*7/2 PN*)7,0.0,
(*8/2 PN*)8,-5 (1263141 (\[Chi]1L+\[Chi]2L+\[Chi]1L \[Delta]-\[Chi]2L \[Delta])-2 (794075 (\[Chi]1L+\[Chi]2L)+178533 (\[Chi]1L-\[Chi]2L) \[Delta]) \[Eta]+94344 (\[Chi]1L+\[Chi]2L) \[Eta]^2) \[Pi]^(11/3)/9072,
(*9/2 PN*)9,0.0,
(*DEFAULT*)_,0.0
]/.
\[Delta]->\[Delta][\[Eta]]


(* ::Input::Initialization:: *)



(* ::Input::Initialization:: *)
Alpha[j_Integer,\[Eta]_,\[Chi]1L_,\[Chi]2L_]:=AlphaNS[j,\[Eta],\[Chi]1L,\[Chi]2L]+AlphaS[j,\[Eta],\[Chi]1L,\[Chi]2L]
AlphaL[j_Integer,\[Eta]_,\[Chi]1L_,\[Chi]2L_]:=AlphaLNS[j,\[Eta],\[Chi]1L,\[Chi]2L]+AlphaLS[j,\[Eta],\[Chi]1L,\[Chi]2L]
Alpha[
jstr_/;(StringMatchQ[jstr,RegularExpression["^[0-9]l$"]]),
\[Eta]_,
\[Chi]1L_/;(-1<\[Chi]1L<1),
\[Chi]2L_/;(-1<\[Chi]2L<1)
]:=AlphaL[ToExpression[StringTake[jstr,1]],\[Eta],\[Chi]1L,\[Chi]2L]


(* ::Input::Initialization:: *)
LALSimPhi[j_Integer,\[Eta]_,\[Chi]1L_,\[Chi]2L_]:=(Alpha[j,\[Eta],\[Chi]1L,\[Chi]2L]+AlphaL[j,\[Eta],\[Chi]1L,\[Chi]2L]Log[\[Pi]])PiOneThirds[j-5]
LALSimPhiL[j_Integer,\[Eta]_,\[Chi]1L_,\[Chi]2L_]:=(AlphaLNS[j,\[Eta],\[Chi]1L,\[Chi]2L])PiOneThirds[j-5]
LALSimPhi[
jstr_/;(StringMatchQ[jstr,RegularExpression["^[0-9]l$"]]),
\[Eta]_,
\[Chi]1L_/;(-1<\[Chi]1L<1),
\[Chi]2L_/;(-1<\[Chi]2L<1)
]:=LALSimPhiL[ToExpression[StringTake[jstr,1]],\[Eta],\[Chi]1L,\[Chi]2L]


(* ::Input::Initialization:: *)
Psi[j_Integer,\[Eta]_,M_,\[Chi]1L_,\[Chi]2L_]:=TimeMetricOneThirds[j-5]CE[\[Eta]](Alpha[j,\[Eta],\[Chi]1L,\[Chi]2L]+AlphaL[j,\[Eta],\[Chi]1L,\[Chi]2L]Log[\[Pi] timeMetric M])PiMOneThirds[j-5,M]
PsiL[j_Integer,\[Eta]_,M_,\[Chi]1L_,\[Chi]2L_]:=TimeMetricOneThirds[j-5]CE[\[Eta]](AlphaL[j,\[Eta],\[Chi]1L,\[Chi]2L])PiMOneThirds[j-5,M]
Psi[
jstr_/;(StringMatchQ[jstr,RegularExpression["^[0-9]l$"]]),
\[Eta]_,
M_,
\[Chi]1L_/;(-1<\[Chi]1L<1),
\[Chi]2L_/;(-1<\[Chi]2L<1)
]:=PsiL[ToExpression[StringTake[jstr,1]],\[Eta],M,\[Chi]1L,\[Chi]2L]


(* ::Input:: *)
(*Psi[3,27,69,0,0]*)


(* ::Input:: *)
(*CE[\[Eta][27,69]]*)


(* ::Input:: *)
(*PiMOneThirds[3-5,27+69]*)


(* ::Section::Initialization:: *)
(*Re-express $\psi$*)


(* ::Input::Initialization:: *)
InsertGetArgs[args_,extras_,extraIndices_]:=Module[
{
(* \:521d\:59cb\:5316args *)
temp=args
},
Do[
(* \:6309\:63d2\:5165\:53c2\:6570\:957f\:5ea6\:ff0c\:6839\:636eindex\:5faa\:73af\:63d2\:5165\:5bf9\:5e94\:53c2\:6570\:ff1b\:6ce8\:610findex\:4e0d\:80fd\:9519\:4f4d *)
temp=Insert[temp,extras[[n]],extraIndices[[n]]],
{n,1,Length[extras]}
];
Sequence@@temp
]


(* ::Input::Initialization:: *)
SolutionToPoint[F_,F0_,gridPoint_,axesIndices_,FArgs_,h0_,hIndex_]:=Module[{
(* \:9664\:4e86\:63d2\:5165\:4e00\:4e9bstatic\:53c2\:6570\:ff0c\:8fd8\:8981\:5c06\:89e3\:65b9\:7a0b\:7684\:53d8\:91cfh\:63d2\:5165\:5230\:51fd\:6570F\:7684\:53c2\:6570\:5217\:8868\:4e2d\:ff1b *)
(* FindRoot\:53ea\:4f1a\:627e\:6700\:8fd1\:7684\:4e00\:4e2a\:6839\:ff0c\:6ce8\:610f\:6839\:7ed3\:6784 *)
solution=FindRoot[{F[Sequence@@Insert[List[InsertGetArgs[FArgs,gridPoint,axesIndices]](*[[1;;3]]*),h,hIndex]]==F0},
{h,h0}]},
(* \:8f93\:51fa\:70b9\:5750\:6807 *)
Flatten[{gridPoint,solution[[1]][[2]]}]
]


(* ::Input::Initialization:: *)

HypersurfaceOverGrid[F_,F0_,axes_,axesIndices_,FArgs_,h0_,hIndex_]:=Module[
{
hypersurface=Switch[Length[axes],
1,
(*
Table[
SolutionToPoint[F,F0,{axes[[1]][[i]]},axesIndices,FArgs,h0,hIndex],
{i,1,Length[axes[[1]]]}
]
*)
FoldList[
Function[{point,i},SolutionToPoint[F,F0,{axes[[1]][[i]]},axesIndices,FArgs,point[[2]],hIndex]],
{0,h0},
Range[Length[axes[[1]]]]
][[2;;]]
,
2,
Table[
SolutionToPoint[F,F0,{axes[[1]][[i]],axes[[2]][[k]]},axesIndices,FArgs,h0,hIndex],
{i,1,Length[axes[[1]]]},{k,1,Length[axes[[2]]]}
]
]
},
hypersurface
]


(* ::Input::Initialization:: *)
HypersurfacePancake[F_,F0s_,axes_,axesIndices_,FArgs_,h0s_,hIndex_]:=Module[
{
pancake=Table[
HypersurfaceOverGrid[F,F0s[[i]],axes,axesIndices,FArgs,h0s[[i]],hIndex],
{i,1,Length[F0s]}
]
},
pancake
]


(* ::Input::Initialization:: *)
(* transform solutions into dots for plot *)
SolutionsToDots[xindex_,y0_,yindex_,yconstraint_,F_,F0_,Fargs_]:=Module[
{
(*solutions=Normal[NSolve[{Normal[Series[F[Sequence@@Insert[Fargs,y,yindex]],{y,y0,4}]]==F0,yconstraint},y,Reals]]*)
solutions=Normal[FindRoot[{F[Sequence@@Insert[Fargs,y,yindex]]==F0},{y,y0}]]
},
(*Table[{Fargs[[xindex]],solutions[[j]][[1]][[2]]},{j,1,Length[solutions]}]*)
Table[{Fargs[[xindex]],solutions[[j]][[2]]},{j,1,Length[solutions]}]
]



(* ::Input::Initialization:: *)
PsiSolutionsToDots[xs_,xindex_,y0_,yindex_,yconstraint_,psij0_,Psiargs_]:=Flatten[
Table[
SolutionsToDots[
xindex,
y0,yindex,yconstraint,
Function[{j,m1,m2,\[Chi]1L,\[Chi]2L},Psi[j,\[Eta][m1,m2],M[m1,m2],\[Chi]1L,\[Chi]2L]],
psij0,List[Sequence@@Insert[Psiargs,x,xindex]]
],
{x,xs}],1]




(* ::Input::Initialization:: *)
AlphaSolutionsToDots[xs_,xindex_,y0_,yindex_,yconstraint_,psij0_,Psiargs_]:=Flatten[
Table[
SolutionsToDots[
xindex,
y0,yindex,yconstraint,
Function[{j,eta,\[Chi]1L,\[Chi]2L},Alpha[j,eta,\[Chi]1L,\[Chi]2L]],
psij0,List[Sequence@@Insert[Psiargs,x,xindex]]
],
{x,xs}],1]



(* ::Input::Initialization:: *)
AlphaASSolutionsToDots[xs_,xindex_,y0_,yindex_,yconstraint_,psij0_,Psiargs_]:=Flatten[
Table[
SolutionsToDots[
xindex,
y0,yindex,yconstraint,
Function[{j,eta,\[Chi]s,\[Chi]a},Alpha[j,eta,\[Chi]1L[\[Chi]s,\[Chi]a],\[Chi]2L[\[Chi]s,\[Chi]a]]],
psij0,List[Sequence@@Insert[Psiargs,x,xindex]]
],
{x,xs}],1]


(* ::Input:: *)
(*PointsFromGrid[xGrid_,yGrid_,hFunc_,h0_,hArgs_,zStartFunc_,zIndex_]:=Table[*)
(*{xGrid[[i]],yGrid[[j]],*)
(*FindRoot[*)
(*{hFunc[]==h0},*)
(*{z,zStartFunc[]}*)
(*][[1]][[2]]*)
(*},*)
(*{i,Length[xGrid]},{j,Length[yGrid]}*)
(*]*)


(* ::Section::Initialization:: *)
(*Plot Functions*)


(* ::Input::Initialization:: *)
PlotParameterWithBCI[parSamples_,parLatex_,BCI_:{0.1,0.5,0.9},bins_:100]:=Module[
{mass1BCI=Quantile[parSamples,{0.1,0.5,0.9}],peak=Max[HistogramList[parSamples,bins,"PDF"][[2]]]},
Print[Histogram[parSamples,bins,"PDF",
PlotLabel->parLatex,
AxesLabel->{parLatex,"Probability Density"},
Epilog->{
Red,
Dashed,
Line[{{mass1BCI[[1]],0},{mass1BCI[[1]],peak}}],(*\:4e0b\:9650\:5782\:7ebf*)
Line[{{mass1BCI[[3]],0},{mass1BCI[[3]],peak}}],(*\:4e0a\:9650\:5782\:7ebf*)
Black,
Dashed,
Line[{{mass1BCI[[2]],0},{mass1BCI[[2]],peak}}],(*Median*)
Black,
(* Text["95% CI",{mass1BCI[[2]],0.12}], *)
Text[Style[NumberForm[mass1BCI[[1]],3],FontSize->10],{mass1BCI[[1]],peak/2}],
Text[Style[NumberForm[mass1BCI[[2]],3],FontSize->10],{mass1BCI[[2]],peak/2}],Text[Style[NumberForm[mass1BCI[[3]],3],FontSize->10],{mass1BCI[[3]],peak/2}]
}
]];
{mass1BCI}
]


(* ::Input::Initialization:: *)
Plot3DCredibleKernel[x_,y_,z_,quantile_:{0.1,0.5,0.9}]:=Module[{
BCI,position
},
BCI["x"]=Quantile[x,quantile];
BCI["y"]=Quantile[y,quantile];
BCI["z"]=Quantile[z,quantile];
position[{"x","left"}]={
{BCI["x"][[1]],BCI["y"][[1]],BCI["z"][[1]]},
{BCI["x"][[1]],BCI["y"][[3]],BCI["z"][[1]]},
{BCI["x"][[1]],BCI["y"][[3]],BCI["z"][[3]]},
{BCI["x"][[1]],BCI["y"][[1]],BCI["z"][[3]]}
};
position[{"x","right"}]={
{BCI["x"][[3]],BCI["y"][[1]],BCI["z"][[1]]},
{BCI["x"][[3]],BCI["y"][[3]],BCI["z"][[1]]},
{BCI["x"][[3]],BCI["y"][[3]],BCI["z"][[3]]},
{BCI["x"][[3]],BCI["y"][[1]],BCI["z"][[3]]}
};
position[{"y","left"}]={
{BCI["x"][[1]],BCI["y"][[1]],BCI["z"][[1]]},
{BCI["x"][[1]],BCI["y"][[1]],BCI["z"][[3]]},
{BCI["x"][[3]],BCI["y"][[1]],BCI["z"][[3]]},
{BCI["x"][[3]],BCI["y"][[1]],BCI["z"][[1]]}
};
position[{"y","right"}]={
{BCI["x"][[1]],BCI["y"][[3]],BCI["z"][[1]]},
{BCI["x"][[1]],BCI["y"][[3]],BCI["z"][[3]]},
{BCI["x"][[3]],BCI["y"][[3]],BCI["z"][[3]]},
{BCI["x"][[3]],BCI["y"][[3]],BCI["z"][[1]]}
};
position[{"z","left"}]={
{BCI["x"][[1]],BCI["y"][[1]],BCI["z"][[1]]},
{BCI["x"][[3]],BCI["y"][[1]],BCI["z"][[1]]},
{BCI["x"][[3]],BCI["y"][[3]],BCI["z"][[1]]},
{BCI["x"][[1]],BCI["y"][[3]],BCI["z"][[1]]}
};
position[{"z","right"}]={
{BCI["x"][[1]],BCI["y"][[1]],BCI["z"][[3]]},
{BCI["x"][[3]],BCI["y"][[1]],BCI["z"][[3]]},
{BCI["x"][[3]],BCI["y"][[3]],BCI["z"][[3]]},
{BCI["x"][[1]],BCI["y"][[3]],BCI["z"][[3]]}
};
Sequence[{
FaceForm[LightRed],EdgeForm[Black],Polygon[position[{"x","left"}]],
FaceForm[LightRed],EdgeForm[Black],Polygon[position[{"x","right"}]],
FaceForm[LightYellow],EdgeForm[Black],Polygon[position[{"y","left"}]],
FaceForm[LightYellow],EdgeForm[Black],Polygon[position[{"y","right"}]],
FaceForm[LightBlue],EdgeForm[Black],Polygon[position[{"z","left"}]],
FaceForm[LightBlue],EdgeForm[Black],Polygon[position[{"z","right"}]]
}]
]


(* ::Section::Initialization:: *)
(*Value Domains*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Manipulate[*)
(*Plot3D[*)
(*Function[{m1,m2},Psi[3,\[Eta][m1,m2],M[m1,m2],spin1z,spin2z]][m1,m2],*)
(*{m1,20,80},*)
(*{m2,20,80},*)
(*PlotPoints->{20,20},*)
(*RegionFunction->Function[{m1,m2,psi},m1>m2],*)
(*ColorFunction->(If[#3<=-800,Red,Blue]&),*)
(*ColorFunctionScaling->False*)
(*],*)
(*(*{spin1z,BCI["spin1_z"]},{spin2z,BCI["spin2_z"]}*)*)
(*{spin1z,-1,1},{spin2z,-1,1}*)
(**)
(*]*)


(* ::Input:: *)
(*Manipulate[*)
(*Plot3D[*)
(*Function[{m1,m2},Psi["5l",\[Eta][m1,m2],M[m1,m2],spin1z,spin2z]][m1,m2],*)
(*{m1,20,80},*)
(*{m2,20,80},*)
(*PlotPoints->{20,20},*)
(*RegionFunction->Function[{m1,m2,psi},m1>m2],*)
(*ColorFunction->(If[#3<=-800,Red,Blue]&),*)
(*ColorFunctionScaling->False*)
(*],*)
(*(*{spin1z,BCI["spin1_z"]},{spin2z,BCI["spin2_z"]}*)*)
(*{spin1z,-1,1},{spin2z,-1,1}*)
(**)
(*]*)


(* ::Input:: *)
(*Manipulate[*)
(*Plot3D[*)
(*Function[{m1,m2},Psi[5,\[Eta][m1,m2],M[m1,m2],spin1z,spin2z]][m1,m2],*)
(*{m1,20,80},*)
(*{m2,20,80},*)
(*PlotPoints->{20,20},*)
(*RegionFunction->Function[{m1,m2,psi},m1>m2],*)
(*ColorFunction->(If[#3<=-800,Red,Blue]&),*)
(*ColorFunctionScaling->False*)
(*],*)
(*(*{spin1z,BCI["spin1_z"]},{spin2z,BCI["spin2_z"]}*)*)
(*{spin1z,-1,1},{spin2z,-1,1}*)
(**)
(*]*)


(* ::Input::Initialization:: *)
EndPackage[];
