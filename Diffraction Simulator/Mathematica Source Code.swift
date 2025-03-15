//
//  Mathematica Source Code.swift
//  Diffraction Simulator
//
//  Created by Marco Gonzalez on 3/15/25.
//

/* The following is the source code for the project. THis ocde


(* Based on PhD thesis of McMorran *)
(* Beam parameters *)

w0 = 5*10^-6; (*initial beam width, 30.*10^-6, 0.5*10^-3 *)
r0 = -9.99 *10^20; (*initial radius of wavefront curvature *)
ell0 = w0;(* initial transverse coherence width *)
lambda = 0.003*10^-6(*with Mu: 0.58*10^-9;*)

(* Grating parameters *)

d = 1*10^-6; (*period of grating*)
LT = d^2/lambda; (*Talbot length*)
NZ = Floor[20*10^-3/LT];
nu1 = .5; (*grating 1 open fraction *)
nu2 = .5; (*grating 2 open fraction *)
Z1 = 5*10^-3 ; (* position in mm *)
Z2 = Z1 + NZ*LT; (* position in mm *)
Z3 = Z2 + (Z2 - Z1); (* position in mm *)
X2 = 0;  (*initial offset of G2*)
theta = 0;(*twist between gratings, in degrees*)

(* Limits of the plot *)

Zmin = 1*10^-3;(*1*10^-3;*)
Zmax = Z3 + 0.5*Z1;(*1*10^-3;*)
Xmin = -100*d; (*-400*d*)
Xmax = 100*d; (*1900*d*)

xpoints = 500;(*resolution in the x direction*)
zpoints = 500;(*resolution in the z direction*)

dZ = (Zmax - Zmin)/zpoints;(*step resolution used in computation*)
dX = (Xmax - Xmin)/xpoints;(*step resolution used in computation*)

image = False;(*switch to turn on or off the effects of image charges on an el\
ectron beam*)

xticks = {{1, Zmin}, {1 + 0.2*(Zmax - Zmin)/dZ,
    Zmin + 0.2*(Zmax - Zmin)}, {1 + 0.4*(Zmax - Zmin)/dZ,
    Zmin + 0.4*(Zmax - Zmin)}, {1 + 0.6*(Zmax - Zmin)/dZ,
    Zmin + 0.6*(Zmax - Zmin)} , {1 + 0.8*(Zmax - Zmin)/dZ,
    Zmin + 0.8*(Zmax - Zmin)}, {1 + (Zmax - Zmin)/dZ, Zmax}} ;
(*yticks = \
{{1,Xmin},{1+0.25*(Xmax-Xmin)/dX,Xmin+0.25*(Xmax-Xmin)},{1+0.5*(Xmax-Xmin)/dX,\
Xmin+0.5*(Xmax-Xmin)},{1+0.75*(Xmax-Xmin)/dX,Xmin+0.75*(Xmax-Xmin)},  \
{1+(Xmax-Xmin)/dX,Xmax}} ;*)
yticks = {{1, Xmin}, {1 + 0.25*(Xmax - Xmin)/dX,
    Xmin + 0.25*(Xmax - Xmin)}, {1 + 0.5*(Xmax - Xmin)/dX,
    Xmin + 0.5*(Xmax - Xmin)}, {1 + 0.75*(Xmax - Xmin)/dX,
    Xmin + 0.75*(Xmax - Xmin)},  {1 + (Xmax - Xmin)/dX, Xmax}} ;


Grating Fourier Components


width = nu2*d;
thick = 0;
wedgeangle = 0;
tilt = 0;
ReT = Transpose[{-20 + (# - 1) & /@ Range[41], ConstantArray[0, 41]}];
ImT = ReT;
res = 1000;
nu = width/d;
alpha = wedgeangle*Pi/180;
beta = tilt*Pi/180;
minposX =
  If[beta >= 0, -((width Cos[beta])/2) + width/res,
   If[Abs[beta] <= alpha, -((width Cos[beta])/2) +
     width/res, -((width Cos[beta])/2) + width/res -
     thick*(Tan[alpha] - Tan[beta])]];
maxposX =
  If[beta >= 0,
   If[beta <= alpha, (width Cos[beta])/2 - width/res, (width Cos[beta])/2 -
     width/res + thick*(Tan[alpha] - Tan[beta])], (width Cos[beta])/2 -
    width/res];
x2pnt[arr__, n_] :=
 Position[arr[[All, 1]], n][[1]][[
  1]](*return index where n occurs in the index of ImT or ReT*)

Do[
 fc = 2*Pi*n*ex/d;
 ReT[[x2pnt[ReT, n]]][[2]] += Cos[fc];
 ImT[[x2pnt[ImT, n]]][[2]] += Sin[fc],
 
 {n, -20, 20, 1},
 {ex, minposX, maxposX - width/res, width/res}]

ReT[[All, 2]] /= res;
ImT[[All, 2]] /= res;

MatrixForm[ReT]

ListPlot[ReT, Frame -> True, FrameLabel -> {"index", "ReT"}]
ListPlot[ImT, Frame -> True, FrameLabel -> {"index", "ImT"}]



\!\(
TagBox[
RowBox[{"(", "", GridBox[{
{
RowBox[{"-", "20"}],
RowBox[{"-", "0.0019987369566058987`"}]},
{
RowBox[{"-", "19"}],
RowBox[{"-", "0.040397700198427026`"}]},
{
RowBox[{"-", "18"}],
RowBox[{"-", "0.024337593503117455`"}]},
{
RowBox[{"-", "17"}], "0.029102394792033692`"},
{
RowBox[{"-", "16"}], "0.04663927820419638`"},
{
RowBox[{"-", "15"}],
RowBox[{"-", "0.001999289472640462`"}]},
{
RowBox[{"-", "14"}],
RowBox[{"-", "0.054637860548511226`"}]},
{
RowBox[{"-", "13"}],
RowBox[{"-", "0.0343403530197752`"}]},
{
RowBox[{"-", "12"}], "0.040575707166525754`"},
{
RowBox[{"-", "11"}], "0.0681538421877789`"},
{
RowBox[{"-", "10"}],
RowBox[{"-", "0.0019996841892832213`"}]},
{
RowBox[{"-", "9"}],
RowBox[{"-", "0.08468471769082793`"}]},
{
RowBox[{"-", "8"}],
RowBox[{"-", "0.056836422740706`"}]},
{
RowBox[{"-", "7"}], "0.06842650810990264`"},
{
RowBox[{"-", "6"}], "0.1255030579262745`"},
{
RowBox[{"-", "5"}],
RowBox[{"-", "0.0019999210442038133`"}]},
{
RowBox[{"-", "4"}],
RowBox[{"-", "0.1898135461186309`"}]},
{
RowBox[{"-", "3"}],
RowBox[{"-", "0.15429169925061792`"}]},
{
RowBox[{"-", "2"}], "0.23548689777502743`"},
{
RowBox[{"-", "1"}], "0.7562059069869881`"},
{"0", "0.998`"},
{"1", "0.7562059069869881`"},
{"2", "0.23548689777502743`"},
{"3",
RowBox[{"-", "0.15429169925061792`"}]},
{"4",
RowBox[{"-", "0.1898135461186309`"}]},
{"5",
RowBox[{"-", "0.0019999210442038133`"}]},
{"6", "0.1255030579262745`"},
{"7", "0.06842650810990264`"},
{"8",
RowBox[{"-", "0.056836422740706`"}]},
{"9", "0.08468471769082793`"},
{"10",
RowBox[{"-", "0.0019996841892832213`"}]},
{"11", "0.0681538421877789`"},
{"12", "0.040575707166525754`"},
{"13",
RowBox[{"-", "0.0343403530197752`"}]},
{"14",
RowBox[{"-", "0.054637860548511226`"}]},
{"15",
RowBox[{"-", "0.001999289472640462`"}]},
{"16", "0.04663927820419638`"},
{"17", "0.029102394792033692`"},
{"18",
RowBox[{"-", "0.024337593503117455`"}]},
{"19",
RowBox[{"-", "0.040397700198427026`"}]},
{"20",
RowBox[{"-", "0.0019987369566058987`"}]}
},
GridBoxAlignment->{"Columns" -> {{Center}}, "Rows" -> {{Baseline}}},
GridBoxSpacings->{"Columns" -> {
Offset[0.27999999999999997`], {
Offset[0.7]},
Offset[0.27999999999999997`]}, "Rows" -> {
Offset[0.2], {
Offset[0.4]},
Offset[0.2]}}], "", ")"}],
Function[BoxForm`e$,
MatrixForm[BoxForm`e$]]]\)


izx = ConstantArray[ConstantArray[0, zpoints], xpoints];
fmap = ConstantArray[ConstantArray[0, zpoints], Ceiling[xpoints/2]];
ix = Transpose[{Xmin + (# - 1) ((Xmax - Xmin)/(xpoints - 1)) & /@
     Range[xpoints], ConstantArray[0, xpoints]}];
g1 = Transpose[{Xmin + (# - 1) ((Xmax - Xmin)/(2 xpoints - 1)) & /@
     Range[xpoints*2], ConstantArray[NaN, xpoints*2]}];
g2 = g1;
vis = Transpose[{(Zmin - Z1) + (# - 1) ((Zmax - Zmin)/(zpoints - 1)) & /@
     Range[zpoints], ConstantArray[0, zpoints]}];

f1[x_] :=
  If[((x > d*Round[x/d] - d*nu1/2) && (x < d*Round[x/d] + d*nu1/2)), NaN, Z1];
If[Z1 >= Zmin && Z1 <= Zmax, g1[[All, 2]] = Map[f1, g1[[All, 1]]]];
f2[x_] :=
  If[((x - X2 > d*Round[x/d] - d*nu2/2) && (x - X2 < d*Round[x/d] + d*nu2/2)),
    NaN, Z2];
If[Z2 >= Zmin && Z2 <= Zmax, g2[[All, 2]] = Map[f2, g2[[All, 1]]]];
ListPlot[{g1, g2}]


Calculating GSM parameter at first grating

wz[z_, r0_, ell0_, w0_] :=
  w0*Abs[z/zp[z, r0]]*
   Sqrt[1 + lambda^2 zp[z, r0]^2/(ell0^2 w0^2)];(*GSM beam width*)
zp[z_, v_] := (v z)/(v + z);(*Magnification factor due to wavefront curvature*)
ellz[z_, r0_, ell0_, w0_] :=
  ell0*Abs[z/zp[z, r0]]*
   Sqrt[1 + lambda^2 zp[z, r0]^2/(ell0^2 w0^2)];(*GSM coherence length*)
rz[z_, r0_, ell0_, w0_] := z/(
  1 - zp[z,
     r0]/(z*(1 +
        lambda^2*zp[z, r0]^2/(ell0^2*w0^2))));(*radius of wavefront curvature*)
r1 = rz[Z1, r0, ell0, w0];
ell1 = ellz[Z1, r0, ell0, w0];
w1 = wz[Z1, r0, ell0, w0];


Grating function 0


gp0[z01_, r0_, ell0_, w0_] := Block[
   {w1 = wz[z01, r0, ell0, w0],
    ix = Transpose[{Xmin + (# - 1) ((Xmax - Xmin)/(xpoints - 1)) & /@
        Range[xpoints], ConstantArray[0, xpoints]}]
    },
   
   ix[[All, 2]] = Map[Exp[-Pi*(#/w1)^2] &, ix[[All, 1]]];
   
   Return[ix];
   ];
gp0[Zmin + dZ*50, r0, ell0, w0];
ListPlot[%]

Grating function 1

gp1[z12_, r1_, ell1_, w1_] := Block[
   {ell2 = ellz[z12, r1, ell1, w1],
    w2 = wz[z12, r1, ell1, w1],
    r2 = rz[z12, r1, ell1, w1],
    ix = Transpose[{Xmin + (# - 1) ((Xmax - Xmin)/(xpoints - 1)) & /@
        Range[xpoints], ConstantArray[0, xpoints]}],
    cutoff = 1*10^-3,
    lim = 4,
    dn = 0, dm = 0, coef = 1},
   
   Do[
    dn = n - m;
    dm = (m + n)/2;
    
    If[ Not[image],(*If true, ignore effect of image charge*)
     coef = Sinc[nu1*Pi*n]*Sinc[nu1*Pi*m]*nu1^2,
     coef = ReT[[x2pnt[ReT, n]]][[2]]*ReT[[x2pnt[ReT, m]]][[2]] +
       ImT[[x2pnt[ReT, n]]][[2]]*ImT[[x2pnt[ReT, m]]][[2]]];
    
    coef *= Exp[-Pi*(dn*lambda*z12/(d*ell2))^2];
    
    If[coef >= cutoff,
     ix[[All, 2]] +=
       Map[(coef*Exp[-Pi*((# - dm*lambda*z12/d)/w2)^2]*
           Cos[2*Pi*(dn/d)*(# - dm*lambda*z12/d)*(1 - z12/r2)]) &, ix[[All, 1]]];
     ],
    
    {n, -lim, lim, 1}, {m, -lim, lim, 1}];
   Return[ix];
   ];



gp1[Z1 + LT, r1, ell1, w1];

ListPlot[%, Joined -> True]

Grating function 2

gp2[z12_, z23_, mytheta_, ell1x_, w1x_, r1x_, ell1y_, w1y_, r1y_, X2_] :=
  Block[
   {theta = Pi*mytheta/180,
    d1 = d,
    d2 = d,
    
    z13 = z12 + z23,
    ell3x = ellz[z12 + z23, r1x, ell1x, w1x],
    w3x = wz[z12 + z23, r1x, ell1x, w1x],
    r3x = rz[z12 + z23, r1x, ell1x, w1x],
    
    ell3y = ellz[z12 + z23, r1y, ell1y, w1y],
    w3y = wz[z12 + z23, r1y, ell1y, w1y],
    r3y = rz[z12 + z23, r1y, ell1y, w1y],
    
    ix = Transpose[{Xmin + (# - 1) ((Xmax - Xmin)/(xpoints - 1)) & /@
        Range[xpoints], ConstantArray[0, xpoints]}],
    phix = ix, phi = 0,
    cutoff = 1*10^-3,
    lim = 4,
    coef = 1, dn = 0, dm = 0, m = 0, n = 0, a = 0, b = 0, c = 0, d = 0},
   
   Do[
    dn = n1 - n2;
    n = (n1 + n2)/2;
    dm = m1 - m2;
    m = (m1 + m2)/2;
    a = x2pnt[ReT, m1];
    b = x2pnt[ReT, m2];
    c = x2pnt[ReT, n1];
    d = x2pnt[ReT, n2];
    
    
    If[Not[image],(*If true, ignore effect of image charge*)
     {coef = Sinc[nu1*Pi*m1] + 0 I,
      coef *= Sinc[nu1*Pi*m2 + 0 I]},
     {coef = ReT[[a]][[2]] + ImT[[a]][[2]] I,
      coef *= ReT[[b]][[2]] - ImT[[b]][[2]] I}];
    
    coef *= ReT[[c]][[2]] + ImT[[c]][[2]] I;
    coef *= ReT[[d]][[2]] - ImT[[d]][[2]] I;
    
    (*Accounting for twist dependence of visibility*)
    coef *= Exp[-Pi*(dn*Sin[theta]*lambda*z23/(d2*ell3y))^2];
    coef *= Exp[-Pi*(lambda*z23*(dn*Cos[theta] + dm*z13/z23)/(d1*ell3x))^2];
    
    
    (*Print[ix[[All,1]]];
    ix[[All,2]]+=Map[(Re[coef]*Cos[#[[2]]]-Im[coef]*Sin[#[[2]]])*Exp[-Pi*((#[[
    1]]-(lambda*z23/d1)*(n*Cos[theta]+m*z13/z23))/w3x)^2]&,phi];
    Print[ix]*)
    If[Re[coef] >= cutoff || Im[coef] >= cutoff,
     {
      phi =
       dn*n*(1 - z23/r3x)*Cos[theta]^2 + dn*n*(1 - z23/r3y)*Sin[theta]^2 +
        dn*m*(1 - z13/r3x)*Cos[theta];
      phi += dm*n*(1 - z13/r3x)*Cos[theta] + dm*m*(z13/z23)*(1 - z13/r3x);
      phi *= 2*Pi*lambda*z23/(d1^2); phi -= 2*Pi*dn*X2/d2;
      phix[[All, 2]] =
       Map[(phi - (2*Pi*#/d2)*(dn*Cos[theta]*(1 - z23/r3x) +
              dm*(1 - z13/r3x)) &), phix[[All, 1]]];
      
      ix[[All, 2]] +=
       Map[(Re[coef]*Cos[#[[2]]] - Im[coef]*Sin[#[[2]]])*
          Exp[-Pi*((#[[1]] - (lambda*z23/d1)*(n*Cos[theta] + m*z13/z23))/
               w3x)^2] &, phix]
      }
     ]
    ,
    
    {m1, -lim, lim, 1}, {m2, -lim, lim, 1}, {n1, -lim, lim, 1}, {n2, -lim,
     lim, 1}];
   
   Return[ix]];

Testing:

gp2[Z2 - Z1, Zmin + 200*dZ, theta, ell1, w1, r1, ell1, w1, r1, X2];


ListPlot[%, Joined -> True]


Compute intensity at the position of 3rd grating:


gp2dat = gp2[Z2 - Z1, Z2 - Z1, theta, ell1, w1, r1, ell1, w1, r1, X2];
gp2dat[[All, 1]] *= 10^9;
Dimensions[gp2dat]
gp2min = Min[gp2dat[[All, 2]]];
gp2max = Max[gp2dat[[All, 2]]];
gp2mean = (gp2max + gp2min)/2;
gp2amp = gp2max - gp2mean;
bb := Style[#, {Bold, Larger, Black}] &;
marker = Graphics[];
Show[Plot[
  gp2mean + gp2amp*Sin[2* \[Pi]*x/100 + \[Pi]/2], {x, -200 \[Pi], 200 \[Pi]},
  PlotStyle -> LightGray, Frame -> True, AxesOrigin -> {-610, 0.2},
  FrameLabel -> {bb@" y (nm)", bb@"Intensity(arb. units)"},
  TargetUnits -> {1*10^-9,}],
 ListPlot[gp2dat, PlotMarkers -> {"\[FilledSmallCircle]", 11}, Joined -> True]
 ]
(*Print["gp2 Min, Max=",MinMax[gp2dat[[All,2]]]];*)

Text["gp2 Min, Max, Contrast="]
gp2min
gp2max
(gp2max - gp2min)/(gp2max + gp2min)



Part 9: Compute intensity along z axis


tbefore = TimeUsed[]

Do[
 zpos = Zmin + i*dZ;
 If[zpos > Z2,
  ix = gp2[Z2 - Z1, zpos - Z2, theta, ell1, w1, r1, ell1, w1, r1, X2]
  ,
  If[zpos > Z1,
   ix = gp1[zpos - Z1, r1, ell1, w1]
   ,
   ix = gp0[zpos, r0, ell0, w0]
   ]
  ];
 filter1[x_] := If[x < 0, 0., x];
 (*filter2[x_]:=If[x==NaN,0.,x];*)
 ix[[All, 2]] = Map[filter1, ix[[All, 2]]];
 (*ix[[All,2]]=Map[filter2,ix[[All,2]]];*)
 ix[[All, 2]] /= Max[ix[[All, 2]]];
 izx[[All, i + 1]] = ix[[All, 2]];
 ,
 {i, 0, zpoints - 1, 1}]

tafter = TimeUsed[]
(tafter - tbefore)/60


(*Zoomout Plot COLORs: SunsetColors or LakeColors *)

ticx = 0.02*10^-3;
ticz =  4*10^-3;
scalex = 10^-3;
scalez = 10^-3;

ListDensityPlot[
 {izx}, ColorFunction -> ColorData["SunsetColors"]  , AspectRatio -> 0.75,
 ImageSize -> 420, PlotRange -> All, PlotRangePadding -> .001,
 FrameLabel -> {Switch[scalez, LT , "Z distance [Talbot length]", 1,
    "Z distance [m]",  10^-3, "Z distance [mm]", 10^-6 , "Z distance [um]",
    10^-9, "Z distance [nm]"],
   Switch[scalex, d, "X distance [slit spacing]", 1, "X distance [m]",  10^-3,
    "X distance [mm]", 10^-6 , "X distance [um]", 10^-9, "X distance [nm]"]},
 LabelStyle -> (FontSize -> 12),
 FrameTicks -> {{{#, Round[(Xmin + #*dX)/scalex, 0.0001] } & /@
     Range[0, xpoints, ticx/dX],
    None}, {{#, NumberForm@((Zmin + #*dZ)/scalez)} & /@
     Range[0, zpoints, ticz/dZ], None}}
 ]

*/
