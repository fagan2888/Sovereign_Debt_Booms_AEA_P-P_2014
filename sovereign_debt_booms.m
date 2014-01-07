(* ::Package:: *)

(* ::Title:: *)
(*Sovereign Debt Booms in Monetary Unions*)


(* ::Subtitle:: *)
(*Mark Aguiar, Manuel Amador, Emmanuel Farhi and Gita Gopinath *)


(* ::Subsubtitle:: *)
(*AER Papers and Proceedings, 2014*)


(* ::Text:: *)
(*Note: This code has only been tested with Log utility. It may not work more generally. Tested with Mathematica 9.*)


Block[{
		Vaut, bmax, bmaxNoPi, bPi, c, HJB, HJBBorrowingRoot, HJBSavingRoot, Vpart1, 
		Vpart2, Vpart3, bstar, Vsolution, VsolNoPi, figs, V, b, p, pi
	}, 
	Block[{(* Parameters *) 
			\[Chi]=.1547, y =1, \[Rho]=.07, r=.06, \[Psi]=.2, u=Log, piBar=.2
		},

		(* Basic definitions *)
		Vaut =u[(1-\[Chi]) y]/\[Rho]; (* -2.4*)
		bmax= (y - InverseFunction[u][(\[Rho] Vaut + \[Psi] piBar )])/r;
		bmaxNoPi = (y - InverseFunction[u][(\[Rho] Vaut)])/r ;
		bPi = b/.First@FindRoot[u'[y - r b]b==\[Psi], {b,2,  0.0001, 10}];
		c[p_]=InverseFunction[u'][- p];

		(* Defining the HJB root functions *)
		HJB[pi_,p_,V_,b_]=(\[Rho] V-u[c[p]]+\[Psi] pi -p(c[p]+r  b -y) //FullSimplify);

		HJBBorrowingRoot[pi_, V_?NumberQ,b_]:=p/.FindRoot[HJB[pi, p, V, b],
			{p, -10^(-15), -u'[y - r b]}];

		HJBSavingRoot[pi_, V_?NumberQ,b_]:=p/.FindRoot[HJB[pi, p, V, b],
			{p, -u'[y - r b], -10^15}];
 
		(* Solving the HJBs in the different regions *)
		Vpart1 = NDSolveValue[{
				V'[b]==HJBBorrowingRoot[0, V[b],b], 
				V[bPi]==u[y-r bPi]/\[Rho]
			}, 
			V, {b, 0.00001, bPi}
		];

		Vpart2 = NDSolveValue[{
				V'[b]==HJBBorrowingRoot[piBar, V[b],b], 
				V[bmax]==(u[y - r bmax]-\[Psi] piBar)/\[Rho]
			},
			V, {b,  bPi, bmax}
		];

		Vpart3 = NDSolveValue[{
				V'[b]==HJBSavingRoot[piBar, V[b],b], 
				V[bPi]==u[y - r bPi]/\[Rho], 
				WhenEvent[Abs[V'[b]+u'[y - r b ]]<0.01,"StopIntegration"]
			},
			V, {b,  bPi, bmax}
		];

		(* Finding bstar *)
		bstar=b/. FindRoot[Vpart3[b]-Vpart2[b],
			{b,(bPi + bmax)/2,bPi+0.001,bmax-0.001}
		];

		(* Creating the Piecewise value function *)
		Vsolution[b_]:=Piecewise[{
			{Vpart1[b], b<=bPi},
			{Vpart3[b], bstar>b>bPi}}, 
			Vpart2[b]
		];

		(* Constructing the value function with no inflation *) 
		VsolNoPi = NDSolveValue[{
				V'[b]==HJBBorrowingRoot[0, V[b], b], 
				V[bmaxNoPi]==u[y - r bmaxNoPi]/\[Rho]
			},
			V, {b, 0.001, bmaxNoPi}
		];

		(* Results *)
		Print["Vaut -> ", Vaut, ",  bPi -> ", bPi, ",  bstar -> ", 
			bstar, ",  bmax -> ", bmax, ",  bmaxNopi -> ", bmaxNoPi];

		figs = {
			Show[{
					Plot[{y  - r b, c[VsolNoPi'[b]]}, 
						{b, 0.001, bmaxNoPi},
						GridLines->{{bPi, bmax, bstar}, None}, 
						GridLinesStyle->Dashed,  
						PlotStyle->{{Dotted}, 
							{Gray,Dashed,  AbsoluteThickness[3]}}], 
					Plot[c[Vsolution'[b]],
						{b, 0.001, bmax}, 
						PlotStyle->{Black, AbsoluteThickness[3]}]
				},
				PlotRange->All, PlotLabel->"Consumption Policy", AxesOrigin-> {0, .7},
				Epilog-> {
					Text[" \!\(\*SubscriptBox[\(b\), \(\[Pi], \*SubscriptBox[\(\[Psi]\), \(1\)]\)]\)", {bPi, 1.05},{-1,0}],
					Text[" \!\(\*SubsuperscriptBox[\(b\), SubscriptBox[\(\[Psi]\), \(1\)], \(*\)]\)",{bstar, 1.05}, {-1,0}],
					Text[" \!\(\*SubscriptBox[\(b\), \(max, \*SubscriptBox[\(\[Psi]\), \(1\)]\)]\)",{bmax,1.05}, {-1,0}]
				}, 
				ImageSize->500
			], 
			Show[{
					Plot[VsolNoPi[b],
						{b, 0.001, bmaxNoPi}, 
						PlotStyle->{AbsoluteThickness[3], Dashed,  Gray}],
					Plot[ Vsolution[b],
						{b, 0.001, bmax}, PlotStyle->{Black, AbsoluteThickness[3]}, 
						PlotPoints->100 ],
					Plot[Vpart2[b],
						{b, bPi, bstar}, PlotStyle->{Black, Dotted}]
				}, 
				PlotRange->All, GridLines->{{bPi, bmax, bstar}, None}, 
				Epilog-> {
					Text[" \!\(\*SubscriptBox[\(b\), \(\[Pi], \*SubscriptBox[\(\[Psi]\), \(1\)]\)]\)", {bPi, -1},{-1,0}],
					Text[" \!\(\*SubsuperscriptBox[\(b\), SubscriptBox[\(\[Psi]\), \(1\)], \(*\)]\)",{bstar, -1}, {-1,0}],
					Text[" \!\(\*SubscriptBox[\(b\), \(max, \*SubscriptBox[\(\[Psi]\), \(1\)]\)]\)",{bmax,-1}, {-1,0}]
				},
				AxesOrigin->{0,-2.5}, GridLinesStyle->Dashed, 
				PlotLabel->"Value Function", ImageSize->500
		]}
	]
]
