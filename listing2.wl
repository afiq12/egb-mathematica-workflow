(* Step 3: slow-roll observables for arbitrary V(p) and xi(p), kappa = 1 *)
gbSlowRoll[V_, xi_, p_] := Module[{Q, e1, d1, e2, d2, ns},
  Q = D[V, p]/V + 4/3 D[xi, p] V;
  e1 = Q D[V, p]/(2 V); d1 = -4/3 D[xi, p] Q V;
  e2 = -Q D[e1, p]/e1;
  ns = If[PossibleZeroQ[D[xi, p]], 1 - 2 e1 - e2,
    d2 = -Q D[d1, p]/d1; 1 - 2 e1 - (2 e1 e2 - d1 d2)/(2 e1 - d1)];
  <|"Q" -> Q, "eps1" -> e1, "delta1" -> d1, "ns" -> ns, "r" -> 8 (2 e1 - d1)|>];

gbPivot[V_, xi_, p_, nstar_, {pe0_, ps0_}] := Module[{sr, pe, ps, q},
  sr = gbSlowRoll[V, xi, p];
  pe = p /. FindRoot[sr["eps1"] == 1, {p, pe0}];
  ps = q /. FindRoot[NIntegrate[1/sr["Q"], {p, pe, q}] == nstar, {q, ps0},
     Evaluated -> False];
  <|"phiEnd" -> pe, "phiStar" -> ps, "ns" -> (sr["ns"] /. p -> ps),
    "r" -> (sr["r"] /. p -> ps)|>];
