(* Step 4: full background evolution and flow parameters at the pivot *)
gbBackground[V_, xi_, p_, p0_, tmax_] :=
 Module[{dV, dX, ddX, rhs, hd, pd, h0, w0, sol, tend, h, w, f, g, k, l, t},
  dV = D[V, p]; dX = D[xi, p]; ddX = D[xi, {p, 2}];
  rhs[q_?NumericQ, u_?NumericQ, y_?NumericQ] :=
   Module[{vp = dV /. p -> q, x1 = dX /. p -> q, x2 = ddX /. p -> q, aa, bb, cc, r1, r2},
    aa = 2 (1 - 4 x1 u y); bb = -4 x1 y^2; cc = 12 x1 y^2;
    r1 = -u^2 + 4 x2 u^2 y^2 - 4 x1 u y^3; r2 = -3 y u - vp - 12 x1 y^4;
    {(r1 - bb r2), (aa r2 - cc r1)}/(aa - bb cc)];
  hd[q_?NumericQ, u_?NumericQ, y_?NumericQ] := rhs[q, u, y][[1]];
  pd[q_?NumericQ, u_?NumericQ, y_?NumericQ] := rhs[q, u, y][[2]];
  {h0, w0} = {h, w} /. FindRoot[{6 h^2 == w^2 + 2 V + 24 dX w h^3,
       w == -(dV/V + 4/3 dX V) h} /. p -> p0, {h, Sqrt[V/3] /. p -> p0}, {w, -0.01}];
  sol = First@NDSolve[{f'[t] == g[t], g'[t] == pd[f[t], g[t], k[t]],
      k'[t] == hd[f[t], g[t], k[t]], l'[t] == k[t],
      f[0] == p0, g[0] == w0, k[0] == h0, l[0] == 0,
      WhenEvent[hd[f[t], g[t], k[t]] + k[t]^2 < 0, "StopIntegration"]},
     {f, g, k, l}, {t, 0, tmax}, AccuracyGoal -> 10, PrecisionGoal -> 10];
  tend = (f /. sol)["Domain"][[1, 2]];
  <|"phi" -> (f /. sol), "phidot" -> (g /. sol), "H" -> (k /. sol), "lna" -> (l /. sol),
    "tEnd" -> tend, "exit" -> (tend < tmax), "dX" -> Function[q, Evaluate[dX /. p -> q]]|>];

gbFlow[bg_, nstar_] := Module[{f, g, k, l, te, ts, eH, dl, e2, d2, tt, ns},
  {f, g, k, l, te} = bg /@ {"phi", "phidot", "H", "lna", "tEnd"};
  ts = tt /. FindRoot[l[tt] == l[te] - nstar, {tt, te/4, 0, te}];
  eH[u_] := -k'[u]/k[u]^2; dl[u_] := 4 bg["dX"][f[u]] g[u] k[u];
  e2 = eH'[ts]/(eH[ts] k[ts]);
  ns = If[dl[ts] == 0, 1 - 2 eH[ts] - e2,
    d2 = dl'[ts]/(dl[ts] k[ts]);
    1 - 2 eH[ts] - (2 eH[ts] e2 - dl[ts] d2)/(2 eH[ts] - dl[ts])];
  <|"efolds" -> l[te], "phiEnd" -> f[te], "phiStar" -> f[ts], "epsH" -> eH[ts],
    "delta1" -> dl[ts], "ns" -> ns, "r" -> 8 (2 eH[ts] - dl[ts])|>];
