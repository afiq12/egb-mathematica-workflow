(* Step 1: curvature of the flat FLRW metric with lapse n(t) *)
coords = {t, x, y, z};
g = DiagonalMatrix[{-n[t]^2, a[t]^2, a[t]^2, a[t]^2}]; gi = Inverse[g];
Gam = Table[Sum[gi[[i, l]] (D[g[[l, j]], coords[[k]]] + D[g[[l, k]], coords[[j]]]
      - D[g[[j, k]], coords[[l]]]), {l, 4}]/2, {i, 4}, {j, 4}, {k, 4}];
Riem = Table[D[Gam[[i, j, l]], coords[[k]]] - D[Gam[[i, j, k]], coords[[l]]]
    + Sum[Gam[[i, k, m]] Gam[[m, j, l]] - Gam[[i, l, m]] Gam[[m, j, k]], {m, 4}],
   {i, 4}, {j, 4}, {k, 4}, {l, 4}];
Ric = Table[Sum[Riem[[i, j, i, l]], {i, 4}], {j, 4}, {l, 4}];
Rs = Simplify[Sum[gi[[j, l]] Ric[[j, l]], {j, 4}, {l, 4}]];
RicSq = Simplify[Sum[gi[[p, r]] gi[[q, s]] Ric[[p, q]] Ric[[r, s]],
    {p, 4}, {q, 4}, {r, 4}, {s, 4}]];
RiemDn = Table[Sum[g[[i, m]] Riem[[m, j, k, l]], {m, 4}], {i, 4}, {j, 4}, {k, 4}, {l, 4}];
RiemUp = Table[Sum[gi[[j, p]] gi[[k, q]] gi[[l, r]] Riem[[i, p, q, r]],
    {p, 4}, {q, 4}, {r, 4}], {i, 4}, {j, 4}, {k, 4}, {l, 4}];
Kret = Simplify[Total[Flatten[RiemDn RiemUp]]];
GB = Simplify[Rs^2 - 4 RicSq + Kret];
check1 = Simplify[{Rs, GB} /. n -> (1 &)];
check2 = Simplify[n[t] a[t]^3 GB - 8 D[a'[t]^3/n[t]^3, t]];

(* Step 2: minisuperspace Lagrangian and Euler-Lagrange equations *)
L = n[t] a[t]^3 (Rs/2 + ph'[t]^2/(2 n[t]^2) - V[ph[t]] - xi[ph[t]] GB/2);
el[q_] := D[L, q[t]] - D[D[L, q'[t]], t] + D[D[L, q''[t]], {t, 2}];
{eN, eA, eP} = Simplify[{el[n], el[a], el[ph]} /. n -> (1 &)];
H = a'[t]/a[t]; Hd = D[H, t];
xid = D[xi[ph[t]], t]; xidd = D[xi[ph[t]], {t, 2}];
c1 = 6 H^2 - ph'[t]^2 - 2 V[ph[t]] - 24 xid H^3;
c2 = 2 Hd + ph'[t]^2 - 4 xidd H^2 - 4 xid H (2 Hd - H^2);
c3 = ph''[t] + 3 H ph'[t] + V'[ph[t]] + 12 xi'[ph[t]] H^2 (Hd + H^2);
check3 = Simplify[{eN/c1, eP/c3}];
check4 = Simplify[eA /. First@Solve[c1 == 0, V[ph[t]]] /. First@Solve[c2 == 0, a''[t]]];
