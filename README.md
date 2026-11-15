# egb-mathematica-workflow

Wolfram Language (Mathematica) code for the tutorial note

**"From action to observables: a reproducible Mathematica workflow for Einstein-Gauss-Bonnet inflation"**
by Afiq Agung (Universitas Sulawesi Barat, Indonesia).

Preprint: [doi.org/10.5281/zenodo.24242179](https://doi.org/10.5281/zenodo.24242179)
Code archive: [doi.org/10.5281/zenodo.24239065](https://doi.org/10.5281/zenodo.24239065) (v1.0.0)

The three files take the Einstein-Gauss-Bonnet inflation action to checked background
equations, slow-roll predictions for an arbitrary potential and coupling, and a numerical
test of those predictions. They use only built-in functions and were tested with
Wolfram Language 15.0.1.

## Files

| File | Content |
|---|---|
| `listing1.wl` | Curvature invariants of flat FLRW with a lapse function, minisuperspace Lagrangian, Euler-Lagrange equations, and four symbolic checks (`check1` to `check4`) |
| `listing2.wl` | `gbSlowRoll` and `gbPivot`: slow-roll parameters, spectral index and tensor-to-scalar ratio for any `V` and `xi` |
| `listing3.wl` | `gbBackground` and `gbFlow`: numerical integration of the full background equations and flow parameters at the pivot |

## Quick start

Evaluate `listing2.wl` and `listing3.wl` in a notebook, then run the worked example
(E-model potential with an exponential Gauss-Bonnet coupling, `kappa = 1`, `V0 = 1`):

```wolfram
b = Sqrt[2/3]; pot = (1 - Exp[-b p])^2; be = 0.6;

gbPivot[pot, be Exp[-b p], p, 60, {0.8, 5.5}]
(* ns = 0.96890, r = 0.00255 *)

bg = gbBackground[pot, be Exp[-b p], p, 6.2, 2000];
bg["exit"]        (* True: inflation ends *)
gbFlow[bg, 60]
(* ns = 0.96858, r = 0.00262 *)
```

`gbFlow` may print the message `FindRoot::lstol`. The root is already converged and the
results are not affected.

## Conventions

Signature (-,+,+,+), `kappa^2 = 8 Pi G = 1`, and the action

    S = Integral d^4x Sqrt[-g] [ R/2 - (1/2) (d phi)^2 - V(phi) - (1/2) xi(phi) G ]

where `G` is the Gauss-Bonnet invariant.

## Citation

If you use this code, please cite both the note and the code archive:

> Agung, A. (2026). From action to observables: a reproducible Mathematica workflow for
> Einstein-Gauss-Bonnet inflation. Zenodo. https://doi.org/10.5281/zenodo.24242179

> Agung, A. (2026). egb-mathematica-workflow: Wolfram Language code for
> Einstein-Gauss-Bonnet inflation (Version v1.0.0) [Computer software]. Zenodo.
> https://doi.org/10.5281/zenodo.24239065

See `CITATION.cff` for machine-readable metadata.

## License

MIT License. See `LICENSE`.
