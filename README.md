# Schur Functions at (1,−1,t,t⁻¹) — Part IV

[![DOI](https://img.shields.io/badge/DOI-pending-lightgrey?logo=doi&logoColor=white)](https://doi.org/10.5281/zenodo.21438226)

**For every partition λ with at most four rows, `s_λ(1,−1,t,t⁻¹)` is zero, or ± a product of exactly three SU(2) characters read off the 2‑quotient of λ — with no hypothesis on the shape of λ.**

**📄 Paper (EN + ES), the three figures with their generators, and a verification script → Zenodo (DOI pending, links from concept 10.5281/zenodo.21438226).**

## What this is

The companion to Part III (*A Centre‑Charge Selection Rule for the Wilson‑Line Potential*, concept DOI 10.5281/zenodo.21438226). Part III classified which SU(4) representations have identically‑vanishing character on the non‑identity (determinant −1) component of O(4) and exported the underlying symmetric‑function question; Part IV answers it with a closed form, of which that classification is exactly the zero locus.

The point is not the identity but the **reduction it makes possible**: conditionally on the closed form, every quantity that survives Part III's leading cancellation — the whole tower of even moments of the twist‑imbalance distribution — factors through three integers read off the 2‑quotient, `Φ = F ∘ Q`. An arbitrarily large four‑row partition is compressed to `(p,q,r)` with no loss of anything the physics can see, and weight enumeration is replaced by explicit polynomial formulae. The normalised moment tower lives in the invariant ring `ℚ[C_p,C_q,C_r]^{S₃}` of the three SU(2) Casimirs.

## Status, stated exactly

The closed form is an **Observation**: verified exhaustively by two independent implementations on all 3060 partitions with four parts ≤ 14 and all 4845 with parts ≤ 16, against the bialternant in exact arithmetic, with no mismatch — and **not proved**. The zero‑locus statement is a Proposition, proved independently in Part III. Rank two is open, and the paper says exactly where and why.

## Honesty ledger

Classical, and not claimed here: log‑concavity preserved under convolution (Hoggar 1974) giving a unimodal spectrum (Stanley 1989); the coefficient sequence as a threefold discrete‑uniform convolution with polynomial moments (Neuschel; Bradley–Gupta); the piecewise‑polynomial closed form of the distribution by inclusion–exclusion (De Moivre; André 1876; Comtet; a discrete box spline, Dahmen–Micchelli 1988); the involution character value (Karmakar, arXiv:2412.17324). Ours: the closed form, the fixed‑point‑pair repair, the identification of the zero locus, the reduction `Φ = F ∘ Q`, and the Casimir reading of the moments.

## Reproduce every number

```
python verify_secondfixed.py
```

Self‑contained (plain CPython, no SciPy/SymPy/Sage). It regenerates every displayed number and checks it against the bialternant `det(x_i^{μ_j})/det(x_i^{n−j})` at `(1,−1,t,1/t)`, evaluated in exact rational arithmetic at many values of `t` — a proof of equality per partition, sharing no code with the routines the closed form was derived from — and ends with a control built to fail.

The figures regenerate with `python make_fig_secondfixed.py` (families, zero locus) and `python make_fig_spectrum3d.py` (the twist‑imbalance spectrum).

## Contents

| file | what |
|---|---|
| `ghu_secondfixed.pdf`, `ghu_secondfixed_es.pdf` | the paper, English and Spanish |
| `ghu_secondfixed.tex`, `ghu_secondfixed_es.tex` | sources |
| `fig_families.pdf`, `fig_zerolocus.pdf`, `fig_spectrum3d.pdf` | the three figures |
| `make_fig_secondfixed.py`, `make_fig_spectrum3d.py` | figure generators |
| `verify_secondfixed.py` | reproduces and checks every displayed number |

## Citation

Carles Marín, *Schur Functions at (1,−1,t,t⁻¹): A Closed Form on the Non‑Identity Component of O(4)* (Part IV), Zenodo (2026).

## License

Apache‑2.0. Author: Carles Marín (with Claude, Anthropic, as AI research assistant).
