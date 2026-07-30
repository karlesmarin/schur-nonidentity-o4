# 🔭 Schur Functions at $(1,-1,t,t^{-1})$

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21463000-1B6F8C?logo=doi&logoColor=white)](https://doi.org/10.5281/zenodo.21463000)
[![License](https://img.shields.io/badge/License-Apache_2.0-B5530F)](LICENSE)
[![Verified](https://img.shields.io/badge/verified-exact_vs_bialternant-2C2C2C)](verify_secondfixed.py)
[![Language](https://img.shields.io/badge/paper-EN_%2B_ES-1B6F8C)](.)

**For every partition $\lambda$ with at most four rows, $s_\lambda(1,-1,t,t^{-1})$ is zero, or $\pm$ a product of exactly three $SU(2)$ characters read off the 2-quotient of $\lambda$ — with no hypothesis on the shape of $\lambda$.**

**📄 Paper (EN + ES), each with its own four plates and their generators, and seven verification runs with their archived output, on Zenodo → https://doi.org/10.5281/zenodo.21463000**

> ### 📚 Part **IV** of a series
> - **Part I — *Anomaly- and Tadpole-Compatible Fermion Completion of 6D SU(4) GHU***
>   → [github.com/karlesmarin/ghu-su4-completion](https://github.com/karlesmarin/ghu-su4-completion) · [Zenodo 10.5281/zenodo.21432625](https://doi.org/10.5281/zenodo.21432625)
> - **Part II — *Three Gates to a Quark Generation***
>   → [github.com/karlesmarin/su4-sm-cell-criterion](https://github.com/karlesmarin/su4-sm-cell-criterion) · [Zenodo 10.5281/zenodo.21432627](https://doi.org/10.5281/zenodo.21432627)
> - **Part III — *A Centre-Charge Selection Rule for the Wilson-Line Potential***
>   → [github.com/karlesmarin/centre-parity-selection](https://github.com/karlesmarin/centre-parity-selection) · [Zenodo 10.5281/zenodo.21438226](https://doi.org/10.5281/zenodo.21438226)
> - **Part IV — *Schur Functions at $(1,-1,t,t^{-1})$*** (this repo)

## 🎯 The closed form

> **Observation (verified, not proved).** For a partition $\lambda$ with at most four rows,
> $$s_\lambda(1,-1,t,t^{-1}) \;=\; 0 \quad\text{or}\quad \pm\,\chi_p(t)\,\chi_q(t)\,\chi_r(t),$$
> where $\chi_k$ is the character of the $(k{+}1)$-dimensional $SU(2)$ irreducible and
> $(p,q,r)$ are three integers read off the **2-quotient** of $\lambda$. There is **no
> hypothesis on the shape of $\lambda$.**

The alphabet $\{1,-1,t,t^{-1}\}$ is the generic semisimple element of the determinant $-1$
component of $O(4)$. It sits between the reciprocal-pair factorisation programme of
Ciucu–Krattenthaler and Ayyer–Behrend and the root-of-unity line of Littlewood, Prasad,
Ayyer–Kumari and Karmakar, and **belongs to neither**. The one unstated hypothesis that kept it
untreated: the block argument needs every inversion-fixed point to contribute a *constant* row
of the bialternant, which $+1$ does and $-1$ does not. Replacing the two fixed-point rows by
even/odd column indicators repairs it at rank one — and not at rank two, which is where we stop.

## 🧭 The actual result is a reduction, not the identity

Part III (this series) classified which $SU(4)$ representations have identically-vanishing
character on the non-identity component of $O(4)$; that classification is exactly the **zero
locus** of the closed form above. But the point is what the identity makes possible:

$$\Phi \;=\; F \circ Q.$$

Conditionally on the closed form, every quantity that survives Part III's leading cancellation
— the whole tower of even moments of the twist-imbalance distribution $\delta(m)$ — factors
through the three integers $Q(\lambda)=(p,q,r)$. An arbitrarily large four-row partition is
**compressed to $(p,q,r)$ with no loss of anything the physics can see**, and weight enumeration
is replaced by explicit polynomial formulae. The normalised moment tower lives in the invariant
ring $\mathbb{Q}[C_p,C_q,C_r]^{S_3}$ of the three $SU(2)$ Casimirs; in particular the variance
is one third of the sum of the three Casimirs.

## ✅ Reproduce every number

```bash
python verify_secondfixed.py
```

Self-contained (plain CPython — no SciPy, SymPy or Sage). It regenerates every displayed number
and checks it against the **bialternant** $\det(x_i^{\mu_j})/\det(x_i^{n-j})$ at $(1,-1,t,1/t)$,
evaluated in exact rational arithmetic at many values of $t$ — a proof of equality per partition,
sharing no code with the routines the closed form was derived from — and ends with a control
built to fail.

```
  PASS eq.(6) closed form = bialternant  --  495 partitions, 0 mismatches
  PASS eq.(7) zero-locus predicate       --  0 mismatches
  PASS eq.(10) M0, M2/M0 = e1/3, M4/M0 Casimir form
  PASS counts 3060 / 4845 / 7905
  PASS rank-two table (3 closed rows)
  PASS block identity det[[A,B],[B,A]] = det(A-B) det(A+B)
  PASS control: shifted partition disagrees (0 of 3)
```

Figures: `python make_fig_secondfixed.py` (families, zero locus) and
`python make_fig_spectrum3d.py` (the twist-imbalance spectrum).

## 📌 Status, stated exactly

The closed form is an **Observation**: verified exhaustively on all 3060 partitions with four
parts $\le 14$ and all 4845 with parts $\le 16$, against the bialternant in exact arithmetic,
with no mismatch — and **not proved**. The zero-locus statement is a **Proposition**, proved
independently in Part III. **Rank two is open**, and the paper says exactly where and why.

**Version 4** withdraws one claim, forces one, and adds one. The rank-two partition offered in
versions 1–3 as a third cause of vanishing is the second cause written with six parts, and that
cause holds at every number of reciprocal pairs. The dichotomy is forced: two of the three factors
can never vanish, so the size-three profile never does. And the sign $\zeta$, which earlier
versions could not interpret, is the sign of the parity index $n_+-n_-$ of the orbifold involution
on $V_\lambda$, whose magnitude is the volume of the box — so the vanishing criterion says the
orbifold projection is exactly balanced, and the reduction $\Phi=F\circ Q$ inverts in closed form.

## 📂 Contents

```
paper/     ghu_secondfixed.tex/.pdf (EN), ghu_secondfixed_es.tex/.pdf (ES), four plates in each
           language, the three figure generators, and outputs/ with every archived run
```

| file | what |
|---|---|
| `verify_secondfixed.py` | reproduces and checks every displayed number vs the bialternant |
| `verify_v4_audit.py` | the counts, the closed form on all 3060, rigidity, and the rank-two rows |
| `verify_v4_inverse.py` | the inversion: three moments return the three indices |
| `verify_zeta_root.py` | the index identity `M₀ = n₊ − n₋` and `(p+1)(q+1)(r+1) = |n₊ − n₋|` |
| `verify_v4_legacy.py` | the four counts that had no archived run before v4, plus the sign rules |
| `check_structure.py` · `check_numbers.py` · `check_layout.py` | the audit pass: references, number backing, page layout |
| `make_fig_secondfixed.py` · `make_fig_spectrum3d.py` · `make_fig_index.py` | the plates; `--es` for the Spanish edition |
| `outputs/` | the saved stdout of every run, so each number in the paper is greppable |

## 📖 Citing

```bibtex
@misc{Marin2026SchurNonIdentity,
  author = {Mar\'in, Carles},
  title  = {Schur Functions at $(1,-1,t,t^{-1})$: A Closed Form on the
            Non-Identity Component of $O(4)$},
  year   = {2026},
  doi    = {10.5281/zenodo.21463000},
  note   = {Part IV of a series}
}
```

## ⚖️ Honesty ledger

**Classical, and not claimed here:** that log-concavity is preserved under convolution (Hoggar
1974) and gives a unimodal spectrum (Stanley 1989); that the coefficient sequence is a threefold
discrete-uniform convolution with polynomial moments, because cumulants of a discrete uniform
are polynomial in its width and add (Neuschel; Bradley–Gupta); that the distribution has a
piecewise-polynomial closed form by inclusion–exclusion (De Moivre; André 1876; Comtet; a
discrete box spline, Dahmen–Micchelli 1988); and the involution character value reproduced along
the way (Karmakar, arXiv:2412.17324).

**What is ours:** the closed form, the fixed-point-pair repair, the identification of the zero
locus as Part III's class, the reduction $\Phi = F \circ Q$, and the Casimir reading of the
moments.

---

Carles Marín · `karlesmarin@gmail.com` · Claude (Anthropic) as AI research assistant
