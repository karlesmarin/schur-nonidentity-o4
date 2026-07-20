#!/usr/bin/env python3
"""make_fig_spectrum3d.py - the third figure of Part IV: the spectrum of twist imbalances.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

delta(m), the coefficient of t^m in s_lambda(1,-1,t,1/t), is by the closed form the triple
convolution chi_p * chi_q * chi_r. This plate stacks the aggregate profile over every partition
with four parts, |lambda| <= 20, grouped by s = p+q+r, as a ridgeline surface in three
dimensions. Two things are meant to be visible and both are claims the paper makes: the support
is exactly [-s, s] and never wider, and the profile is a smooth single-peaked tent -- the
convolution's shape, not merely its low moments.

Nothing is hand-placed; the coefficients are expanded exactly from eq. (6) in integer arithmetic.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)
import numpy as np
from itertools import product
from collections import defaultdict

HDR = "#1F4E79"
NPARTS, MAXN = 4, 20


def chi(k):
    """chi_k = t^k + t^(k-2) + ... + t^-k, as {exponent: 1}; chi_{-1} = 0."""
    return {} if k < 0 else {k - 2 * i: 1 for i in range(k + 1)}


def conv(a, b):
    c = defaultdict(int)
    for i, u in a.items():
        for j, v in b.items():
            c[i + j] += u * v
    return dict(c)


def triple(lam):
    """Return (p, q, r) of eq. (6), or None when the character vanishes."""
    mu = [lam[i] + NPARTS - i - 1 for i in range(NPARTS)]
    E = [m // 2 for m in mu if m % 2 == 0]
    O = [(m - 1) // 2 for m in mu if m % 2 == 1]
    if not E or not O:
        return None
    if len(E) == 2:                                    # |E| = |O| = 2
        a, b = sorted(E, reverse=True), sorted(O, reverse=True)
        l0 = [a[0] - 1, a[1]]                          # beta-set -> partition
        l1 = [b[0] - 1, b[1]]
        z = sum(l0) - sum(l1) - 1
        return (l0[0] - l0[1], l1[0] - l1[1], abs(z) - 1)
    nu = sorted(E if len(E) == 3 else O, reverse=True)  # three-part component
    v = [nu[0] - 2, nu[1] - 1, nu[2]]
    return (v[0] - v[1], v[1] - v[2], v[0] - v[2] + 1)


def spectra():
    """Aggregate |delta| profile per s = p + q + r."""
    acc = defaultdict(lambda: defaultdict(int))
    for lam in product(range(MAXN + 1), repeat=NPARTS):
        if sum(lam) > MAXN or not all(lam[i] >= lam[i + 1] for i in range(NPARTS - 1)):
            continue
        pqr = triple(lam)
        if pqr is None or min(pqr) < 0:
            continue
        d = conv(conv(chi(pqr[0]), chi(pqr[1])), chi(pqr[2]))
        s = sum(pqr)
        for m, c in d.items():
            acc[s][m] += abs(c)
    return acc


def fig(path):
    acc = spectra()
    ss = sorted(acc)
    fig = plt.figure(figsize=(7.8, 5.2))
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.get_cmap("plasma")

    # the support claim, checked before anything is drawn: the profile lives exactly on [-s, s]
    for s in ss:
        assert min(acc[s]) == -s and max(acc[s]) == s, "support is not [-s, s] at s = %d" % s

    # each stratum is normalised to unit MASS, not unit height: the widening support then has to
    # pay for itself in amplitude, and the decay of the ridge is the log-concavity being claimed
    # s <= 1 is a one- or two-point stratum whose spike would set the whole vertical scale;
    # it is checked by the assertion above and omitted from the plate only
    ss = [s for s in ss if s > 1]
    S, verts, colors = max(ss), [], []
    for i, s in enumerate(ss):
        ms = sorted(acc[s])
        tot = sum(acc[s].values())
        xs = [ms[0]] + ms + [ms[-1]]
        zs = [0.0] + [acc[s][m] / tot for m in ms] + [0.0]
        verts.append(list(zip(xs, zs)))
        colors.append(cmap(0.10 + 0.80 * i / max(1, len(ss) - 1)))

    poly = PolyCollection(verts, facecolors=colors, edgecolors="white", linewidths=.6, alpha=.92)
    ax.add_collection3d(poly, zs=ss, zdir="y")
    top = max(z for v in verts for _, z in v)
    # the two support edges: the boundary the spectrum is claimed never to cross
    ax.plot([-s for s in ss], ss, [0] * len(ss), color=HDR, lw=1.2, ls="--", zorder=6)
    ax.plot([s for s in ss], ss, [0] * len(ss), color=HDR, lw=1.2, ls="--", zorder=6)

    ax.set_xlim(-S - 1, S + 1)
    ax.set_ylim(S + 1, min(ss) - 1)   # small s at the back
    ax.set_zlim(0, top * 1.05)
    ax.set_xlabel(r"Wilson charge $m$", labelpad=8)
    ax.set_ylabel(r"$s=p+q+r$", labelpad=8)
    ax.set_zlabel(r"$|\delta(m)|$  (unit mass)", labelpad=10)
    ax.set_title("The spectrum of twist imbalances, stratified by the $2$-quotient",
                 color=HDR, pad=6)
    ax.view_init(elev=24, azim=-118)
    ax.set_box_aspect((1.5, 1.15, .75))
    ax.grid(True, alpha=.25)
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_facecolor("white")
        pane.pane.set_alpha(.6)
    # matplotlib mis-measures 3d axis labels, so the margin is set by hand rather than
    # left to bbox_inches="tight", which clips the z label
    fig.subplots_adjust(left=.13, right=.99, top=.95, bottom=.05)
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path, "|", len(ss), "strata, s =", ss[0], "..", ss[-1])


if __name__ == "__main__":
    fig("paper/fig_spectrum3d.pdf")
