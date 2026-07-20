#!/usr/bin/env python3
"""make_fig_secondfixed.py - the two figures of Part IV.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

fig_zerolocus  the vanishing set over partitions with four parts, coloured by which of the two
               branches kills it: alternating beta-parity, or self-complementary with odd
               constant. The point of the plate is that the two branches tile the zero locus
               with no residue -- every vanishing partition is one colour or the other.
fig_families   where our alphabet sits in the Ciucu-Krattenthaler / Ayyer-Behrend sequence, and
               what changes at the step that was never taken.

Both are regenerated from scratch; nothing is hand-placed.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
from itertools import product

HDR = "#1F4E79"
BLUE = "#3B7DD8"
ORANGE = "#E2711D"
GREY = "#EEF2F7"


def beta(lam):
    return [lam[i] + 4 - i - 1 for i in range(4)]


def branch(lam):
    """0 = does not vanish, 1 = alternating parity, 2 = self-complementary with odd constant."""
    L = list(lam)
    if all((L[i] - L[i + 1]) % 2 == 1 for i in range(3)):
        return 1
    s = {L[i] + L[3 - i] for i in range(4)}
    if len(s) == 1 and list(s)[0] % 2 == 1:
        return 2
    return 0


def fig_zerolocus(path, M=14):
    """Plot by (|lambda|, lambda_1): each cell says whether SOME partition of that shape
    vanishes, and by which branch. Cells with both are split."""
    cells = {}
    for l in product(range(M + 1), repeat=4):
        if not all(l[i] >= l[i + 1] for i in range(3)) or sum(l) == 0:
            continue
        b = branch(l)
        if b == 0:
            continue
        key = (sum(l), l[0])
        cells.setdefault(key, set()).add(b)

    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    for (n, l1), bs in cells.items():
        if bs == {1}:
            ax.add_patch(plt.Rectangle((n - .5, l1 - .5), 1, 1, fc=BLUE, ec="white", lw=.5))
        elif bs == {2}:
            ax.add_patch(plt.Rectangle((n - .5, l1 - .5), 1, 1, fc=ORANGE, ec="white", lw=.5))
        else:
            ax.add_patch(plt.Rectangle((n - .5, l1 - .5), .5, 1, fc=BLUE, ec="white", lw=.5))
            ax.add_patch(plt.Rectangle((n, l1 - .5), .5, 1, fc=ORANGE, ec="white", lw=.5))

    ns = [k[0] for k in cells]
    l1s = [k[1] for k in cells]
    ax.set_xlim(min(ns) - 1, max(ns) + 1)
    ax.set_ylim(min(l1s) - 1, max(l1s) + 1)
    ax.set_xlabel(r"$|\lambda|$")
    ax.set_ylabel(r"$\lambda_1$")
    ax.set_title("Where the character vanishes, and by which branch", color=HDR, fontsize=11)
    ax.legend(handles=[Patch(fc=BLUE, label=r"alternating $\beta$-parity"),
                       Patch(fc=ORANGE, label=r"self-complementary, $c$ odd"),
                       Patch(fc="0.6", label="both occur at this $(|\\lambda|,\\lambda_1)$")],
              loc="upper left", fontsize=8, framealpha=.95)
    ax.set_axisbelow(True)
    ax.grid(color="0.9", lw=.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return len(cells)


def fig_families(path):
    """The sequence of alphabets, and the step nobody took."""
    fig, ax = plt.subplots(figsize=(7.8, 3.5))
    ax.axis("off")

    rows = [
        (r"$(X,\bar X)$", "$+1$", "none", "Ciucu--Krattenthaler 2008", True),
        (r"$(X,\bar X,1)$", "$+1$", r"$+1$: constant row", "Ayyer--Behrend 2018", True),
        (r"$(X,\bar X,1,-1)$", "$-1$", r"$-1$: parity row", "this paper", False),
    ]
    y = 2.4
    ax.text(0.02, y + .75, "alphabet", fontsize=9, color=HDR, weight="bold")
    ax.text(0.30, y + .75, r"$\det$", fontsize=9, color=HDR, weight="bold")
    ax.text(0.42, y + .75, "inversion-fixed points", fontsize=9, color=HDR, weight="bold")
    ax.text(0.78, y + .75, "who did it", fontsize=9, color=HDR, weight="bold")
    for alph, det, fixed, who, done in rows:
        fc = GREY if done else "#FDEDE3"
        ax.add_patch(plt.Rectangle((0.0, y - .28), 1.0, .66, fc=fc, ec="0.75", lw=.7,
                                   transform=ax.transData, clip_on=False))
        ax.text(0.02, y, alph, fontsize=11, va="center")
        ax.text(0.30, y, det, fontsize=10, va="center",
                color="0.25" if done else "#B22222", weight="normal" if done else "bold")
        ax.text(0.42, y, fixed, fontsize=9, va="center")
        ax.text(0.78, y, who, fontsize=9, va="center",
                style="normal" if done else "italic",
                color="0.25" if done else "#B22222")
        y -= .95
    ax.text(0.02, y + .25,
            "the block identity needs every fixed point to give a CONSTANT row;\n"
            r"$-1$ gives $((-1)^{\mu_j})_j$, and the triangularisation stops",
            fontsize=8.5, color="#B22222", va="top")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(y - .1, 3.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


if __name__ == "__main__":
    n = fig_zerolocus("paper/fig_zerolocus.pdf")
    print("fig_zerolocus.pdf   %d occupied cells" % n)
    fig_families("paper/fig_families.pdf")
    print("fig_families.pdf    written")
