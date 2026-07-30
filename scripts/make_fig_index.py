#!/usr/bin/env python3
"""make_fig_index.py - the parity-index plate of Part IV, version 4.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Two panels, one statement each, both recomputed here from the bialternant rather than asserted:

  LEFT   every partition with four parts <= 12, plotted by size against its parity index
         n_+ - n_-, the trace of the orbifold involution on V_lambda.  The character vanishes
         identically exactly on the line index = 0, and the plate's content is which partitions
         sit on that line: the two branches of the criterion, and nothing else.

  RIGHT  the same index, drawn as what it is.  For four representative partitions the three
         integers of the 2-quotient are the edges of a box whose VOLUME is |n_+ - n_-|.  The
         partitions are large and their representations are large; the box is small, because it
         measures how far the orbifold projection is from perfectly balanced.

Palette: #2E6DA4 / #C55A11, validated (lightness band, chroma floor, CVD separation, normal-vision
floor, contrast) before use.  Sign is carried by a text label and never by colour alone.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from verify_v4_audit import partitions_le, beta, branch_a, branch_b
from verify_zeta_root import box, schur_at_1y11, peval

ES = "--es" in sys.argv          # Spanish edition of the plate; a lookup, not a source rewrite
T = lambda en, es: es if ES else en
SUF = "_es" if ES else ""

BLUE, ORANGE = "#2E6DA4", "#C55A11"
INK, MUTED, GRID = "#1a1a1a", "#5b6470", "#d8dde3"
FAINT = "#9aa4b0"

# ---------------------------------------------------------------- data
rows = []
for lam in partitions_le(12):
    poly = schur_at_1y11(lam)
    dim, idx = peval(poly, 1), peval(poly, -1)
    bx = box(lam)
    rows.append(dict(lam=tuple(lam), size=sum(lam), dim=dim, idx=idx,
                     bx=bx[0] if bx else None,
                     a=branch_a(lam), b=branch_b(lam)))

zero = [r for r in rows if r["idx"] == 0]
nonzero = [r for r in rows if r["idx"] != 0]
assert all(r["a"] or r["b"] for r in zero), "a zero outside the two branches"

fig = plt.figure(figsize=(6.38, 2.95))
gs = fig.add_gridspec(1, 2, width_ratios=[1.17, 1.0], wspace=0.16,
                      left=0.118, right=0.995, bottom=0.255, top=0.855)

# ---------------------------------------------------------------- LEFT: the balance line
ax = fig.add_subplot(gs[0, 0])
ax.axhline(0, color=INK, lw=1.0, zorder=2)
ax.scatter([r["size"] for r in nonzero], [r["idx"] for r in nonzero],
           s=5, c=FAINT, alpha=.55, linewidths=0, zorder=3)

both = [r for r in zero if r["a"] and r["b"]]
only_a = [r for r in zero if r["a"] and not r["b"]]
only_b = [r for r in zero if r["b"] and not r["a"]]

def counts(grp):
    """how many vanishing partitions of each size -- the marker carries the count"""
    c = {}
    for r in grp:
        c[r["size"]] = c.get(r["size"], 0) + 1
    xs = sorted(c)
    return np.array(xs, float), np.array([c[x] for x in xs], float)

for grp, col, lab in ((only_a, BLUE, T(r"alternating $\beta$-parity",
                                       r"paridad $\beta$ alternante")),
                      (only_b, ORANGE, T("self-complementary, odd width",
                                         "autocomplementaria, anchura impar"))):
    xs, n = counts(grp)
    ax.scatter(xs, np.zeros_like(xs), s=10 + 5.5 * n, c=col, marker="o",
               edgecolors="white", linewidths=.6, zorder=5, label=lab)
xs, n = counts(both)
ax.scatter(xs, np.zeros_like(xs), s=10 + 5.5 * n, facecolors="white", edgecolors=INK,
           linewidths=.9, zorder=6, label=T("both", "ambas"))

ax.set_yscale("symlog", linthresh=1, linscale=.45)
ax.set_xlabel(r"$|\lambda|$", color=INK)
ax.set_ylabel(T(r"parity index $\;n_+-n_-$", r"índice de paridad $\;n_+-n_-$"), color=INK)
ax.set_xlim(-1.5, 49)
ax.set_ylim(-260, 260)
ax.set_yticks([-100, -10, 0, 10, 100])
ax.set_yticklabels(["$-100$", "$-10$", "$0$", "$10$", "$100$"])
ax.grid(True, which="major", color=GRID, lw=.6, zorder=0)
ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.tick_params(colors=MUTED, labelsize=8, length=3)
ax.text(48.2, 0.36, r"$D\equiv 0$", color=INK, fontsize=8.0, ha="right", va="bottom",
        zorder=7, bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=.85))
ax.legend(loc="upper center", bbox_to_anchor=(.5, -0.235), ncol=3, fontsize=7.0,
          frameon=False, handletextpad=.3, columnspacing=1.1, labelcolor=INK,
          scatterpoints=1)
ax.set_title(T(r"the character vanishes exactly where the projection is balanced",
               r"el carácter se anula exactamente donde la proyección está equilibrada"),
             fontsize=8.4 if not ES else 7.9, color=INK, pad=7, loc="left")

# ---------------------------------------------------------------- RIGHT: the box IS the index
ax3 = fig.add_subplot(gs[0, 1], projection="3d")
# (8,6,2,0) and (9,7,3,1) differ by the determinant twist lambda -> lambda+(1^4): same box,
# opposite sign, same dimension.  (12,11,6,0) is the plate's headline: dimension 10920, box 12.
picks = [(12, 11, 6, 0), (8, 6, 2, 0), (9, 7, 3, 1), (8, 5, 2, 0)]
chosen = []
for p in picks:
    r = next(x for x in rows if x["lam"] == p)
    chosen.append(r)

def cuboid(ax, o, d, face, edge):
    x, y, z = o
    dx, dy, dz = d
    v = np.array([[x, y, z], [x + dx, y, z], [x + dx, y + dy, z], [x, y + dy, z],
                  [x, y, z + dz], [x + dx, y, z + dz], [x + dx, y + dy, z + dz],
                  [x, y + dy, z + dz]])
    f = [[v[0], v[1], v[2], v[3]], [v[4], v[5], v[6], v[7]], [v[0], v[1], v[5], v[4]],
         [v[2], v[3], v[7], v[6]], [v[1], v[2], v[6], v[5]], [v[0], v[3], v[7], v[4]]]
    ax.add_collection3d(Poly3DCollection(f, facecolors=face, edgecolors=edge,
                                         linewidths=.7, alpha=.92))

# one hue only: in this panel colour carries nothing, the sign is written out.  Blue and
# orange keep the single meaning they have on the left -- the two branches of the criterion.
SLOT = 7.6
cols = []
for i, r in enumerate(chosen):
    p, q, s = r["bx"]
    d = (p + 1, q + 1, s + 1)
    ox = i * SLOT + (SLOT - d[0]) / 2
    cuboid(ax3, (ox, 0, 0), d, BLUE, "white")
    cols.append((ox + d[0] / 2, d, r))

ax3.set_xlim(0, 4 * SLOT)
ax3.set_ylim(0, 5)
ax3.set_zlim(0, 5.4)
ax3.set_box_aspect((3.3, 1.0, 0.95))
ax3.view_init(elev=16, azim=-60)
ax3.set_axis_off()

# label rows in figure coordinates, one column per box, so nothing can collide
bb = ax3.get_position()
for i, (cx, d, r) in enumerate(cols):
    fx = bb.x0 + bb.width * ((i + 0.5) / 4)
    sign = "+" if r["idx"] > 0 else "-"
    fig.text(fx, .795, r"$%d\!\times\!%d\!\times\!%d$" % d,
             ha="center", va="bottom", fontsize=6.8, color=INK)
    fig.text(fx, .715, r"$=\,%s%d$" % (sign, abs(r["idx"])),
             ha="center", va="bottom", fontsize=6.8, color=INK)
    fig.text(fx, .195, r"$%s$" % str(r["lam"]).replace(" ", "").strip("()").replace(",", "{,}"),
             ha="center", va="bottom", fontsize=6.2, color=MUTED)
    fig.text(fx, .125, r"$%s$" % format(r["dim"], ",").replace(",", "{,}"),
             ha="center", va="bottom", fontsize=6.2, color=MUTED)

fig.text(bb.x0 + bb.width / 2, .915,
         T(r"the index is a box: volume $=|n_+-n_-|$",
           r"el índice es una caja: volumen $=|n_+-n_-|$"),
         ha="center", va="bottom", fontsize=8.0 if not ES else 7.6, color=INK)
fig.text(bb.x0 + bb.width / 2, .045, T(r"$\lambda$ and $\dim V_\lambda$",
                                       r"$\lambda$ y $\dim V_\lambda$"),
         ha="center", va="bottom", fontsize=5.9, color=MUTED, style="italic")

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig_index" + SUF)
fig.savefig(out + ".png", dpi=220)
fig.savefig(out + ".pdf")
print("wrote", out + ".pdf / .png")
print("zeros %d = a-only %d + b-only %d + both %d ; non-zero %d"
      % (len(zero), len(only_a), len(only_b), len(both), len(nonzero)))
print("examples:", [(r["lam"], r["bx"], r["idx"], r["dim"]) for r in chosen])
