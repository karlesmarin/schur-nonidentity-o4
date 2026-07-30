#!/usr/bin/env python3
"""verify_zeta_root.py - what IS the sign zeta of Part IV?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Part IV says of its sign: "a sorting sign of a beta-set deciding the sign of a physical curvature
is either a coincidence of bookkeeping or the shadow of something, and we cannot tell which."
This script tests one candidate root, in the order the questions were asked:

  Q1  zeta multiplies every moment, and delta = zeta * (non-negative), so zeta = sign(M_0).
  Q2  M_0 = sum_m delta(m) = D(1), the character at t = 1.
  Q3  D(1) = s_lambda(1,-1,1,1) = Tr(sigma | V_lambda) with sigma = diag(1,-1,1,1), i.e. the
      PARITY INDEX n_+ - n_- of the orbifold involution on the representation.

If Q3 holds then zeta is not bookkeeping: it is the sign of an index, and the "sorting sign" is
sorting because the Weyl formula's split by residue class IS the block split by the eigenvalues
of sigma.  Each step is checked exactly, and two controls are included: the determinant twist
lambda -> lambda + (1^4), which must flip the sign and preserve the box, and the dimensions
n_+ + n_- = dim V_lambda, which must come out integral and positive.

Also audited here, because it fell out of Q1: Part IV writes the sign as `zeta` in section 3 (the
beta-set statistic) and again as `zeta` in section 5 (the moment formulas), but the total sign of
the closed form is zeta * sgn(z).  If the two are not the same symbol, the M_0 formula is wrong on
half of the balanced branch.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement

from verify_secondfixed import det, closed_at
from verify_v4_audit import beta, zeta_v3, jacobi_trudi, partitions_le


# ---------------------------------------------------------------- polynomials in one variable y
def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pscale(a, c):
    return [c * x for x in a]


def schur_at_1y11(lam, rows=4):
    """s_lambda(1, y, 1, 1) as a polynomial in y, by Jacobi-Trudi over Z[y].

    The coefficient of y^k is the number of basis vectors of V_lambda on which the second
    torus coordinate has weight k, so evaluating at y = -1 gives n_+ - n_- and at y = 1 the
    dimension.  Nothing here knows about the closed form.
    """
    L = list(lam) + [0] * (rows - len(lam))
    hmax = max(L) + rows
    # e_i of the alphabet (1, y, 1, 1)
    e = [[1]]
    for x in ([1], [0, 1], [1], [1]):            # the four letters as polynomials
        new = [[0]] * (len(e) + 1)
        for k in range(len(e) + 1):
            t = e[k] if k < len(e) else [0]
            if k > 0:
                t = padd(t, pmul(x, e[k - 1]))
            new[k] = t
        e = new
    h = [[1]] + [[0]] * hmax
    for k in range(1, hmax + 1):
        acc = [0]
        for i in range(1, min(k, 4) + 1):
            acc = padd(acc, pscale(pmul(e[i], h[k - i]), (-1) ** (i + 1)))
        h[k] = acc
    M = [[h[L[i] - i + j] if 0 <= L[i] - i + j <= hmax else [0] for j in range(rows)]
         for i in range(rows)]
    # 4x4 determinant over Z[y], expanded by permutations (small and exact)
    from itertools import permutations
    total = [0]
    for perm in permutations(range(rows)):
        sgn = 1
        p = list(perm)
        for i in range(rows):
            for j in range(i + 1, rows):
                if p[i] > p[j]:
                    sgn = -sgn
        term = [1]
        for i in range(rows):
            term = pmul(term, M[i][perm[i]])
        total = padd(total, pscale(term, sgn))
    return total


def peval(p, y):
    return sum(c * y ** i for i, c in enumerate(p))


def box(lam):
    """(p,q,r), the three indices of the closed form, and the total sign epsilon"""
    b = beta(lam)
    E = [j for j in range(4) if b[j] % 2 == 0]
    O = [j for j in range(4) if b[j] % 2 == 1]
    if not E or not O:
        return None
    al = sorted([b[j] // 2 for j in E], reverse=True)
    be = sorted([(b[j] - 1) // 2 for j in O], reverse=True)
    part = lambda bs: [bs[i] - (len(bs) - 1 - i) for i in range(len(bs))]
    if len(E) == 2:
        l0, l1 = part(al), part(be)
        z = sum(l0) - sum(l1) - 1
        idx = (l0[0] - l0[1], l1[0] - l1[1], abs(z) - 1)
    else:
        nu = part(al if len(E) == 3 else be)
        idx = (nu[0] - nu[1], nu[1] - nu[2], nu[0] - nu[2] + 1)
    return idx, zeta_v3(lam)


def main():
    ok = True

    def check(name, good, detail=""):
        nonlocal ok
        ok = ok and good
        print(("  PASS " if good else "  FAIL ") + name + ("  --  " + detail if detail else ""))

    print("What is zeta?  A descent, one computation per question.\n")
    P = list(partitions_le(14))

    # ---- Q2: M_0 = D(1) --------------------------------------------------------------
    bad = []
    for lam in P:
        bx = box(lam)
        if bx is None:
            continue
        (p, q, r), eps = bx
        M0_closed = eps * (p + 1) * (q + 1) * (r + 1) if min(p, q, r) >= 0 else 0
        D1 = jacobi_trudi(lam, [F(1), F(-1), F(1), F(1)])
        if F(M0_closed) != D1:
            bad.append((lam, M0_closed, D1))
    check("Q2: M_0 from the closed form equals D(1) = s_lambda(1,-1,1,1)", not bad,
          "%d partitions, %d failures" % (len(P), len(bad)))

    # ---- the notation audit: is section 5's zeta section 3's zeta? --------------------
    mism = []
    for lam in P:
        bx = box(lam)
        if bx is None:
            continue
        (p, q, r), eps = bx
        if eps == 0:
            continue
        b = beta(lam)
        E = [j for j in range(4) if b[j] % 2 == 0]
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if b[i] % 2 == 1 and b[j] % 2 == 0)
        zeta_sec3 = (-1) ** (inv + len(E) * (4 - len(E)))
        if zeta_sec3 != eps:
            mism.append(lam)
    print("  INFO  section 3's zeta differs from the total sign on %d of %d partitions "
          "(the balanced branch with z < 0)" % (len(mism), len(P)))

    # ---- Q3: D(1) is the parity index of sigma = diag(1,-1,1,1) ----------------------
    bad3, examples = [], []
    for lam in P[:400]:
        poly = schur_at_1y11(lam)
        trace = peval(poly, -1)
        dim = peval(poly, 1)
        D1 = jacobi_trudi(lam, [F(1), F(-1), F(1), F(1)])
        nplus = F(dim + trace, 2)
        nminus = F(dim - trace, 2)
        if F(trace) != D1 or nplus < 0 or nminus < 0 or nplus.denominator != 1:
            bad3.append(lam)
        if len(examples) < 3 and trace != 0:
            examples.append((lam, dim, int(nplus), int(nminus), trace))
    check("Q3: D(1) = Tr(sigma|V_lambda) = n_+ - n_-, counted over the weight basis", not bad3,
          "%d partitions, %d failures" % (len(P[:400]), len(bad3)))
    for lam, dim, np_, nm_, tr in examples:
        print("        %-16s dim %-6d n_+ %-6d n_- %-6d  n_+-n_- = %d"
              % (str(tuple(lam)), dim, np_, nm_, tr))

    # ---- control 1: the determinant twist must flip the sign and keep the box --------
    bad4 = []
    for lam in P:
        if max(lam) > 12:
            continue
        tw = [x + 1 for x in lam]
        a, b_ = box(lam), box(tw)
        if a is None or b_ is None:
            continue
        # the shift moves every beta by one, so the two parity classes swap: the box is
        # preserved as a MULTISET (d_1 and d_2 trade places), not as an ordered triple
        if sorted(a[0]) != sorted(b_[0]) or a[1] != -b_[1]:
            bad4.append((lam, a, b_))
    check("control: lambda -> lambda+(1^4) permutes the box and flips the sign", not bad4,
          "%d failures" % len(bad4))

    # ---- control 1b: and the value itself is exactly negated ------------------------
    bad4b = []
    for lam in P[:300]:
        d = jacobi_trudi(lam, [F(1), F(-1), F(3), F(1) / 3])
        dt = jacobi_trudi([x + 1 for x in lam], [F(1), F(-1), F(3), F(1) / 3])
        if dt != -d:
            bad4b.append(lam)
    check("control: D_{lambda+(1^4)} = -D_lambda, i.e. the twist is det(A) = -1", not bad4b,
          "%d failures" % len(bad4b))

    # ---- the identity the descent produces ------------------------------------------
    bad6 = []
    for lam in P:
        bx = box(lam)
        if bx is None:
            continue
        (p, q, r), eps = bx
        D1 = jacobi_trudi(lam, [F(1), F(-1), F(1), F(1)])
        if abs(D1) != (p + 1) * (q + 1) * (r + 1) * (0 if eps == 0 else 1):
            bad6.append((lam, D1, (p, q, r)))
    check("THE IDENTITY: (p+1)(q+1)(r+1) = |n_+ - n_-|, the box volume is the parity index",
          not bad6, "%d partitions, %d failures" % (len(P), len(bad6)))

    # ---- control 2: the sign is a property of the element, not of our ordering -------
    bad5 = []
    for lam in P[:300]:
        d1 = jacobi_trudi(lam, [F(1), F(-1), F(1), F(1)])
        d2 = jacobi_trudi(lam, [F(-1), F(1), F(1), F(1)])
        d3 = jacobi_trudi(lam, [F(1), F(1), F(-1), F(1)])
        if not (d1 == d2 == d3):
            bad5.append(lam)
    check("control: the index does not depend on which coordinate carries the -1", not bad5,
          "%d failures" % len(bad5))

    print("\n" + ("ALL CHECKS PASS" if ok else "*** SOME CHECK FAILED ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
