#!/usr/bin/env python3
"""verify_secondfixed.py - regenerate and check every displayed number of Part IV.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Part IV states an identity it calls an Observation: for a partition with at most four rows,
s_lambda(1,-1,t,1/t) is +-(a product of three SU(2) characters read off the 2-quotient) or zero.
This script regenerates every number the paper quotes and checks it against an INDEPENDENT
ground truth -- the bialternant det(x_i^{mu_j})/det(x_i^{n-j}), evaluated in exact rational
arithmetic at many values of t. Two Laurent polynomials of degree < N that agree at more than
2N+1 points are equal, so the exact-evaluation check is a proof of equality for each partition,
and it shares no code with the chi-convolution routines the closed form was derived from.

Runs in a few seconds on plain CPython; no SciPy, SymPy or Sage needed. Every check prints PASS
or the failing case, and a control that is built to fail confirms the tests can fail.
"""
from fractions import Fraction as F
from itertools import product, combinations_with_replacement
from collections import defaultdict
from math import comb

# ---------------------------------------------------------------- exact linear algebra
def det(M):
    M = [row[:] for row in M]; n = len(M); d = F(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if M[r][i] != 0), None)
        if piv is None: return F(0)
        if piv != i: M[i], M[piv] = M[piv], M[i]; d = -d
        d *= M[i][i]; inv = M[i][i]
        for r in range(i + 1, n):
            if M[r][i] != 0:
                f = M[r][i] / inv; M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return d

def schur(lam, xs):
    """The bialternant s_lambda(xs), exact, at a numeric alphabet xs."""
    n = len(xs); L = list(lam) + [0] * (n - len(lam)); mu = [L[i] + n - 1 - i for i in range(n)]
    return det([[x**m for m in mu] for x in xs]) / det([[x**(n - 1 - j) for j in range(n)] for x in xs])

def chi_at(k, tv): return F(0) if k < 0 else sum(tv**(k - 2 * i) for i in range(k + 1))

# ---------------------------------------------------------------- the closed form, eq. (6)
def closed_at(lam, tv):
    n = 4; L = list(lam) + [0] * (n - len(lam)); mu = [L[i] + n - 1 - i for i in range(n)]
    E = [j for j in range(n) if mu[j] % 2 == 0]; O = [j for j in range(n) if mu[j] % 2 == 1]
    if not E or not O: return F(0)
    inv = sum(1 for i in range(n) for j in range(i + 1, n) if mu[i] % 2 == 1 and mu[j] % 2 == 0)
    zeta = (-1)**(inv + len(E) * len(O))
    al = sorted([mu[j] // 2 for j in E], reverse=True); be = sorted([(mu[j] - 1) // 2 for j in O], reverse=True)
    part = lambda bs: [bs[i] - (len(bs) - 1 - i) for i in range(len(bs))]
    if len(E) == 2:
        l0, l1 = part(al), part(be); z = sum(l0) - sum(l1) - 1
        if z == 0: return F(0)
        return zeta * (1 if z > 0 else -1) * chi_at(l0[0] - l0[1], tv) * chi_at(l1[0] - l1[1], tv) * chi_at(abs(z) - 1, tv)
    nu = part(al if len(E) == 3 else be)
    return zeta * chi_at(nu[0] - nu[1], tv) * chi_at(nu[1] - nu[2], tv) * chi_at(nu[0] - nu[2] + 1, tv)

TS = [F(2), F(3), F(5), F(7), F(3, 2), F(5, 2), F(5, 3), F(7, 3), F(9, 4), F(11, 5),
      F(4, 3), F(8, 5), F(13, 6), F(6, 5), F(10, 7), F(14, 9)]
ALPH = lambda tv: [F(1), F(-1), tv, 1 / tv]

def check(name, ok, detail=""):
    print(("  PASS " if ok else "  FAIL ") + name + (("  --  " + detail) if detail else ""))
    return ok

def main():
    all_ok = True
    print("Part IV -- verification of every displayed number\n")

    # (1) eq. (6): closed form vs bialternant, exact, on all partitions with parts <= 8
    bad = []; n = 0
    for lam in product(range(9), repeat=4):
        if not all(lam[i] >= lam[i + 1] for i in range(3)): continue
        n += 1
        if any(schur(lam, ALPH(tv)) != closed_at(lam, tv) for tv in TS): bad.append(lam)
    all_ok &= check("eq.(6) closed form = bialternant", not bad, "%d partitions, %d mismatches" % (n, len(bad)))

    # (2) eq. (7): the zero locus is exactly the two branches
    bad2 = 0
    for lam in product(range(11), repeat=4):
        if not all(lam[i] >= lam[i + 1] for i in range(3)): continue
        L = list(lam)
        b1 = all((L[i] - L[i + 1]) % 2 == 1 for i in range(3))
        s = {L[i] + L[3 - i] for i in range(4)}
        b2 = (len(s) == 1 and list(s)[0] % 2 == 1 and sum(1 for x in [L[i] + 3 - i for i in range(4)] if x % 2 == 0) == 2)
        if (all(schur(lam, ALPH(tv)) == 0 for tv in TS)) != (b1 or b2): bad2 += 1
    all_ok &= check("eq.(7) zero-locus predicate (parts<=10)", not bad2, "%d mismatches" % bad2)

    # (3) moments in the SU(2) Casimirs, eq. (10)
    def conv(a, b):
        c = defaultdict(int)
        for i, u in a.items():
            for j, v in b.items(): c[i + j] += u * v
        return dict(c)
    def chi(k): return {} if k < 0 else {k - 2 * i: 1 for i in range(k + 1)}
    def mom(d, e): return sum(m**e * c for m, c in d.items())
    bad3 = 0
    for p, q, r in product(range(11), repeat=3):
        d = conv(conv(chi(p), chi(q)), chi(r))
        M0, M2, M4 = mom(d, 0), mom(d, 2), mom(d, 4)
        C = [p * (p + 2), q * (q + 2), r * (r + 2)]
        e1 = sum(C); e2 = C[0] * C[1] + C[0] * C[2] + C[1] * C[2]; s2 = sum(x * x for x in C)
        if M0 != (p + 1) * (q + 1) * (r + 1) \
           or F(M2, M0) != F(e1, 3) \
           or F(M4, M0) != F(2, 3) * e2 + F(1, 5) * s2 - F(4, 15) * e1:
            bad3 += 1
    all_ok &= check("eq.(10) M0, M2/M0=e1/3, M4/M0 Casimir form", not bad3, "1331 triples, %d mismatches" % bad3)

    # (4) the partition counts quoted in the paper (include the empty partition)
    counts_ok = (comb(18, 4) == 3060 and comb(20, 4) == 4845 and 3060 + 4845 == 7905)
    all_ok &= check("counts 3060 / 4845 / 7905 (=3060+4845)", counts_ok)

    # (5) the rank-two table (six-letter alphabet), exact at sample points
    def alph6(a, b): return [F(1), F(-1), a, 1 / a, b, 1 / b]
    PTS = [(F(2), F(3)), (F(3), F(5)), (F(5, 2), F(7, 3)), (F(4), F(7)), (F(7, 5), F(9, 4))]
    rows = {
        (2, 2, 2, 0, 0, 0): lambda a, b: (a**2 * b**2 + a**2 + a * b + b**2 + 1)**2 / (a * b)**2,
        (2, 2, 1, 1, 0, 0): lambda a, b: -(a**2 - a + 1) * (a**2 + a + 1) * (b**2 - b + 1) * (b**2 + b + 1) / (a * b)**2,
        (3, 2, 2, 1, 1, 0): lambda a, b: F(0),
    }
    rt = [lam for lam, fn in rows.items() if any(schur(lam, alph6(a, b)) != fn(a, b) for a, b in PTS)]
    all_ok &= check("rank-two table (3 closed rows)", not rt, "mismatch " + str(rt) if rt else "")

    # (6) the block identity, det[[A,B],[B,A]] = det(A-B) det(A+B)
    bad6 = 0
    for seed in (1, 2, 3, 5, 7, 11):
        m = 3
        A = [[F((seed * (i + 1) + 3 * j + 1) % 11) for j in range(m)] for i in range(m)]
        B = [[F((seed * (2 * i + j) + 5) % 13) for j in range(m)] for i in range(m)]
        M = [A[i] + B[i] for i in range(m)] + [B[i] + A[i] for i in range(m)]
        AmB = [[A[i][j] - B[i][j] for j in range(m)] for i in range(m)]
        ApB = [[A[i][j] + B[i][j] for j in range(m)] for i in range(m)]
        if det(M) != det(AmB) * det(ApB): bad6 += 1
    all_ok &= check("block identity det[[A,B],[B,A]]=det(A-B)det(A+B)", not bad6)

    # control: the closed form must DISAGREE with a shifted partition -- the tests can fail
    ctrl = sum(1 for lam in [(3, 1, 0, 0), (4, 2, 1, 0), (5, 3, 2, 1)]
               if all(schur(lam, ALPH(tv)) == closed_at((lam[0] + 1,) + lam[1:], tv) for tv in TS))
    all_ok &= check("control: shifted partition disagrees (0 of 3)", ctrl == 0, "agreed in %d of 3" % ctrl)

    print("\n" + ("ALL CHECKS PASS" if all_ok else "*** SOME CHECK FAILED ***"))
    return all_ok


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
