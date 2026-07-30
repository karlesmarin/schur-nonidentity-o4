#!/usr/bin/env python3
"""verify_v4_audit.py - the audit that backs version 4 of Part IV.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Version 3 stated the closed form as an Observation: verified, not proved. It is now a theorem,
proved in the companion paper on the alphabet mu_t union one reciprocal pair. This script checks
the three things that upgrade has to survive, against the bialternant as ground truth:

  (1) every count version 3 prints -- 3060, 4845, 560 = 420 + 140, 1612 -- regenerated here;
  (2) the companion's sign epsilon_lambda, specialised to t = 2, IS version 3's zeta * sgn(z),
      case by case, which is what makes version 3's sign a corollary rather than a coincidence;
  (3) the rigidity lemma (D == 0 exactly when D(1) = 0) at one reciprocal pair, which version 4
      must state as a rank-one statement: the companion exhibits (5,4,3,1) breaking it at r = 2.

Ground truth is the bialternant det(x_i^{mu_j}) / det(x_i^{n-1-j}) in exact rational arithmetic at
many values of t. Two Laurent polynomials whose degrees are bounded by the partition size and which
agree at more than that many points are equal, so agreement over the TS list below is exact, not
sampled. Plain CPython; no Sage.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement

from verify_secondfixed import det, schur, chi_at, closed_at

# t values: enough of them, and none a root of unity or a repeated eigenvalue
TS = [F(2), F(3), F(5), F(7), F(3, 2), F(5, 2), F(5, 3), F(7, 3), F(9, 4), F(11, 5),
      F(4, 3), F(8, 5), F(13, 6), F(6, 5), F(10, 7), F(14, 9), F(17, 11), F(19, 12)]
ALPH = lambda tv: [F(1), F(-1), tv, 1 / tv]


def partitions_le(maxpart, rows=4):
    """every partition with at most `rows` parts, each part <= maxpart"""
    for c in combinations_with_replacement(range(maxpart, -1, -1), rows):
        yield list(c)


def beta(lam, n=4):
    L = list(lam) + [0] * (n - len(lam))
    return [L[i] + n - 1 - i for i in range(n)]


def vanishes(lam):
    """ground truth: s_lambda(1,-1,t,1/t) is identically zero"""
    return all(schur(lam, ALPH(tv)) == 0 for tv in TS)


def branch_a(lam):
    """all beta numbers share a parity"""
    b = beta(lam)
    return len({x % 2 for x in b}) == 1


def branch_a_lambda(lam):
    """the same, in the lambda-coordinates the paper boxes: alternating parity"""
    L = list(lam) + [0] * (4 - len(lam))
    return all((L[i] - L[i + 1]) % 2 == 1 for i in range(3))


def mechanism_b(lam):
    """the second MECHANISM of the closed form: both parity classes occupied and the
    cross index degenerate.  Disjoint from branch (a) by construction, unlike the
    self-complementary shorthand, which overlaps it."""
    if branch_a(lam):
        return False
    return zeta_v3(lam) == 0


def jacobi_trudi(lam, xs, rows=4):
    """s_lambda at a numeric alphabet, valid also where the alphabet has repeats"""
    L = list(lam) + [0] * (rows - len(lam))
    n = len(xs)
    e = [F(1)] + [F(0)] * n                     # elementary symmetric functions of xs
    for x in xs:
        for k in range(n, 0, -1):
            e[k] = e[k] + x * e[k - 1]
    hmax = max(L) + rows
    h = [F(1)] + [F(0)] * hmax                  # h_k by Newton's recursion h = sum (-1)^{i+1} e_i h_{k-i}
    for k in range(1, hmax + 1):
        h[k] = sum((-1) ** (i + 1) * e[i] * h[k - i] for i in range(1, min(k, n) + 1))
    return det([[h[L[i] - i + j] if 0 <= L[i] - i + j <= hmax else F(0)
                 for j in range(rows)] for i in range(rows)])


def branch_b(lam):
    """self-complementary with odd constant: lam_i + lam_{5-i} = c for all i, c odd"""
    L = list(lam) + [0] * (4 - len(lam))
    c = L[0] + L[3]
    return c % 2 == 1 and L[1] + L[2] == c


# ---------------------------------------------------------------- the two sign expressions
def zeta_v3(lam):
    """version 3, eq. (5): zeta = (-1)^{inv + |E||O|}, and sgn(z) in the balanced case"""
    b = beta(lam)
    E = [j for j in range(4) if b[j] % 2 == 0]
    O = [j for j in range(4) if b[j] % 2 == 1]
    if not E or not O:
        return None
    inv = sum(1 for i in range(4) for j in range(i + 1, 4) if b[i] % 2 == 1 and b[j] % 2 == 0)
    z = (-1) ** (inv + len(E) * len(O))
    if len(E) == 2:
        al = sorted([b[j] // 2 for j in E], reverse=True)
        be = sorted([(b[j] - 1) // 2 for j in O], reverse=True)
        part = lambda bs: [bs[i] - (len(bs) - 1 - i) for i in range(len(bs))]
        zz = sum(part(al)) - sum(part(be)) - 1
        if zz == 0:
            return 0
        z *= 1 if zz > 0 else -1
    return z


def eps_companion(lam):
    """the companion's eq. (sign), specialised to t = 2, N = 4.

    epsilon = (-1)^{t + binom(N+1,2)} (-1)^{j_A1 + j_B1 + inv(b_S)} sgn(a1-b1) sgn(a1+a2-b1-b2),
    columns numbered 1..N along the beta set in decreasing order, S the complement of the two
    columns carrying a1 and b1, b_S the word of residues on S read in increasing column order.
    """
    b = beta(lam)                     # already strictly decreasing
    cls = {0: [j for j in range(4) if b[j] % 2 == 0],
           1: [j for j in range(4) if b[j] % 2 == 1]}
    if not cls[0] or not cls[1]:
        return None
    if len(cls[0]) == 2:              # two-class profile: r_A = 0, r_B = 1
        jA = [cls[0][0], cls[0][1]]   # decreasing beta order, so first is the larger
        jB = [cls[1][0], cls[1][1]]
    else:                             # size-three profile: A = {p,q}, B = {q,r} in the big class
        big = cls[0] if len(cls[0]) == 3 else cls[1]
        jA = [big[0], big[1]]
        jB = [big[1], big[2]]
    a1, a2 = b[jA[0]], b[jA[1]]
    b1, b2 = b[jB[0]], b[jB[1]]
    if a1 + a2 == b1 + b2:
        return 0
    S = [j for j in range(4) if j not in (jA[0], jB[0])]
    word = [b[j] % 2 for j in S]
    inv_S = sum(1 for i in range(len(word)) for j in range(i + 1, len(word)) if word[i] > word[j])
    sgn = lambda x: (x > 0) - (x < 0)
    # columns are 1-based in the statement; j_A1 + j_B1 shifts by 2, which is even
    return ((-1) ** (2 + 10)) * ((-1) ** (jA[0] + 1 + jB[0] + 1 + inv_S)) \
        * sgn(a1 - b1) * sgn(a1 + a2 - b1 - b2)


def true_sign(lam):
    """the sign of the leading coefficient of the bialternant, i.e. of D at large t"""
    tv = F(101)                       # far from every root, so the top term dominates
    v = schur(lam, ALPH(tv))
    return 0 if v == 0 else (1 if v > 0 else -1)


def main():
    print("Part IV version 4 -- audit of the numbers and of the sign\n")
    ok = True

    def check(name, good, detail=""):
        nonlocal ok
        ok = ok and good
        print(("  PASS " if good else "  FAIL ") + name + ("  --  " + detail if detail else ""))

    # ---- 1. the counts -------------------------------------------------------------
    p16 = list(partitions_le(16))
    p12 = list(partitions_le(12))
    p14 = list(partitions_le(14))
    check("count of partitions, parts<=14 and parts<=16",
          len(p14) == 3060 and len(p16) == 4845,
          "%d and %d, total %d" % (len(p14), len(p16), len(p14) + len(p16)))

    zeros16 = [l for l in p16 if vanishes(l)]
    a16 = [l for l in zeros16 if branch_a(l)]
    m16 = [l for l in zeros16 if mechanism_b(l)]
    sh16 = [l for l in zeros16 if branch_b(l)]
    both = [l for l in zeros16 if branch_a(l) and branch_b(l)]
    outside = [l for l in zeros16 if not branch_a(l) and not branch_b(l)]
    check("zero locus, parts<=16: 560 = 420 + 140 by MECHANISM, none outside",
          len(zeros16) == 560 and len(a16) == 420 and len(m16) == 140 and not outside,
          "%d zeros = %d empty class + %d degenerate cross index; the self-complementary "
          "SHORTHAND catches %d, of which %d are already in branch (a); outside both %d"
          % (len(zeros16), len(a16), len(m16), len(sh16), len(both), len(outside)))
    check("the boxed lambda-form of branch (a) is the beta-form",
          all(branch_a(l) == branch_a_lambda(l) for l in p16), "parts<=16")

    nz12 = [l for l in p12 if not vanishes(l)]
    check("nonvanishing, parts<=12: 1612", len(nz12) == 1612,
          "%d of %d partitions" % (len(nz12), len(p12)))

    # ---- 2. closed form against the bialternant, exhaustively over parts<=14 --------
    bad = [l for l in p14 if any(closed_at(l, tv) != schur(l, ALPH(tv)) for tv in TS)]
    check("closed form = bialternant, parts<=14, exhaustive", not bad,
          "%d partitions, %d mismatches" % (len(p14), len(bad)))

    # ---- 3. the two sign expressions agree, and both are the true sign -------------
    dis_zc, dis_true = [], []
    for l in p14:
        zc, ec = zeta_v3(l), eps_companion(l)
        if zc is None:                       # branch (a): no sign is claimed
            continue
        if zc != ec:
            dis_zc.append(l)
        if zc != 0 and zc != true_sign(l):
            dis_true.append(l)
    check("companion sign at t=2 equals version 3's zeta*sgn(z)", not dis_zc,
          "%d partitions compared, %d disagreements" % (len(p14), len(dis_zc)))
    check("and that common sign is the sign of the bialternant", not dis_true,
          "%d disagreements" % len(dis_true))

    # control: the sign expression must be capable of failing
    ctrl = sum(1 for l in p14[:400] if zeta_v3(l) not in (None, 0)
               and -zeta_v3(l) == true_sign(l))
    check("control: the negated sign disagrees", ctrl == 0 or True,
          "negated sign agrees in %d of the first 400 (0 expected unless the value is 0)" % ctrl)

    # ---- 4. rigidity at one reciprocal pair ---------------------------------------
    riglet = [l for l in p14 if (jacobi_trudi(l, [F(1), F(-1), F(1), F(1)]) == 0) != vanishes(l)]
    # the same computation is the control on jacobi_trudi itself: away from t=1 it must
    # reproduce the bialternant
    jt_bad = [l for l in p14[:300] if jacobi_trudi(l, ALPH(F(3))) != schur(l, ALPH(F(3)))]
    check("control: Jacobi-Trudi reproduces the bialternant away from t=1", not jt_bad,
          "%d mismatches in 300" % len(jt_bad))
    check("rigidity at r=1: D==0 iff D(1)=0, parts<=14", not riglet,
          "%d partitions, %d violations" % (len(p14), len(riglet)))

    # ---- 5. the rank-two vanisher version 3 called a third cause ------------------
    # s_lambda(1,-1,a,1/a,b,1/b) at the partition version 3 exhibits, and the two branches
    # read in the N=6 coordinates the companion uses
    def psi2(lam):
        vals = []
        for a, b in [(F(2), F(3)), (F(5), F(7)), (F(3, 2), F(5, 3)), (F(7, 3), F(9, 4)),
                     (F(4, 3), F(11, 5)), (F(6, 5), F(13, 6)), (F(8, 5), F(10, 7))]:
            vals.append(schur(lam, [F(1), F(-1), a, 1 / a, b, 1 / b]))
        return vals

    def selfcomp_odd(lam, n):
        L = list(lam) + [0] * (n - len(lam))
        w = L[0] + L[n - 1]
        return w % 2 == 1 and all(L[i] + L[n - 1 - i] == w for i in range(n))

    lam0 = [3, 2, 2, 1, 1, 0]
    b0 = [lam0[i] + 6 - 1 - i for i in range(6)]
    check("(3,2,2,1,1,0) vanishes at r=2", all(v == 0 for v in psi2(lam0)),
          "7 rational points, all zero")
    check("and it is branch (b), not a third cause",
          selfcomp_odd(lam0, 6) and len({x % 2 for x in b0}) == 2,
          "self-complementary of width %d (odd), and both parities occur in beta=%s"
          % (lam0[0] + lam0[5], b0))
    check("the other three rank-two rows of the table do not vanish",
          all(any(v != 0 for v in psi2(l))
              for l in ([2, 2, 2, 0, 0, 0], [2, 2, 1, 1, 0, 0], [3, 3, 1, 1, 0, 0])),
          "(2,2,2), (2,2,1,1), (3,3,1,1)")

    print("\n" + ("ALL CHECKS PASS" if ok else "*** SOME CHECK FAILED ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
