#!/usr/bin/env python3
"""verify_v4_inverse.py - is the moment tower of Part IV invertible?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

Part IV proves Phi = F o Q: everything the physics sees depends on lambda only through the three
integers (p,q,r) read off the 2-quotient.  The question this script settles is the converse --
whether the map can be run backwards, i.e. whether finitely many moments of the twist-imbalance
sequence determine (p,q,r).  Two claims are tested, the first algebraically and the second by
exhaustion:

  (1) M_2/M_0 and M_4/M_0 give e_1 and e_2 of the three Casimirs C_k = k(k+2), in closed form;
  (2) M_6/M_0 has a NON-ZERO coefficient on e_3, so the elementary symmetric functions are
      complete, and therefore (M_0, M_2, M_4, M_6) separates unordered triples.

Claim (2) is what makes the tower invertible: (C_p, C_q, C_r) are then the roots of a cubic whose
coefficients are read off four moments, and k = sqrt(C+1) - 1 recovers each index.  Controls: the
same test with C_k replaced by k(k+3), and injectivity re-tested with M_6 dropped.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement


def chi(k):
    """the coefficient sequence of the SU(2) character, as a dict m -> multiplicity"""
    return {k - 2 * i: 1 for i in range(k + 1)}


def conv(a, b):
    out = {}
    for m, x in a.items():
        for n, y in b.items():
            out[m + n] = out.get(m + n, 0) + x * y
    return out


def moments(p, q, r, upto=6):
    d = conv(conv(chi(p), chi(q)), chi(r))
    return [sum(m ** (2 * j) * v for m, v in d.items()) for j in range(upto // 2 + 1)]


def esym(vals):
    a, b, c = vals
    return a + b + c, a * b + a * c + b * c, a * b * c


def main():
    ok = True

    def check(name, good, detail=""):
        nonlocal ok
        ok = ok and good
        print(("  PASS " if good else "  FAIL ") + name + ("  --  " + detail if detail else ""))

    print("Part IV -- is the moment tower invertible?\n")

    C = lambda k: k * (k + 2)
    TRIPLES = [t for t in combinations_with_replacement(range(0, 25), 3)]

    # ---- (1) e_1 and e_2 in closed form from M_2/M_0 and M_4/M_0 --------------------
    bad1 = []
    for p, q, r in TRIPLES:
        M0, M2, M4, _ = moments(p, q, r)
        e1, e2, _ = esym([C(p), C(q), C(r)])
        got_e1 = F(M2, M0) * 3
        got_e2 = F(15, 4) * F(M4, M0) - F(3, 4) * got_e1 ** 2 + got_e1
        if got_e1 != e1 or got_e2 != e2:
            bad1.append((p, q, r, got_e1, e1, got_e2, e2))
    check("e_1 = 3 M_2/M_0 and e_2 = (15/4)M_4/M_0 - (3/4)e_1^2 + e_1", not bad1,
          "%d triples, %d failures" % (len(TRIPLES), len(bad1)))

    # ---- (2) e_3 needs no new moment: M_0 already carries it -----------------------
    # (k+1)^2 = k(k+2)+1 = C_k+1, so M_0^2 = prod (C_k+1) = 1 + e_1 + e_2 + e_3
    bad0 = []
    for p, q, r in TRIPLES:
        M0, M2, M4, _ = moments(p, q, r)
        e1, e2, e3 = esym([C(p), C(q), C(r)])
        if M0 ** 2 != 1 + e1 + e2 + e3:
            bad0.append((p, q, r))
    check("M_0^2 = 1 + e_1 + e_2 + e_3", not bad0,
          "%d triples, %d failures" % (len(TRIPLES), len(bad0)))

    # ---- (2b) the full inversion, run backwards on every triple --------------------
    def invert(M0, M2, M4):
        """recover the multiset {p,q,r} from the first three even moments"""
        e1 = F(M2, M0) * 3
        e2 = F(15, 4) * F(M4, M0) - F(3, 4) * e1 ** 2 + e1
        e3 = M0 ** 2 - 1 - e1 - e2
        # roots of x^3 - e1 x^2 + e2 x - e3, sought among the integers C_k
        out = []
        for k in range(0, 400):
            x = C(k)
            if x ** 3 - e1 * x ** 2 + e2 * x - e3 == 0:
                out.append(k)
        return out

    badinv = []
    for p, q, r in TRIPLES:
        M0, M2, M4, _ = moments(p, q, r)
        got = invert(M0, M2, M4)
        if set(got) != {p, q, r}:
            badinv.append((p, q, r, got))
    check("the three indices are recovered from (M_0,M_2,M_4), as a multiset", not badinv,
          "%d triples, %d failures" % (len(TRIPLES), len(badinv)))

    # ---- (3) the coefficient of e_3 in M_6/M_0, on a non-degenerate sample ---------
    # fit M_6/M_0 = A e_1^3 + B e_1 e_2 + D e_3 + E e_1^2 + G e_2 + H e_1 + I  (7 unknowns)
    # NOTE: the sample must not be degenerate.  The first 60 entries of TRIPLES all have p = 0,
    # hence e_3 = 0 throughout, and the fit then reports coeff(e_3) = 0 for lack of information.
    FIT = [t for t in TRIPLES if min(t) >= 1][:60]
    HOLD = [t for t in TRIPLES if min(t) >= 1][60:]
    rows, rhs = [], []
    for p, q, r in FIT:
        M0, _, _, M6 = moments(p, q, r)
        e1, e2, e3 = esym([C(p), C(q), C(r)])
        rows.append([e1 ** 3, e1 * e2, e3, e1 ** 2, e2, e1, 1])
        rhs.append(F(M6, M0))

    def solve(rows, rhs):
        n = len(rows[0])
        A = [[F(x) for x in row] + [v] for row, v in zip(rows, rhs)]
        piv = 0
        where = [-1] * n
        for col in range(n):
            sel = next((i for i in range(piv, len(A)) if A[i][col] != 0), None)
            if sel is None:
                continue
            A[piv], A[sel] = A[sel], A[piv]
            f = A[piv][col]
            A[piv] = [x / f for x in A[piv]]
            for i in range(len(A)):
                if i != piv and A[i][col] != 0:
                    g = A[i][col]
                    A[i] = [x - g * y for x, y in zip(A[i], A[piv])]
            where[col] = piv
            piv += 1
        sol = [A[where[c]][n] if where[c] >= 0 else F(0) for c in range(n)]
        # consistency over ALL rows, not only the ones used
        for row, v in zip(rows, rhs):
            if sum(F(x) * s for x, s in zip(row, sol)) != v:
                return None
        return sol

    sol = solve(rows, rhs)
    check("M_6/M_0 is a polynomial in e_1,e_2,e_3", sol is not None,
          "" if sol is None else "coefficients " + ", ".join(str(s) for s in sol))
    if sol is not None:
        check("and its coefficient on e_3 is non-zero", sol[2] != 0, "coeff(e_3) = %s" % sol[2])
        # out of sample
        outs = []
        for p, q, r in HOLD:
            M0, _, _, M6 = moments(p, q, r)
            e1, e2, e3 = esym([C(p), C(q), C(r)])
            v = [e1 ** 3, e1 * e2, e3, e1 ** 2, e2, e1, 1]
            if sum(F(x) * s for x, s in zip(v, sol)) != F(M6, M0):
                outs.append((p, q, r))
        check("the fit holds out of sample", not outs,
              "%d triples held out, %d failures" % (len(HOLD), len(outs)))

    # ---- (3) injectivity of (M_0, M_2, M_4, M_6) on unordered triples ---------------
    seen = {}
    coll = []
    for t in TRIPLES:
        key = tuple(moments(*t))
        if key in seen:
            coll.append((seen[key], t))
        seen[key] = t
    check("(M_0,M_2,M_4,M_6) separates unordered triples", not coll,
          "%d triples, %d collisions" % (len(TRIPLES), len(coll)))

    # control: three moments are not enough to be sure a priori -- report, do not assert
    seen3, coll3 = {}, []
    for t in TRIPLES:
        key = tuple(moments(*t)[:3])
        if key in seen3:
            coll3.append((seen3[key], t))
        seen3[key] = t
    print("  INFO  dropping M_6: %d collisions among %d triples%s"
          % (len(coll3), len(TRIPLES), (" e.g. " + str(coll3[0])) if coll3 else ""))

    # control: the wrong invariant
    Cbad = lambda k: k * (k + 3)
    badrows, badrhs = [], []
    for p, q, r in FIT:
        M0, _, _, M6 = moments(p, q, r)
        e1, e2, e3 = esym([Cbad(p), Cbad(q), Cbad(r)])
        badrows.append([e1 ** 3, e1 * e2, e3, e1 ** 2, e2, e1, 1])
        badrhs.append(F(M6, M0))
    check("control: with C_k = k(k+3) no polynomial fits", solve(badrows, badrhs) is None)

    print("\n" + ("ALL CHECKS PASS" if ok else "*** SOME CHECK FAILED ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
