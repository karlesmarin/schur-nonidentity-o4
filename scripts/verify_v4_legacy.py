#!/usr/bin/env python3
"""verify_v4_legacy.py - the four numbers Part IV has printed since version 1 with no archived run.

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

The audit rule adopted for the companion paper is that a number in the text must be locatable in an
archived output.  Four of Part IV's survived three versions without one, because the scripts that
produced them were run and not saved:

  105   the negative control: s_lambda(1,-1,t,1/t) is NOT +- a constant times
        s_{lambda^(0)}(1,t,1/t), so the free-t object is not the centraliser character whose
        dimension Karmakar's Thm 1B computes
  729   the moment identity M_{2r}/M_0 in the invariant ring, solved for r = 1..8 on p,q,r <= 8
  216   the Casimir forms checked by symbolic (t d/dt)^{2r} at t = 1, with two controls
  3375  |delta| symmetric, gap-free on its parity class, log-concave and unimodal

Each is recomputed here from the definitions.  Where the range that produces the printed count was
not recorded, the range is stated and the count reported, so the text can be corrected to match the
run rather than the run bent to match the text.
"""
from fractions import Fraction as F
from itertools import combinations_with_replacement, product

from verify_secondfixed import schur
from verify_v4_audit import partitions_le, jacobi_trudi, vanishes, ALPH, TS


def chi_seq(k):
    return {k - 2 * i: 1 for i in range(k + 1)}


def conv(a, b):
    out = {}
    for m, x in a.items():
        for n, y in b.items():
            out[m + n] = out.get(m + n, 0) + x * y
    return out


def delta(p, q, r):
    return conv(conv(chi_seq(p), chi_seq(q)), chi_seq(r))


def quotient0(lam):
    """lambda^(0): the even-beta component of the 2-quotient, as a partition"""
    b = [lam[i] + 4 - 1 - i for i in range(4)]
    E = sorted([x // 2 for x in b if x % 2 == 0], reverse=True)
    if not E:
        return None
    return [E[i] - (len(E) - 1 - i) for i in range(len(E))]


def main():
    ok = True

    def check(name, good, detail=""):
        nonlocal ok
        ok = ok and good
        print(("  PASS " if good else "  FAIL ") + name + ("  --  " + detail if detail else ""))

    print("Part IV -- the four legacy numbers, recomputed\n")

    # ---- 105: the free-t object is not the centraliser character --------------------
    # the range was not recorded; report it for every plausible bound so the text can name one
    # The count printed since version 1 is 105 and no range reproduces it (parts <= 5,6,7,8 give
    # 110, 186, 290, 439).  The claim itself survives: the only alphabet-wide proportional case in
    # any range tried is the EMPTY partition, where both sides are the constant 1, and that is not
    # a test case.  Version 4 therefore states the range and the count this run produces.
    print("  [105->109] s_lambda(1,-1,t,1/t) vs +- const * s_{lambda^(0)}(1,t,1/t)")
    for M in (5, 6, 7, 8):
        cases = matches = 0
        for lam in partitions_le(M):
            if vanishes(lam) or sum(lam) == 0:
                continue
            q0 = quotient0(lam)
            if q0 is None:
                continue
            cases += 1
            # compare at several t: proportional means one constant works at every point
            ratios = []
            good = True
            for tv in (F(2), F(3), F(5), F(7, 3), F(9, 4)):
                lhs = schur(lam, ALPH(tv))
                rhs = schur(q0 + [0] * (3 - len(q0)), [F(1), tv, 1 / tv])
                if rhs == 0:
                    good = False
                    break
                ratios.append(lhs / rhs)
            if good and ratios and all(x == ratios[0] for x in ratios) and abs(ratios[0]) != 0:
                matches += 1
        print("        parts<=%d : %3d non-vanishing non-trivial cases, %d proportional"
              % (M, cases, matches))
        if M == 5:
            check("109 cases (four rows, parts <= 5, trivial rep excluded), 0 proportional",
                  cases == 109 and matches == 0, "%d cases, %d proportional" % (cases, matches))

    # ---- 729: the moment tower lies in Q[e1,e2,e3], r = 1..8 -----------------------
    TRI8 = list(product(range(9), repeat=3))
    print("\n  [729] moment identity over p,q,r <= 8: %d ordered triples" % len(TRI8))
    C = lambda k: k * (k + 2)

    def esym(v):
        a, b, c = v
        return a + b + c, a * b + a * c + b * c, a * b * c

    def moments(p, q, r, upto):
        d = delta(p, q, r)
        return [sum(m ** (2 * j) * v for m, v in d.items()) for j in range(upto + 1)]

    def basis(e1, e2, e3, deg):
        """monomials e1^a e2^b e3^c with a+2b+3c <= deg, weighted degree"""
        out = []
        for a in range(deg + 1):
            for b in range(deg + 1):
                for c in range(deg + 1):
                    if a + 2 * b + 3 * c <= deg:
                        out.append(e1 ** a * e2 ** b * e3 ** c)
        return out

    def solve(rows, rhs):
        n = len(rows[0])
        A = [[F(x) for x in row] + [F(v)] for row, v in zip(rows, rhs)]
        piv, where = 0, [-1] * n
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
        for row, v in zip(rows, rhs):
            if sum(F(x) * s for x, s in zip(row, sol)) != v:
                return None
        return sol

    allr = True
    for r in range(1, 9):
        rows, rhs, lead = [], [], None
        for (p, q, s) in TRI8:
            M = moments(p, q, s, r)
            e1, e2, e3 = esym([C(p), C(q), C(s)])
            rows.append(basis(e1, e2, e3, r))
            rhs.append(F(M[r], M[0]))
        sol = solve(rows, rhs)
        # the coefficient of e1^r is the first basis entry with a=r,b=c=0
        idx = [i for i, (a, b, c) in enumerate((a, b, c) for a in range(r + 1)
               for b in range(r + 1) for c in range(r + 1) if a + 2 * b + 3 * c <= r)
               if (a, b, c) == (r, 0, 0)]
        good = sol is not None and (not idx or sol[idx[0]] == F(1, 2 * r + 1))
        allr = allr and good
        print("        r=%d : %s%s" % (r, "in the ring" if sol is not None else "NO FIT",
              "" if not idx or sol is None else ", coeff(e1^%d) = %s" % (r, sol[idx[0]])))
    check("729 triples, r = 1..8: the tower lies in Q[e1,e2,e3] and coeff(e1^r) = 1/(2r+1)", allr)

    # ---- 216: the Casimir forms by a second route, with the two controls -----------
    TRI5 = list(product(range(6), repeat=3))
    print("\n  [216] Casimir forms over p,q,r <= 5: %d ordered triples" % len(TRI5))
    agree = shifted = perturbed = 0
    for (p, q, r) in TRI5:
        M = moments(p, q, r, 2)
        e1, e2, _ = esym([C(p), C(q), C(r)])
        m2 = F(M[1], M[0])
        m4 = F(M[2], M[0])
        want4 = F(2, 3) * e2 + F(1, 5) * (C(p) ** 2 + C(q) ** 2 + C(r) ** 2) - F(4, 15) * e1
        if m2 == F(e1, 3) and m4 == want4:
            agree += 1
        e1s, e2s, _ = esym([C(p + 1), C(q), C(r)])       # control: shift one index
        if m2 == F(e1s, 3) and m4 == F(2, 3) * e2s + F(1, 5) * (
                C(p + 1) ** 2 + C(q) ** 2 + C(r) ** 2) - F(4, 15) * e1s:
            shifted += 1
        bad4 = F(2, 3) * e2 + F(1, 5) * (C(p) ** 2 + C(q) ** 2 + C(r) ** 2) - F(3, 15) * e1
        if m2 == F(e1, 3) and m4 == bad4:
            perturbed += 1
    check("216 triples agree with the Casimir forms", agree == len(TRI5),
          "%d of %d" % (agree, len(TRI5)))
    check("control: shifting one index agrees in 0 of 216", shifted == 0, "%d" % shifted)
    check("control: 4/15 -> 3/15 agrees in 1 of 216", perturbed == 1, "%d" % perturbed)

    # ---- 3375: the four shape properties of |delta| ---------------------------------
    TRI14 = list(product(range(15), repeat=3))
    print("\n  [3375] shape of |delta| over p,q,r <= 14: %d ordered triples" % len(TRI14))
    bad = {"symmetric": 0, "gap-free": 0, "log-concave": 0, "unimodal": 0}
    for (p, q, r) in TRI14:
        d = delta(p, q, r)
        ms = sorted(d)
        v = [d[m] for m in ms]
        if any(d[m] != d[-m] for m in ms):
            bad["symmetric"] += 1
        if any(ms[i + 1] - ms[i] != 2 for i in range(len(ms) - 1)):
            bad["gap-free"] += 1
        if any(v[i] * v[i] < v[i - 1] * v[i + 1] for i in range(1, len(v) - 1)):
            bad["log-concave"] += 1
        top = max(v)
        first, last = v.index(top), len(v) - 1 - v[::-1].index(top)
        if any(v[i] > v[i + 1] for i in range(first)) or \
           any(v[i] < v[i + 1] for i in range(last, len(v) - 1)):
            bad["unimodal"] += 1
    check("3375 triples: |delta| symmetric, gap-free, log-concave, unimodal",
          all(x == 0 for x in bad.values()), str(bad))

    # ---- section 8: the four sign representatives, and the (a,a,b,b) rule ----------
    def sign_of(lam):
        v = schur(list(lam), ALPH(F(101)))
        return 0 if v == 0 else (1 if v > 0 else -1)

    reps = {(6, 6, 0, 0): +1, (4, 4, 2, 2): +1, (5, 5, 1, 1): -1, (3, 3, 3, 3): -1}
    got = {k: sign_of(k) for k in reps}
    check("section 8: (6,6,0,0),(4,4,2,2) give + and (5,5,1,1),(3,3,3,3) give -",
          got == reps, str({k: ("+" if v > 0 else "-") for k, v in got.items()}))
    bad_rule = [(a, a, b, b) for a in range(13) for b in range(a + 1)
                if sign_of((a, a, b, b)) not in (0,) and sign_of((a, a, b, b)) != (-1) ** a]
    check("section 8: on (a,a,b,b) the sign is (-1)^{lambda_1}, a <= 12", not bad_rule,
          "%d exceptions" % len(bad_rule))

    # ---- section 5: the control that agrees in 0 of 1612 ---------------------------
    from verify_secondfixed import closed_at, chi_at
    from verify_v4_audit import zeta_v3

    def closed_shifted(lam, tv):
        """the closed form with the first index shifted by one -- the control that must fail"""
        b = [lam[i] + 4 - 1 - i for i in range(4)]
        E = [j for j in range(4) if b[j] % 2 == 0]
        O = [j for j in range(4) if b[j] % 2 == 1]
        if not E or not O:
            return F(0)
        eps = zeta_v3(lam)
        al = sorted([b[j] // 2 for j in E], reverse=True)
        be = sorted([(b[j] - 1) // 2 for j in O], reverse=True)
        part = lambda bs: [bs[i] - (len(bs) - 1 - i) for i in range(len(bs))]
        if len(E) == 2:
            l0, l1 = part(al), part(be)
            z = sum(l0) - sum(l1) - 1
            idx = (l0[0] - l0[1] + 1, l1[0] - l1[1], abs(z) - 1)
        else:
            nu = part(al if len(E) == 3 else be)
            idx = (nu[0] - nu[1] + 1, nu[1] - nu[2], nu[0] - nu[2] + 1)
        return eps * chi_at(idx[0], tv) * chi_at(idx[1], tv) * chi_at(idx[2], tv)

    nz = [l for l in partitions_le(12) if not vanishes(l)]
    agree = sum(1 for l in nz
                if all(closed_shifted(l, tv) == schur(l, ALPH(tv)) for tv in TS[:5]))
    check("section 5 control: the shifted index agrees in 0 of 1612", agree == 0 and len(nz) == 1612,
          "%d of %d" % (agree, len(nz)))

    print("\n" + ("ALL CHECKS PASS" if ok else "*** SOME CHECK FAILED ***"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
