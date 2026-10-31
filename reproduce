#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
reproduce.py -- reproduces every table and every numerical claim in

    W. Zhou, "A 2-adic logarithm law and analytic isometries for the
    Collatz odd cores of Mersenne numbers".

Only the Python standard library is used (Python >= 3.8), and all
computations are exact integer computations, except a few explicitly
marked floating-point approximations of constants (h2(alpha), rho, ...).

Usage:
    python reproduce.py                 # CSV files are written to ./output
    python reproduce.py --outdir DIR

Every computed value is compared with the value printed in the paper.
Exit status: 0 if all checks pass, 1 otherwise.
"""

import argparse
import csv
import math
import os
import sys
import time
import traceback
from collections import Counter
from fractions import Fraction
from functools import lru_cache

OUTDIR = "output"
RESULTS = []          # (name, passed, detail)

# ----------------------------------------------------------------------
# Reporting helpers
# ----------------------------------------------------------------------
def check(name, ok, detail=""):
    ok = bool(ok)
    RESULTS.append((name, ok, str(detail)))
    line = "  [{}] {}".format("PASS" if ok else "FAIL", name)
    if not ok and detail:
        line += "\n         " + str(detail)
    print(line)
    return ok

def write_csv(fname, header, rows):
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, fname)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print("  -> wrote " + path)

# ----------------------------------------------------------------------
# Basic arithmetic
# ----------------------------------------------------------------------
def v2(x):
    """2-adic valuation of a nonzero integer."""
    if x == 0:
        raise ValueError("v2(0) is undefined")
    x = abs(x)
    return (x & -x).bit_length() - 1

def kappa(n):
    return v2(3 ** n - 1)

@lru_cache(maxsize=None)
def core(n):
    """Odd core A_n = (3^n - 1) / 2^{kappa_n}."""
    x = 3 ** n - 1
    return x >> v2(x)

def collatz_step(x):
    return x >> 1 if x % 2 == 0 else 3 * x + 1

def collatz_tau(x):
    """Total stopping time: least t >= 0 with T^t(x) = 1."""
    t = 0
    while x != 1:
        x = collatz_step(x)
        t += 1
    return t

def syracuse(x):
    y = 3 * x + 1
    return y >> v2(y)

def terras_iter(x, k):
    """Return (T1^k(x), number of odd steps among the first k steps)."""
    a = 0
    for _ in range(k):
        if x & 1:
            x = (3 * x + 1) >> 1
            a += 1
        else:
            x >>= 1
    return x, a

def parity_vector(x, L):
    """First L Terras parity bits of x (x may be any residue representative)."""
    v = []
    for _ in range(L):
        b = x & 1
        v.append(b)
        x = (3 * x + 1) >> 1 if b else x >> 1
    return tuple(v)

def to_ternary(x, width=None):
    """Ternary string, most significant digit first; optional zero padding."""
    if x < 0:
        raise ValueError("negative input")
    digits = []
    chunk = 3 ** 20
    while x:
        x, r = divmod(x, chunk)
        for _ in range(20):
            r, d = divmod(r, 3)
            digits.append(d)
    s = "".join(map(str, reversed(digits))).lstrip("0") or "0"
    if width is not None:
        if x == 0 and s == "0":
            return "0" * width
        if len(s) > width:
            raise ValueError("number does not fit into {} ternary digits".format(width))
        s = s.zfill(width)
    return s

def ternary_length(x):
    return len(to_ternary(x)) if x > 0 else 0

def floor_k_log3_2(k):
    """Largest j with 3^j < 2^k, i.e. floor(k*log_3 2), by integer comparison."""
    j = 0
    while 3 ** (j + 1) < 2 ** k:
        j += 1
    return j

@lru_cache(maxsize=None)
def B(q):
    """B_q = (3^{2^q} - 1) / 2^{q+2}."""
    x = 3 ** (1 << q) - 1
    assert v2(x) == q + 2
    return x >> (q + 2)

@lru_cache(maxsize=None)
def F(j):
    """F_j = (3^{2^j} + 1) / 2."""
    return (3 ** (1 << j) + 1) // 2

def B_mod(q, L):
    """B_q mod 2^L without computing 3^{2^q} in full."""
    mod = 1 << (L + q + 2)
    t = pow(3, 1 << q, mod)
    return ((t - 1) >> (q + 2)) % (1 << L)

def Phi_mod(q, r, P):
    """Phi_q(r) mod 2^P for an integer r >= 0, using Phi_q(r) = (3^{2^q r} - 1)/2^{q+2}."""
    mod = 1 << (P + q + 2)
    t = pow(3, (1 << q) * r, mod)
    return ((t - 1) >> (q + 2)) & ((1 << P) - 1)

def lambda_mod(K=120):
    """lambda = Log(9)/8 modulo 2^(K-3), from the series Log(1+8) = sum (-1)^{j+1} 8^j / j."""
    M = 1 << K
    S = 0
    for j in range(1, K + 1):          # terms with j > K have valuation >= K
        e = v2(j)
        o = j >> e
        term = (pow(2, 3 * j - e, M) * pow(o, -1, M)) % M
        S = (S + term) % M if j % 2 == 1 else (S - term) % M
    assert S % 8 == 0 and (S >> 3) & 1 == 1      # v2(Log 9) = 3
    return S >> 3

@lru_cache(maxsize=None)
def lam_value():
    return lambda_mod(120)     # lambda mod 2^117

def h2(t):
    return -t * math.log2(t) - (1 - t) * math.log2(1 - t)

def real_root_rho():
    """Real root of x^3 = x^2 + 1 by bisection (float)."""
    lo, hi = 1.0, 2.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid ** 3 - mid ** 2 - 1 < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

ALPHA = math.log(2) / math.log(3)

# ----------------------------------------------------------------------
# Values printed in the paper
# ----------------------------------------------------------------------
EXP_A1 = [(2, 3, 7, 0, 1), (3, 1, 16, 9, 13), (4, 4, 17, 5, 5), (5, 1, 106, 95, 121),
          (6, 3, 107, 92, 91), (7, 1, 46, 31, 1093), (8, 5, 47, 26, 205),
          (9, 1, 61, 42, 9841), (10, 3, 62, 39, 7381), (11, 1, 156, 133, 88573),
          (12, 4, 157, 129, 33215), (13, 1, 158, 131, 797161), (14, 3, 159, 128, 597871),
          (15, 1, 129, 98, 7174453), (16, 6, 130, 92, 672605), (17, 1, 224, 189, 64570081),
          (18, 3, 225, 186, 48427561), (19, 1, 177, 138, 581130733),
          (20, 4, 178, 134, 217924025), (21, 1, 303, 260, 5230176601)]

EXP_T2 = [(2, 1), (3, 5), (4, 13), (5, 29), (6, 61), (7, 61), (8, 61), (9, 61),
          (10, 573), (11, 573), (12, 2621), (13, 6717), (14, 14909), (15, 14909), (16, 47677)]

EXP_LAMBDA_BITS = "101111000101110100100101000100110001110101111011"

EXP_T3 = [(4, 3, 3, "0.3750", "2.1", 14), (6, 6, 8, "0.2500", "4.6", 52),
          (8, 13, 19, "0.1484", "9.9", 194), (10, 28, 64, "0.1250", "21.3", 724),
          (12, 60, 226, "0.1104", "45.7", 2701), (14, 129, 734, "0.0896", "98.2", 10081),
          (16, 277, 2114, "0.0645", "210.9", 37622), (18, 595, 7495, "0.0572", "453.0", 140402),
          (20, 1278, 27328, "0.0521", "973.0", 523965)]

EXP_T4 = [(6, 1, 0, False, "27.400"), (8, 3, 5, True, "0.846"), (10, 7, 5, True, "0.844"),
          (12, 1, 2, False, "6.750"), (16, 3, 6, True, "0.422"), (20, 2, 3, False, "3.375"),
          (24, 1, 3, False, "3.375"), (32, 3, 7, True, "0.211"), (40, 2, 4, False, "1.688"),
          (64, 3, 8, True, "0.105"), (266, 9, 7, True, "0.211")]

EXP_T5 = [(3, "01", 1, 1, 1, 1, 0), (4, "0012", 5, 2, 2, 1, 1),
          (5, "00021121", 205, 3, 3, 3, 2), (6, "0001021011122022", 672605, 3, 6, 5, 5)]

EXP_T6 = ["0012", "0120", "0221", "1022", "1200", "2001", "2102", "2210"]

EXP_A2 = [(1, 1, 1, 1, 0, 0, 1), (2, 2, 2, 1, 1, 0, 3), (3, 5, 3, 3, 2, 0, 7),
          (4, 13, 3, 5, 5, 3, 15), (5, 28, 4, 11, 10, 7, 31), (6, 59, 5, 21, 21, 17, 63),
          (7, 123, 5, 43, 42, 38, 127), (8, 250, 6, 85, 85, 80, 255),
          (9, 506, 6, 171, 170, 165, 511), (10, 1017, 7, 341, 341, 335, 1023),
          (11, 2040, 8, 683, 682, 675, 2047), (12, 4088, 8, 1365, 1365, 1358, 4095)]

EXP_A3 = {1: [5], 2: [41], 3: [17, 193], 4: [21523361], 5: [926510094425921]}

EXP_A4 = [(3, [2], 2), (5, [2], 2), (17, [2], 3), (41, [2, 5], 6), (193, [2, 3], 5),
          (21523361, [2, 5, 17, 41, 193], 6),
          (926510094425921, [2, 5, 17, 41, 193, 21523361], 7)]

EXP_A5 = [(6, 4, 32), (8, 3, 5), (12, 7, 36), (24, 6, 117), (64, 5, 211), (100, 4, 427),
          (150, 7, 563), (200, 3, 777), (256, 4, 818), (300, 3, 1472)]

# ----------------------------------------------------------------------
# Section 2: Mersenne prefix, core reduction, Table A1
# ----------------------------------------------------------------------
def section_prefix():
    check("Lemma 2.1: v2(3^n - 1) = 1 (n odd), q+2 (n = 2^q r), 1 <= n <= 300",
          all(kappa(n) == (1 if n % 2 else v2(n) + 2) for n in range(1, 301)))

    ok22 = ok24 = True
    for n in range(1, 41):
        orbit = [2 ** n - 1]
        for _ in range(2 * n + kappa(n)):
            orbit.append(collatz_step(orbit[-1]))
        for t in range(2 * n + 1):
            j, odd = divmod(t, 2)
            exp = (3 ** j * 2 ** (n - j) - 1) if odd == 0 else 2 * (3 ** (j + 1) * 2 ** (n - j - 1) - 1)
            if orbit[t] != exp:
                ok22 = False
        if orbit[2 * n + kappa(n)] != core(n):
            ok24 = False
    check("Theorem 2.2: Mersenne prefix formulas, 1 <= n <= 40", ok22)
    check("Theorem 2.4: T^{2n+kappa_n}(M_n) = A_n, 1 <= n <= 40", ok24)

    check("Lemma 2.3: A_n = 1 iff n in {1,2}; A_n >= 5 for 3 <= n <= 200",
          core(1) == 1 and core(2) == 1 and all(core(n) >= 5 for n in range(3, 201)))
    check("Lemma 2.3: 4n*A_n >= 3^n - 1 for even 2 <= n <= 200",
          all(4 * n * core(n) >= 3 ** n - 1 for n in range(2, 201, 2)))
    check("Proposition 2.6: 3 does not divide A_n, 1 <= n <= 300",
          all(core(n) % 3 != 0 for n in range(1, 301)))
    check("Proposition 2.7: C(A_n) = A_{n+1} < A_n for odd 3 <= n <= 299",
          all(syracuse(core(n)) == core(n + 1) < core(n) for n in range(3, 300, 2)))

    tauM = {n: collatz_tau(2 ** n - 1) for n in range(1, 101)}
    tauA = {n: collatz_tau(core(n)) for n in range(1, 101)}
    check("Theorem 2.4: tau(M_n) = 2n + kappa_n + tau(A_n), 2 <= n <= 100",
          all(tauM[n] == 2 * n + kappa(n) + tauA[n] for n in range(2, 101)))
    check("Theorem 2.5: tau(M_{n+1}) = tau(M_n) + 1 for odd 3 <= n <= 99",
          all(tauM[n + 1] == tauM[n] + 1 for n in range(3, 100, 2)))
    ok25 = True
    for n in range(3, 40, 2):
        k = kappa(n + 1)
        x = 2 ** n - 1
        for _ in range(2 * n + k + 1):
            x = collatz_step(x)
        if x != core(n + 1):
            ok25 = False
    check("Theorem 2.5: T^{2n+k+1}(M_n) = A_{n+1}, odd 3 <= n <= 39", ok25)

    rows = [(n, kappa(n), tauM[n], tauA[n], core(n)) for n in range(2, 101)]
    write_csv("tableA1_tau_Mn.csv", ["n", "kappa_n", "tau(M_n)", "tau(A_n)", "A_n"], rows)
    check("Table A1 (n <= 21) matches the paper", rows[:20] == EXP_A1)

# ----------------------------------------------------------------------
# Section 3: lambda and the 2-adic functions Phi_q
# ----------------------------------------------------------------------
def section_2adic():
    lam = lam_value()

    rows = [(L, B_mod(L - 1, L), format(B_mod(L - 1, L), "0{}b".format(L))) for L in range(2, 17)]
    write_csv("table2_lambda_mod_2L.csv", ["L", "b_L", "binary (most significant first)"], rows)
    check("Table 2 matches the paper", [(r[0], r[1]) for r in rows] == EXP_T2)
    check("Corollary 3.5: B_{L-1} mod 2^L = lambda mod 2^L for 2 <= L <= 60",
          all(B_mod(L - 1, L) == lam % (1 << L) for L in range(2, 61)))
    check("lambda is a 2-adic unit and lambda = 5 (mod 8)", lam % 8 == 5)
    check("lambda = 573 (mod 1024)", lam % 1024 == 573)
    bits = "".join(str((lam >> i) & 1) for i in range(48))
    check("first 48 binary digits of lambda (least significant first)",
          bits == EXP_LAMBDA_BITS, bits)
    write_csv("lambda_binary_digits.csv", ["bit index i (2^i)", "digit"],
              [(i, (lam >> i) & 1) for i in range(117)])

    # Theorem 3.3(ii): exact interpolation on odd integers
    check("Theorem 3.3(ii): (3^{2^q r}-1)/2^{q+2} = A(2^q r), q <= 4, odd r <= 31",
          all(((3 ** ((1 << q) * r) - 1) >> (q + 2)) == core((1 << q) * r)
              for q in range(1, 5) for r in range(1, 32, 2)))

    # Theorem 3.3(iii): isometry on integers
    P = 24
    iso = True
    for q in range(1, 7):
        vals = [Phi_mod(q, r, P) for r in range(129)]
        for r in range(129):
            for s in range(r):
                d = (vals[r] - vals[s]) % (1 << P)
                if d == 0 or v2(d) != v2(r - s):
                    iso = False
    check("Theorem 3.3(iii): v2(Phi_q(x)-Phi_q(x')) = v2(x-x'), q <= 6, 0 <= x,x' <= 128", iso)

    # Theorem 3.3(iv): permutation of Z/2^L preserving odd residues
    perm = True
    for q in range(1, 6):
        for L in range(1, 13):
            M = 1 << L
            img = [Phi_mod(q, r, L) for r in range(M)]
            if len(set(img)) != M or any((img[r] - r) & 1 for r in range(M)):
                perm = False
    check("Theorem 3.3(iv): Phi_q permutes Z/2^L and preserves parity, q <= 5, L <= 12", perm)

    # Theorem 3.4: two congruences
    c1 = c2 = True
    for q in range(1, 13):
        Bq = B_mod(q, q + 2)
        for r in range(1, 301):
            val = Phi_mod(q, r, q + 2)
            if (val - lam * r) % (1 << (q + 1)):
                c1 = False
            if (val - Bq * r) % (1 << (q + 2)):
                c2 = False
    check("Theorem 3.4: Phi_q(x) = lambda*x (mod 2^{q+1}), q <= 12, 1 <= x <= 300", c1)
    check("Theorem 3.4: Phi_q(x) = B_q*x (mod 2^{q+2}), q <= 12, 1 <= x <= 300", c2)

    # Corollary 3.6
    ok_a = ok_b = True
    for q in range(1, 11):
        for r in range(1, 1024, 2):
            a8 = Phi_mod(q, r, 3)
            if a8 % 4 != r % 4:
                ok_a = False
            if a8 != ((r if q == 1 else 5 * r) % 8):
                ok_b = False
    check("Corollary 3.6(a),(b): A(2^q r) mod 4 and mod 8, q <= 10, odd r < 1024", ok_a and ok_b)
    check("Corollary 3.6(a): C(A_n) < A_n iff r = 1 (mod 4), even 4 <= n <= 400",
          all((syracuse(core(n)) < core(n)) == ((n >> v2(n)) % 4 == 1) for n in range(4, 401, 2)))

    # Proposition 3.7
    ok37 = True
    for q in range(1, 11):
        M = 1 << (q + 2)
        cnt = Counter()
        sols = {}
        for r in range(1, 1 << (q + 1), 2):
            x = (3 * Phi_mod(q, r, q + 2) + 1) % M
            e = v2(x) if x else q + 2
            cnt[e] += 1
            sols.setdefault(e, set()).add(r)
        for e in range(1, q + 1):
            if cnt[e] != 1 << (q - e):
                ok37 = False
            mod = 1 << (e + 1)
            target = ((2 ** e - 1) * pow(3 * lam, -1, mod)) % mod
            if any(r % mod != target for r in sols.get(e, ())):
                ok37 = False
    check("Proposition 3.7: v2(3A+1) = e on exactly one class mod 2^{e+1}, proportion 2^{-e}", ok37)

    # Corollary 3.9(iii)
    check("Corollary 3.9(iii): first q+1 parity bits of B_q equal those of lambda, q <= 60",
          all(parity_vector(B_mod(q, q + 1), q + 1) == parity_vector(lam % (1 << (q + 1)), q + 1)
              for q in range(1, 61)))

# ----------------------------------------------------------------------
# Section 4: Terras bijection, exceptional sets, Table 3
# ----------------------------------------------------------------------
def is_exceptional(x, L):
    a = 0
    for k in range(1, L + 1):
        b = x & 1
        a += b
        x = (3 * x + 1) >> 1 if b else x >> 1
        if 3 ** a < 2 ** k:
            return False
    return True

def exceptional_count(L):
    """|E_L^T| via dynamic programming over the number of odd steps."""
    states = {0: 1}
    for k in range(1, L + 1):
        new = {}
        for a, c in states.items():
            for b in (0, 1):
                a2 = a + b
                if 3 ** a2 >= 2 ** k:
                    new[a2] = new.get(a2, 0) + c
        states = new
    return sum(states.values())

def section_terras():
    lam = lam_value()
    check("Lemma 2.9: parity vectors give a bijection Z/2^L -> {0,1}^L, L <= 14",
          all(len({parity_vector(x, L) for x in range(1 << L)}) == 1 << L for L in range(1, 15)))

    rho = real_root_rho()
    h2a = h2(ALPHA)
    g = [1, 1, 1]
    for L in range(3, 41):
        g.append(g[L - 1] + g[L - 3])

    rows = []
    for L in range(4, 21, 2):
        E = exceptional_count(L)
        rows.append((L, g[L], E, "{:.4f}".format(E / 2 ** (L - 1)),
                     "{:.1f}".format(rho ** (L - 2)), round(2 ** (h2a * L))))
    write_csv("table3_exceptional_sets.csv",
              ["L", "g(L)", "|E_L^T|", "delta_L", "rho^(L-2)", "2^(h2(alpha)L)"], rows)
    check("Table 3 matches the paper", rows == EXP_T3, rows)

    check("|E_L^T| by DP equals brute force over residues, L <= 14",
          all(exceptional_count(L) == sum(1 for x in range(1 << L) if is_exceptional(x, L))
              for L in range(1, 15)))
    check("Theorem 4.4(ii): -1 is exceptional, L <= 40",
          all(is_exceptional((1 << L) - 1, L) for L in range(1, 41)))
    ok_b = True
    for L in range(1, 41):
        E = exceptional_count(L)
        if not (rho ** (L - 2) <= g[L] + 1e-9 and g[L] <= E <= 2 ** (h2a * L)
                and E / 2 ** (L - 1) <= 2 ** (1 - 0.05 * L)):
            ok_b = False
    check("Theorem 4.4(iii): rho^{L-2} <= g(L) <= |E| <= 2^{h2(alpha)L}, delta_L <= 2^{1-0.05L}, L <= 40", ok_b)

    ok_i = True
    for L in range(2, 11):
        q = L - 1
        M = 1 << L
        E_set = {x for x in range(1, M, 2) if is_exceptional(x, L)}
        lhs = {r for r in range(1, M, 2) if Phi_mod(q, r, L) in E_set}
        rhs = {r for r in range(1, M, 2) if (lam * r) % M in E_set}
        if lhs != rhs or len(lhs) != len(E_set):
            ok_i = False
    check("Theorem 4.4(i): E_{q,L} = lambda^{-1} E_L^T for q = L-1, 2 <= L <= 10", ok_i)

    ok45 = True
    for q in range(1, 5):
        for L in range(1, 13):
            vecs = {parity_vector(Phi_mod(q, r, L), L) for r in range(1, 1 << L, 2)}
            if len(vecs) != 1 << (L - 1) or any(v[0] != 1 for v in vecs):
                ok45 = False
    check("Theorem 4.5(i): parity vectors of A(2^q r) run once over {1}x{0,1}^{L-1}, q <= 4, L <= 12", ok45)

    # constants and integer inequalities
    check("3^6309 < 2^10000 (alpha > 0.6309)", 3 ** 6309 < 2 ** 10000)
    check("2^9500 * 6309^6309 * 3691^3691 > 10^40000 (h2(0.6309) < 0.95)",
          2 ** 9500 * 6309 ** 6309 * 3691 ** 3691 > 10 ** 40000)
    check("h2(alpha) = 0.94996... and 1 - h2(alpha) = 0.05004... (float)",
          abs(h2a - 0.94996) < 2e-5 and abs((1 - h2a) - 0.05004) < 2e-5, h2a)
    check("rho = 1.46557..., log2(rho) = 0.5515... (float)",
          abs(rho - 1.46557) < 1e-5 and abs(math.log2(rho) - 0.5515) < 1e-4, rho)
    check("Theorem 4.3: 3^{2L+4} - 1 >= 4(2L+4) 6^L for 1 <= L <= 300",
          all(3 ** (2 * L + 4) - 1 >= 4 * (2 * L + 4) * 6 ** L for L in range(1, 301)))

# ----------------------------------------------------------------------
# Section 5: crossing, Table 4
# ----------------------------------------------------------------------
def crossing_data(n):
    a = core(n)
    s = v2(3 * a + 1)
    E = s + v2(n) - v2(n - 2)
    c = syracuse(a)
    prev = core(n - 2)
    return s, E, c < prev, Fraction(c, prev)

def crossing_predicted(n):
    q = v2(n)
    r = n >> q
    return (q >= 5 or (q == 4 and r % 4 == 1) or (q == 3 and r % 8 == 1)
            or (q == 2 and r % 16 == 1) or (q == 1 and r % 128 == 5))

def section_crossing():
    sel = [6, 8, 10, 12, 16, 20, 24, 32, 40, 64, 266]
    rows, comp = [], []
    for n in sel:
        s, E, cr, ratio = crossing_data(n)
        q = v2(n)
        rstr = "{:.3f}".format(float(ratio))
        rows.append((n, "{}*{}".format(1 << q, n >> q), s, E, "yes" if cr else "no", rstr))
        comp.append((n, s, E, cr, rstr))
    write_csv("table4_one_step_crossing.csv", ["n", "2^q*r", "s", "E", "crosses", "C(A_n)/A_{n-2}"], rows)
    check("Table 4 matches the paper", comp == EXP_T4, comp)

    crossing = [n for n in range(6, 801, 2) if crossing_data(n)[2]]
    write_csv("crossing_n_6_to_800.csv", ["n with C(A_n) < A_{n-2}"], [(n,) for n in crossing])
    check("Theorem 5.1: 67 crossing n among even 6 <= n <= 800", len(crossing) == 67, len(crossing))
    check("Theorem 5.1: crossing set = five congruence conditions, 6 <= n <= 800",
          crossing == [n for n in range(6, 801, 2) if crossing_predicted(n)])
    check("Theorem 5.1: crossing iff E >= 5, even 6 <= n <= 800",
          all((crossing_data(n)[1] >= 5) == crossing_data(n)[2] for n in range(6, 801, 2)))

    dens = (Fraction(1, 16) + Fraction(1, 16) * Fraction(1, 2) + Fraction(1, 8) * Fraction(1, 4)
            + Fraction(1, 4) * Fraction(1, 8) + Fraction(1, 2) * Fraction(1, 64))
    check("Theorem 5.1: density = 21/128 (exact sum)", dens == Fraction(21, 128))
    check("Theorem 5.1: exactly 21 of the 128 even n in [256, 510] satisfy the congruences",
          sum(crossing_predicted(n) for n in range(256, 512, 2)) == 21)

    s4, E4, _, _ = crossing_data(4)
    check("Remark 5.2: n = 4 gives s = 4, E = 5 but C(A_4) = A_2 = 1, and 27*9-3+16 = 256 = 2^5*8",
          s4 == 4 and E4 == 5 and syracuse(core(4)) == core(2) == 1 and 27 * 9 - 3 + 16 == 256 == 32 * 8)

    ok54 = True
    for L in range(1, 301):
        n = 3 * L + 10
        lhs, rhs = 3 ** (n - 2), n * 2 ** (L + 2) * 3 ** (L + 1)
        if lhs < rhs or Fraction(lhs, rhs) != Fraction(9, 2) ** L * 3 ** 7 / (4 * (3 * L + 10)):
            ok54 = False
    check("Theorem 5.4: 3^{n-2} >= n 2^{L+2} 3^{L+1} at n = 3L+10, and the ratio formula, L <= 300", ok54)
    check("Corollary 5.5: 2^1000 < 3^631 (alpha < 0.631)", 2 ** 1000 < 3 ** 631)
    check("Corollary 5.5: 631^40 < 2^31 * 369^40", 631 ** 40 < 2 ** 31 * 369 ** 40)
    check("Corollary 5.5: log2(alpha/(1-alpha)) = 0.7736... (float)",
          abs(math.log2(ALPHA / (1 - ALPHA)) - 0.7736) < 1e-4)

# ----------------------------------------------------------------------
# Section 6: ternary side, Choi's strings, Tables 5, 6, A2, A3, A4
# ----------------------------------------------------------------------
def factor_by_step(n, step):
    """Factor n by trial division over candidates 1 + t*step (valid when every prime factor is 1 mod step)."""
    fac = []
    d = step + 1
    while d * d <= n:
        while n % d == 0:
            fac.append(d)
            n //= d
        d += step
    if n > 1:
        fac.append(n)
    return fac

def least_witness(p, primes):
    g = 2
    while True:
        if pow(g, p - 1, p) == 1 and all(pow(g, (p - 1) // t, p) != 1 for t in primes):
            return g
        g += 1

def fmt_factorization(n, primes):
    parts = []
    for t in primes:
        e = 0
        while n % t == 0:
            n //= t
            e += 1
        parts.append("{}^{}".format(t, e) if e > 1 else str(t))
    return "*".join(parts) if n == 1 else None

def section_ternary():
    # Theorem 6.2 and Example 6.3
    check("Theorem 6.2(1): A_n = (11...1)_3 for odd n <= 101",
          all(to_ternary(core(n)) == "1" * n for n in range(1, 102, 2)))
    check("Theorem 6.2(2): block structure of A(2^q r), q <= 5, odd r <= 15",
          all(to_ternary(core((1 << q) * r), r * (1 << q)) == to_ternary(B(q), 1 << q) * r
              for q in range(1, 6) for r in range(1, 16, 2)))
    check("Example 6.3: A_24 = 8825923015 = (21121 000 21121 000 21121)_3",
          core(24) == 8825923015 and to_ternary(core(24)) == "211210002112100021121")

    # Table A2, Theorems 6.7, 6.15, Corollary 6.16
    rows = []
    ok615 = ok616 = True
    for q in range(1, 13):
        b = B(q)
        d = 1 << q
        std = to_ternary(b)
        fix = to_ternary(b, d)
        ell = len(std)
        J = (2 ** q - (-1) ** q) // 3
        n1, n2, n0fix = fix.count("1"), fix.count("2"), fix.count("0")
        s3 = sum(int(c) for c in std)
        rows.append((q, ell, d - ell, n1, n2, std.count("0"), s3))
        if (n1, n2, n0fix) != (J, (d - 1 - J) // 2, (d + 1 - J) // 2):
            ok615 = False
        if (s3 != d - 1 or (n1 + n2) != (2 ** (q + 1)) // 3
                or std.count("0") != (d + 1 - J) // 2 - floor_k_log3_2(q + 2)):
            ok616 = False
    write_csv("tableA2_block_parameters.csv",
              ["q", "ell(B_q)", "G_q", "#1", "#2", "#0 (standard rep.)", "s3(B_q)"], rows)
    check("Table A2 matches the paper", rows == EXP_A2, rows)
    check("Theorem 6.15: digit counts of B_q are Jacobsthal-type, q <= 12", ok615)
    check("Corollary 6.16(i),(ii): s3(B_q), nonzero digits and zeros, q <= 12", ok616)
    check("Corollary 6.16(iii): digit sum of A(2^q r) = r(2^q - 1), q <= 5, odd r <= 15",
          all(sum(int(c) for c in to_ternary(core((1 << q) * r))) == r * (2 ** q - 1)
              for q in range(1, 6) for r in range(1, 16, 2)))
    check("Theorem 6.7: G_q = floor((q+2) log_3 2), q <= 14",
          all((1 << q) - ternary_length(B(q)) == floor_k_log3_2(q + 2) for q in range(1, 15)))
    check("Remark 6.8: G_q - q*alpha in [2alpha-1, 2alpha), q <= 200",
          all(2 * ALPHA - 1 <= floor_k_log3_2(q + 2) - q * ALPHA < 2 * ALPHA for q in range(1, 201)))

    # Theorem 6.5, Corollary 6.4, Tables A3, A4
    check("Corollary 6.4 / Theorem 6.5: B_q = prod_{j<q} F_j, q <= 12",
          all(B(q) == math.prod(F(j) for j in range(1, q)) for q in range(1, 13)))
    check("Theorem 6.5: F_i, F_j pairwise coprime, i < j <= 10",
          all(math.gcd(F(i), F(j)) == 1 for j in range(2, 11) for i in range(1, j)))
    check("F_j - 1 = 2^{j+1} B_j, j <= 12",
          all(F(j) - 1 == 2 ** (j + 1) * B(j) for j in range(1, 13)))

    rowsA3 = []
    okA3 = True
    for j in range(1, 6):
        fac = factor_by_step(F(j), 1 << (j + 1))
        okA3 &= (fac == EXP_A3[j] and math.prod(fac) == F(j)
                 and all(p % (1 << (j + 1)) == 1 for p in fac)
                 and all(least_witness(p, [t for t in set(factor_by_step(p - 1, 1)) ]) for p in fac))
        rowsA3.append((j, F(j), " x ".join(map(str, fac)), 1 << (j + 1),
                       "yes" if all(p % (1 << (j + 1)) == 1 for p in fac) else "no"))
    write_csv("tableA3_Fj_factorization.csv",
              ["j", "F_j", "factorization", "2^(j+1)", "all prime factors = 1 mod 2^(j+1)"], rowsA3)
    check("Table A3: factorizations of F_j, j <= 5 (trial division over p = 1 mod 2^{j+1})", okA3)

    rowsA4 = []
    okA4 = True
    certified = set()
    for p, primes, g_exp in EXP_A4:
        fstr = fmt_factorization(p - 1, primes)
        g = least_witness(p, primes)
        okA4 &= (fstr is not None and g == g_exp and all(t == 2 or t in certified for t in primes))
        certified.add(p)
        rowsA4.append((p, fstr, g))
    write_csv("tableA4_pratt_certificates.csv", ["prime p", "factorization of p-1", "least witness g"], rowsA4)
    check("Table A4: Pratt certificates (factorization of p-1, least witness, recursive)", okA4, rowsA4)
    check("F_4 and F_5 are the primes certified in Table A4",
          EXP_A4[5][0] == F(4) and EXP_A4[6][0] == F(5))

    # Lemma 6.6
    check("Lemma 6.6: 3^{2^{k-2}} - 1 = 2^k B_{k-2}, 3 <= k <= 16",
          all(3 ** (1 << (k - 2)) - 1 == 2 ** k * B(k - 2) for k in range(3, 17)))
    check("Lemma 6.6: order of 3 mod 2^k is 2^{k-2}, 3 <= k <= 64",
          all(pow(3, 1 << (k - 2), 1 << k) == 1 and pow(3, 1 << (k - 3), 1 << k) != 1
              for k in range(3, 65)))

    # Table 5
    rows5 = []
    for k in range(3, 7):
        s = to_ternary(B(k - 2), 1 << (k - 2))
        lz = len(s) - len(s.lstrip("0"))
        rows5.append((k, s, B(k - 2), lz, s.count("0"), s.count("1"), s.count("2")))
    write_csv("table5_choi_rk.csv", ["k", "r_k", "B_{k-2}", "leading zeros", "#0", "#1", "#2"], rows5)
    check("Table 5 matches the paper", rows5 == EXP_T5, rows5)
    check("Corollary 6.12: r_k has floor(k log_3 2) leading zeros, 3 <= k <= 16",
          all((lambda s: len(s) - len(s.lstrip("0")))(to_ternary(B(k - 2), 1 << (k - 2)))
              == floor_k_log3_2(k) for k in range(3, 17)))

    # Table 6 and Corollary 6.13
    rows6 = []
    for v in range(1, 16, 2):
        bv = to_ternary(5 * v, 4)
        rows6.append((v, v % 8, 5 * v, bv,
                       "cyclic shift of beta_1" if v % 8 in (1, 3) else "cyclic shift of complement",
                       "{}/{}/{}".format(bv.count("0"), bv.count("1"), bv.count("2"))))
    write_csv("table6_beta_v_k4.csv", ["v", "v mod 8", "v*B_2", "beta_v", "type", "#0/#1/#2"], rows6)
    check("Table 6 matches the paper", [r[3] for r in rows6] == EXP_T6)
    ok613 = True
    swap = str.maketrans("012", "210")
    for k in range(3, 11):
        m = 1 << (k - 2)
        b1 = to_ternary(B(k - 2), m)
        comp = b1.translate(swap)
        for v in range(1, 1 << k, 2):
            bv = to_ternary(v * B(k - 2), m)
            target = b1 if v % 8 in (1, 3) else comp
            if bv not in target + target:
                ok613 = False
    check("Corollary 6.13: beta_v is a cyclic shift of beta_1 or of its complement by v mod 8, k <= 10", ok613)

    # Theorem 6.11: exact decomposition, checked string by string
    ok611 = True
    bad = ""
    for k in range(3, 8):
        m = 1 << (k - 2)
        Bk = B(k - 2)
        for y in range(1, 41):
            iota = v2(y)
            yp = y >> iota
            I, v = divmod(yp, 1 << k)
            beta = to_ternary(v * Bk, m)
            ylen = len(to_ternary(y))
            for N in range(0, 3 * m + 8):
                X = yp * 3 ** N
                Q, w = divmod(X, 1 << k)
                z, a = terras_iter(X, k)
                zz, _ = terras_iter(y * 3 ** N, k + iota)
                tw, aw = terras_iter(w, k)
                cond = (zz == z and aw == a and a >= 1 and z == Q * 3 ** a + tw and tw < 3 ** a
                        and Q == I * 3 ** N + (v * 3 ** N) // (1 << k))
                if cond:
                    try:
                        got = to_ternary(z, ylen + N + a)
                    except ValueError:
                        got = None
                    exp = (to_ternary(I, ylen) + beta * (N // m) + beta[:N % m] + to_ternary(tw, a))
                    cond = (got == exp)
                if not cond:
                    ok611 = False
                    bad = "k={}, y={}, N={}".format(k, y, N)
    check("Theorem 6.11: decomposition (6.2)-(6.3), 3 <= k <= 7, 1 <= y <= 40, 0 <= N < 3m+8", ok611, bad)

    okmin = True
    for k in range(3, 9):
        m = 1 << (k - 2)
        b1 = to_ternary(B(k - 2), m)
        for N in range(3 * m, 3 * m + 12):
            z, a = terras_iter(3 ** N, k)
            s = to_ternary(z, 1 + N + a)
            pmin = next(p for p in range(1, N + 1)
                        if s[1:1 + p * (N // p)] == s[1:1 + p] * (N // p))
            if pmin != m or s[1:1 + m] != b1 or s[0] != "0":
                okmin = False
    check("Theorem 6.11(i): shortest repeating string of g^k(10^N) is beta_1 of length 2^{k-2}, k <= 8, N >= 3m", okmin)

# ----------------------------------------------------------------------
# Appendix A: first meetings with earlier cores, Table A5
# ----------------------------------------------------------------------
def section_meetings():
    nmax = 300
    cores = [None] + [core(m) for m in range(1, nmax + 1)]
    res = {}
    for n in range(4, nmax + 1, 2):
        targets = {}
        for m in range(n - 2, 0, -1):          # smaller m overwrites: A_1 = A_2 = 1 -> m = 1
            targets[cores[m]] = m
        x, steps = cores[n], 0
        while x not in targets:
            x = syracuse(x)
            steps += 1
        res[n] = (targets[x], steps)
    write_csv("meetings_all_even_n_4_to_300.csv", ["n", "first core met (index m)", "Syracuse steps"],
              [(n, m, s) for n, (m, s) in sorted(res.items())])
    rows = [(n, res[n][0], res[n][1]) for n, _, _ in EXP_A5]
    write_csv("tableA5_first_meetings.csv", ["n", "first core met (index m)", "steps"], rows)
    check("Table A5 matches the paper", rows == EXP_A5, rows)
    check("all 149 even n in [4, 300] meet an earlier core A_m, m <= n-2",
          len(res) == 149 and all(m <= n - 2 for n, (m, _) in res.items()))
    nmaxstep = max(res, key=lambda n: res[n][1])
    check("largest number of steps: 1472 at n = 300, meeting A_3",
          (nmaxstep, res[nmaxstep]) == (300, (3, 1472)), (nmaxstep, res[nmaxstep]))
    cnt = Counter(m for m, _ in res.values())
    check("most frequent first targets: A_4 (71), A_3 (27), A_7 (24), A_5 (22)",
          (cnt[4], cnt[3], cnt[7], cnt[5]) == (71, 27, 24, 22) and
          [m for m, _ in cnt.most_common(4)] == [4, 3, 7, 5], cnt.most_common(6))

# ----------------------------------------------------------------------
SECTIONS = [
    ("Section 2: Mersenne prefix and cores (Table A1)", section_prefix),
    ("Section 3: lambda and Phi_q (Table 2)", section_2adic),
    ("Section 4: parity vectors and exceptional sets (Table 3)", section_terras),
    ("Section 5: crossing A_{n-2} (Table 4)", section_crossing),
    ("Section 6: ternary side (Tables 5, 6, A2, A3, A4)", section_ternary),
    ("Appendix A: meetings with earlier cores (Table A5)", section_meetings),
]

def main():
    global OUTDIR
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", default="output", help="directory for CSV files (default: output)")
    OUTDIR = ap.parse_args().outdir

    t0 = time.time()
    print("Python", sys.version.split()[0])
    for title, fn in SECTIONS:
        print("\n== " + title + " ==")
        t1 = time.time()
        try:
            fn()
        except Exception:
            check(title + ": finished without exception", False, traceback.format_exc())
        print("  ({:.1f} s)".format(time.time() - t1))

    write_csv("check_report.csv", ["check", "passed", "detail"],
              [(n, "yes" if ok else "no", d if not ok else "") for n, ok, d in RESULTS])
    npass = sum(ok for _, ok, _ in RESULTS)
    nfail = len(RESULTS) - npass
    print("\nSummary: {} PASS, {} FAIL, total time {:.1f} s".format(npass, nfail, time.time() - t0))
    sys.exit(0 if nfail == 0 else 1)

if __name__ == "__main__":
    main()
