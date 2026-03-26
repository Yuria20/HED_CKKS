#!/usr/bin/env python3
# dp_full.py
#
# Full-grid variant derived from dp.py (banded/windowed symbol-mapped Myers):
# - Computes Cha's complex symbol-mapped depth-2 recurrence on the FULL (L+1)x(L+1) DP grid
#   in anti-diagonal order k=i+j, keeping ALL valid cells on each diagonal.
# - No data-dependent branching (no delta threshold); eq is computed by a polynomial.
# - Outputs Hsym/Vsym per diagonal (object arrays) compatible with the existing pack.py style.
# - Optionally emits integer-quantized H/V and a recovered distance (plaintext-only debug).
#
# Usage:
#   python3 dp_full.py -i out/pairsrnd.fa -p full -k 0 -o out
#
# Notes:
# - This is for correctness/benchmarking: compares against banded/windowed approximation.
# - In the full grid, "missing neighbors" should only occur at boundaries, which we handle
#   explicitly (top row / left column). There is no missing-neighbor handling in full DP (neighbors always exist for interior cells).

import argparse
from pathlib import Path
import numpy as np


# ------------------------------------------------------------
# Reference Levenshtein DP (exact) for debugging/validation
def levenshtein_dp(a: str, b: str) -> int:
    """Exact Levenshtein distance (unit costs). O(m*n) time, O(n) memory."""
    m, n = len(a), len(b)
    if m == 0:
        return n
    if n == 0:
        return m
    dp = list(range(n + 1))
    for i, ca in enumerate(a, start=1):
        prev = dp[0]
        dp[0] = i
        for j, cb in enumerate(b, start=1):
            tmp = dp[j]
            cost = 0 if ca == cb else 1
            dp[j] = min(dp[j] + 1, dp[j - 1] + 1, prev + cost)
            prev = tmp
    return dp[n]

# ------------------------------------------------------------
# 1) Pair FASTA parser (4-line blocks)
#    >id
#    READ
#    <refinfo
#    REF
def iter_pairs(path: Path):
    with open(path, "r") as f:
        while True:
            h1 = f.readline()
            if not h1:
                break
            s1 = f.readline()
            h2 = f.readline()
            s2 = f.readline()
            if not (s1 and h2 and s2):
                break

            h1 = h1.strip()
            s1 = s1.strip()
            h2 = h2.strip()
            s2 = s2.strip()

            if not h1.startswith(">") or not h2.startswith("<"):
                continue

            pid = h1[1:]
            yield pid, s1, s2


def get_pair(path: Path, k: int):
    for idx, (pid, a, b) in enumerate(iter_pairs(path)):
        if idx == k:
            return pid, a, b
    raise SystemExit(f"pair index {k} out of range")


# ------------------------------------------------------------
# 2) Equality delta(i,j) as a *polynomial* in (a-b)^2
MAP = {"A": -1.5, "C": -0.5, "G": 0.5, "T": 1.5}


def eq_poly(z):
    # 3rd-degree polynomial with:
    # e(0)=1, e(1)=0, e(4)=0, e(9)=0
    # Horner form:
    return 1.0 + z * (-49.0 / 36.0 + z * (7.0 / 18.0 + z * (-1.0 / 36.0)))


def eq_from_chars(a_ch: str, b_ch: str) -> float:
    if a_ch == "N":
        a_ch = "A"
    if b_ch == "N":
        b_ch = "A"
    if a_ch not in MAP or b_ch not in MAP:
        return 0.0
    d = MAP[a_ch] - MAP[b_ch]
    z = d * d
    val = float(eq_poly(z))
    if val < 0.0:
        val = 0.0
    elif val > 1.0:
        val = 1.0
    return val


# ------------------------------------------------------------
# 3) Full diagonal indices (anti-diagonal k=i+j)
def full_diag_indices(L, k):
    """Return the full i-range (i0,i1) for anti-diagonal k=i+j on (L+1)x(L+1) grid."""
    i_min = max(0, k - L)
    i_max = min(L, k)
    if i_min > i_max:
        return None
    return (i_min, i_max)


# ---------------------------------------------
# 4) Cha symbol map g(u) and depth-2 recurrence
#
# g(-1)=1/2, g(0)=0, g(+1)=i/2
#
def g_map(u: int) -> complex:
    if u == -1:
        return 0.5 + 0.0j
    if u == 0:
        return 0.0 + 0.0j
    if u == 1:
        return 0.0 + 0.5j
    raise ValueError("u must be in {-1,0,1}")


def g_inv(z: complex) -> int:
    # nearest neighbor among {1/2, 0, i/2}
    cands = [(-1, 0.5 + 0.0j), (0, 0.0 + 0.0j), (1, 0.0 + 0.5j)]
    best_u = None
    best_d = None
    for u, w in cands:
        d = abs(z - w)
        if best_d is None or d < best_d:
            best_d = d
            best_u = u
    return int(best_u)


def OR(a, b):
    # boolean OR as polynomial (valid for {0,1})
    # Numeric guard: dy/dz are intended to be real booleans. Round+clip prevents drift
    # and stops the polynomial from exploding when inputs leave {0,1}.
    a = np.clip(np.round(np.real(a)), 0.0, 1.0)
    b = np.clip(np.round(np.real(b)), 0.0, 1.0)
    return a + b - a * b


def H_tilde(x, y, z):
    # x in [-1,0], y,z are symbol-mapped
    dy = y + np.conj(y)
    dz = z + np.conj(z)
    o = OR(dy, dz)
    not_o = 1.0 - o
    I = 1.0j
    return (I * np.conj(y)) * o + (x * (I * np.conj(z) + np.conj(y)) + I * np.conj(z)) * not_o


def V_tilde(x, y, z):
    dy = y + np.conj(y)
    dz = z + np.conj(z)
    o = OR(dy, dz)
    not_o = 1.0 - o
    I = 1.0j
    return (I * np.conj(z)) * o + (x * (np.conj(z) + I * np.conj(y)) + np.conj(y)) * not_o


def myers_step_symbolmapped_sym(delta: float, vin_sym: complex, hin_sym: complex):
    """
    Correct depth-2 recurrence in Cha symbol domain (FULL grid).

    Inputs
    ------
    delta : float
        1.0 if match, 0.0 if mismatch (eq_from_chars output; make sure it's clamped to [0,1])
    vin_sym : complex
        symbol for V(i, j-1)
    hin_sym : complex
        symbol for H(i-1, j)

    Returns
    -------
    (vout_sym, hout_sym) : complex, complex
        symbols for V(i, j) and H(i, j)
    """
    I = 1.0j
    P1 = g_map(1)  # symbol for +1

    y = vin_sym
    z = hin_sym

    # A,B are boolean indicators:
    # g(-1)=0.5 (real) => y+conj(y)=1
    # g(0)=0 or g(+1)=0.5j => y+conj(y)=0
    A = y + np.conj(y)  # 1 if vin == -1 else 0
    B = z + np.conj(z)  # 1 if hin == -1 else 0

    # Match case (delta=1): vout=-hin, hout=-vin in integer domain
    v_match = I * np.conj(z)
    h_match = I * np.conj(y)

    # Mismatch case (delta=0):
    # - if vin==-1: vout=-hin, hout=+1
    # - elif hin==-1: vout=+1, hout=-vin
    # - else (vin,hin in {0,1}): vout=1-hin, hout=1-vin
    v_mis = A * v_match + (1.0 - A) * (B * P1 + (1.0 - B) * (P1 - z))
    h_mis = A * P1     + (1.0 - A) * (B * h_match + (1.0 - B) * (P1 - y))

    # Select by delta (still branchless)
    vout = delta * v_match + (1.0 - delta) * v_mis
    hout = delta * h_match + (1.0 - delta) * h_mis
    return vout, hout


# ------------------------------------------------------------
# 5) FULL Myers DP in anti-diagonal order (symbol domain)
def full_myers_symbol(a: str, b: str):
    """
    Full-grid anti-diagonal Myers in Cha symbol domain.

    Computes Cha's complex symbol-mapped depth-2 recurrence on the FULL (L+1)x(L+1) grid.
    Since this is full DP, every interior cell (i>0, j>0) has both neighbors on k-1,
    so there is no "missing" handling.
    """
    L = len(a)
    assert len(b) == L

    i0s: list[int] = []
    Hdiags: list[np.ndarray] = []
    Vdiags: list[np.ndarray] = []

    Z = 0.0 + 0.0j
    P1 = g_map(1)  # +1 symbol

    for k in range(0, 2 * L + 1):
        rng = full_diag_indices(L, k)
        if rng is None:
            i0s.append(-1)
            Hdiags.append(np.zeros((0,), dtype=np.complex128))
            Vdiags.append(np.zeros((0,), dtype=np.complex128))
            continue

        i0, i1 = rng
        n = i1 - i0 + 1
        Hcur = np.zeros((n,), dtype=np.complex128)
        Vcur = np.zeros((n,), dtype=np.complex128)

        # Previous diagonal (k-1)
        if k - 1 >= 0:
            prev_i0 = i0s[k - 1]
            Vprev = Vdiags[k - 1]
            Hprev = Hdiags[k - 1]
        else:
            prev_i0 = -1
            Vprev = Hprev = None  # type: ignore[assignment]

        for t in range(n):
            i = i0 + t
            j = k - i

            if i == 0 and j == 0:
                Hcur[t] = Z
                Vcur[t] = Z
                continue
            if i == 0:
                Hcur[t] = P1
                Vcur[t] = Z
                continue
            if j == 0:
                Hcur[t] = Z
                Vcur[t] = P1
                continue

            # Interior neighbors are guaranteed to exist on the full grid.
            vin = Vprev[i - prev_i0]          # (i, j-1) on k-1
            hin = Hprev[(i - 1) - prev_i0]    # (i-1, j) on k-1

            delta = eq_from_chars(a[i - 1], b[j - 1])
            vout, hout = myers_step_symbolmapped_sym(delta, vin, hin)
            Vcur[t] = vout
            Hcur[t] = hout

        i0s.append(i0) 
        Hdiags.append(Hcur)
        Vdiags.append(Vcur)

    return i0s, Hdiags, Vdiags


# ------------------------------------------------------------
# 6) Optional integer reconstruction + distance from V last column
def quantize_diags(diags_sym):
    out = []
    for arr in diags_sym:
        out.append(np.array([g_inv(z) for z in arr], dtype=np.int8))
    return out


def recover_dist_from_V_int(L, i0s, Vdiags_int):
    # dist = L + sum_{i=1..L} V[i, L]
    s = 0
    for i in range(1, L + 1):
        k = i + L
        if k < 0 or k >= len(Vdiags_int):
            return None
        i0 = i0s[k]
        if i0 < 0:
            return None
        t = i - i0
        arr = Vdiags_int[k]
        if t < 0 or t >= len(arr):
            return None
        j = k - i
        if j != L:
            return None
        s += int(arr[t])
    return int(L + s)


# ------------------------------------------------------------
# 7) Save NPZ with object arrays of diags
def save_npz(path: Path, L: int, pid: str, i0s, Hsym, Vsym, emit_int: bool, dist_true=None):
    Hsym_obj = np.empty((len(Hsym),), dtype=object)
    Vsym_obj = np.empty((len(Vsym),), dtype=object)
    for k in range(len(Hsym)):
        Hsym_obj[k] = Hsym[k]
        Vsym_obj[k] = Vsym[k]

    payload = dict(
        pid=pid,
        L=np.int32(L),
        mode="full",
        i0=np.array(i0s, dtype=np.int32),
        Hsym=Hsym_obj,
        Vsym=Vsym_obj,
    )

    if dist_true is not None:
        payload["dist_true"] = np.int32(dist_true)

    if emit_int:
        H_int = quantize_diags(Hsym)
        V_int = quantize_diags(Vsym)

        H_obj = np.empty((len(H_int),), dtype=object)
        V_obj = np.empty((len(V_int),), dtype=object)
        for k in range(len(H_int)):
            H_obj[k] = H_int[k]
            V_obj[k] = V_int[k]

        dist_sym = recover_dist_from_V_int(L, i0s, V_int)
        payload["H"] = H_obj
        payload["V"] = V_obj
        payload["dist"] = np.int32(dist_sym if dist_sym is not None else -1)
        payload["dist_sym"] = np.int32(dist_sym if dist_sym is not None else -1)

    np.savez_compressed(path, **payload)


# ------------------------------------------------------------
# 8) CLI
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", type=Path, required=True, help="pair FASTA file (4-line blocks)")
    ap.add_argument("-k", "--pair", type=int, default=0, help="pair index")
    ap.add_argument("-p", "--prefix", type=str, default="out", help="output prefix")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("out"), help="output dir")
    ap.add_argument("--no_emit_int", action="store_true", help="do NOT emit integer H/V and dist in output npz")

    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    pid, a, b = get_pair(args.input, args.pair)
    if len(a) != len(b):
        raise SystemExit(f"Lengths differ: len(a)={len(a)} len(b)={len(b)}")

    L = len(a)
    dist_true = levenshtein_dp(a, b)
    i0s, Hsym, Vsym = full_myers_symbol(a, b)

    # If the recurrence still blows up, you'll get NaNs. Dist recovery from symbols becomes meaningless.
    # Keep going (so NPZ is saved) but surface it in the console output.
    has_nan = False
    for arr in Hsym:
        if np.isnan(arr).any():
            has_nan = True
            break
    if not has_nan:
        for arr in Vsym:
            if np.isnan(arr).any():
                has_nan = True
                break

    out_path = args.outdir / f"{args.prefix}.myersHVsym.full.npz"
    save_npz(out_path, L, pid, i0s, Hsym, Vsym, emit_int=(not args.no_emit_int), dist_true=dist_true)

    if not args.no_emit_int:
        V_int = quantize_diags(Vsym)
        dist_sym = recover_dist_from_V_int(L, i0s, V_int)
        if has_nan:
            print(f"[warn] pid={pid} L={L} FULL has_nan=True dist_sym={dist_sym} dist_true={dist_true}")
        else:
            print(f"[ok] pid={pid} L={L} FULL dist_sym={dist_sym} dist_true={dist_true}")
    else:
        if has_nan:
            print(f"[warn] pid={pid} L={L} FULL (symbol-only) has_nan=True dist_true={dist_true}")
        else:
            print(f"[ok] pid={pid} L={L} FULL (symbol-only) dist_true={dist_true}")

    print(f"wrote: {out_path}")


if __name__ == "__main__":
    main()
