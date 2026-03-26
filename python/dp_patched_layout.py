#!/usr/bin/env python3
# dp.py
#
# CKKS-ready variant of dp.py:
# - Computes Myers H/V using Cha's complex symbol-mapped depth-2 recurrence.
# - NO per-cell inverse mapping (g_inv) and NO data-dependent branching (delta>=0.5).
# - Stores complex-valued symbol diags (Hsym/Vsym) so you can pack/encrypt them.
# - Optionally (default) does a SINGLE inverse mapping at the very end (plaintext-only)
#   to produce integer H/V (in {-1,0,1}) for debugging / compatibility with pack.py.
#
# Usage (random pairs):
#   python3 dp.py -i out/pairsrnd.fa -p rnd -w 3 -k 0 -o out
#
# Notes:
# - This file uses a *windowed anti-diagonal band*: for each anti-diagonal k=i+j,
#   we keep up to `w` consecutive cells centered around i≈j. For w=3 this matches
#   the "3 cells per diagonal, sliding by one" pattern you drew.
# - Some interior cells may miss needed neighbors because left/up dependency can
#   fall outside the kept window. Default CLI is --missing pad (fill with (+1)-symbol).
#   You can also use --missing zero (0-symbol) or --missing error to hard-fail.
#
# Output:
#   out/<pre>.myersHVsym.npz  (contains Hsym/Vsym always; H/V ints if emit_int)

import argparse
from pathlib import Path
import numpy as np

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
    return float(eq_poly(z))


# ------------------------------------------------------------
# 3) Band index helper (windowed anti-diagonal; up to w cells per k=i+j)
def band_diag_indices(L, w, k):
    """Return the i-range (i0,i1) for the selected cells on anti-diagonal k=i+j.

    This file uses a *windowed anti-diagonal band*:
      - For each diagonal k, we keep up to `w` consecutive cells (contiguous i-range),
        centered around i≈j (i≈k/2).
      - If the diagonal itself is shorter than w near the corners, we keep the whole
        available diagonal segment.

    This matches the “3 cells per diagonal, sliding by one” pattern you drew when w=3.
    """
    i_min = max(0, k - L)   # smallest i on this diagonal within the (L+1)x(L+1) grid
    i_max = min(L, k)       # largest  i on this diagonal within the grid
    if i_min > i_max:
        return None

    diag_len = i_max - i_min + 1
    if w is None or w <= 0:
        return (i_min, i_max)
    if diag_len <= w:
        return (i_min, i_max)

    # Choose a length-w window centered around i≈k/2 (equivalently j≈k/2)
    offset = (w - 1) // 2
    center = k // 2
    i0 = center - offset

    # Clamp window to the available diagonal range
    if i0 < i_min:
        i0 = i_min
    if i0 + w - 1 > i_max:
        i0 = i_max - (w - 1)

    i1 = i0 + w - 1
    return (i0, i1)

# ------------------------------------------------------------
# 3b) Packing helpers for HEaaN-friendly even/odd buffers
#
# Important distinction:
#   - band_w   : actual algorithmic band width used by DP recurrence
#   - layout_w : fixed HE slot width reserved for the data region inside one block
#                (can be >= band_w so you can sweep different w values while keeping
#                 the same ciphertext/block layout)
#   - guard    : number of guard slots on each side of the data region
#
# When layout_w == band_w and guard == 0, this reduces to the original behavior.

def _pack_to_block(arr, layout_w: int, guard: int = 0, order: str = "desc", pad_value=0.0 + 0.0j):
    """Pack one variable-length diagonal into a fixed HE block.

    Block layout is:
        [ guard | data(layout_w slots) | guard ]

    The real diagonal payload is left-justified inside the data region.
    This makes the block shape independent of the actual diagonal length n<=band_w,
    so the same rotate schedule can be reused even near the corners.
    """
    n = int(len(arr))
    B = int(layout_w + 2 * guard)
    out = np.full((B,), pad_value, dtype=np.complex128)
    valid = np.zeros((B,), dtype=np.uint8)
    if n <= 0:
        return out, valid
    if n > layout_w:
        raise ValueError(f"diag length n={n} exceeds layout_w={layout_w}")
    if order not in ("asc", "desc"):
        raise ValueError("order must be 'asc' or 'desc'")
    seq = arr if order == "asc" else arr[::-1]
    start = guard
    out[start:start + n] = seq
    valid[start:start + n] = 1
    return out, valid


def _pack_mask_to_block(mask_arr, layout_w: int, guard: int = 0, order: str = "desc"):
    n = int(len(mask_arr))
    B = int(layout_w + 2 * guard)
    out = np.zeros((B,), dtype=np.uint8)
    if n <= 0:
        return out
    if n > layout_w:
        raise ValueError(f"mask length n={n} exceeds layout_w={layout_w}")
    if order not in ("asc", "desc"):
        raise ValueError("order must be 'asc' or 'desc'")
    seq = mask_arr if order == "asc" else mask_arr[::-1]
    start = guard
    out[start:start + n] = seq.astype(np.uint8, copy=False)
    return out


def make_block_masks(layout_w: int, band_w: int | None = None, guard: int = 0):
    """Return canonical per-block masks for HE runtime.

    data_mask_block:
        1 on the reserved data region (layout_w slots), 0 on guards.
        Apply this after block-local ±1 rotations to kill guard spillover.

    active_mask_block:
        1 only on the actually used band width for this run.
        Useful when you keep a fixed layout_w but run with a smaller band_w.

    left_guard_block / right_guard_block:
        Convenience masks for debugging and for explicit zeroing/fill on one side.
    """
    if band_w is None:
        band_w = layout_w
    if band_w > layout_w:
        raise ValueError(f"band_w={band_w} exceeds layout_w={layout_w}")

    B = int(layout_w + 2 * guard)
    data_mask = np.zeros((B,), dtype=np.uint8)
    active_mask = np.zeros((B,), dtype=np.uint8)
    left_guard = np.zeros((B,), dtype=np.uint8)
    right_guard = np.zeros((B,), dtype=np.uint8)
    start = guard
    data_mask[start:start + layout_w] = 1
    active_mask[start:start + band_w] = 1
    if guard > 0:
        left_guard[:guard] = 1
        right_guard[-guard:] = 1
    return {
        "block_size": np.int32(B),
        "block_stride": np.int32(B),
        "layout_w": np.int32(layout_w),
        "band_w": np.int32(band_w),
        "guard": np.int32(guard),
        "data_mask_block": data_mask,
        "active_mask_block": active_mask,
        "left_guard_block": left_guard,
        "right_guard_block": right_guard,
    }


def tile_block_pattern(block_arr, nblocks: int):
    """Repeat one block-shaped vector/mask across multiple HE blocks."""
    if nblocks <= 0:
        raise ValueError("nblocks must be >= 1")
    arr = np.asarray(block_arr)
    return np.tile(arr, int(nblocks))


def tile_object_blocks(obj_arr, nblocks: int):
    """Tile every object-array block independently (useful for valid/mv/mh/etc.)."""
    if nblocks <= 0:
        raise ValueError("nblocks must be >= 1")
    out = np.empty((len(obj_arr),), dtype=object)
    for i, x in enumerate(obj_arr):
        out[i] = tile_block_pattern(x, nblocks)
    return out


def _diag_runtime_masks(k: int, i0: int, n: int,
                        mv_raw=None, mh_raw=None):
    """Build raw per-diagonal masks before block packing.

    Returns raw length-n masks for:
      valid, origin, top, left, interior
      and optionally mv, mh, padv, padh where
        padv = interior & ~mv
        padh = interior & ~mh

    The pad* masks are the important guard-aware part:
    do NOT use (1-mv) directly in HE runtime, because that would also fill guards
    and inactive layout slots. Use padv/padh instead.
    """
    valid = np.ones((n,), dtype=np.uint8)
    origin = np.zeros((n,), dtype=np.uint8)
    top = np.zeros((n,), dtype=np.uint8)
    left = np.zeros((n,), dtype=np.uint8)
    interior = np.zeros((n,), dtype=np.uint8)

    for t in range(n):
        i = i0 + t
        j = k - i
        if i == 0 and j == 0:
            origin[t] = 1
        elif i == 0:
            top[t] = 1
        elif j == 0:
            left[t] = 1
        else:
            interior[t] = 1

    out = {
        "valid": valid,
        "origin": origin,
        "top": top,
        "left": left,
        "interior": interior,
    }

    if mv_raw is not None and mh_raw is not None:
        mv = np.asarray(mv_raw, dtype=np.uint8)
        mh = np.asarray(mh_raw, dtype=np.uint8)
        if len(mv) != n or len(mh) != n:
            raise ValueError("mv/mh length mismatch with diagonal length")
        padv = (interior & (1 - mv)).astype(np.uint8, copy=False)
        padh = (interior & (1 - mh)).astype(np.uint8, copy=False)
        out["mv"] = mv
        out["mh"] = mh
        out["padv"] = padv
        out["padh"] = padh

    return out


def build_packed_even_odd(i0s, Hdiags, Vdiags, mv_diags=None, mh_diags=None,
                          band_w: int = 7, layout_w: int | None = None,
                          guard: int = 0, order: str = "asc", tile_blocks: int = 1):
    """Split and pack per-diagonal arrays into fixed HE blocks, separated by parity of k.

    band_w:
        Actual algorithmic band width used when producing Hdiags/Vdiags.

    layout_w:
        Fixed slot width reserved inside each HE block. If larger than band_w,
        you can reuse the same layout for multiple band widths <= layout_w.

    order:
        For HE runtime, use 'asc'. With 'desc', packed values are still valid for storage,
        but the simple local-rotation schedule (delta_i0 / vin/hin offsets) is no longer
        shape-invariant near the corners.

    tile_blocks:
        Number of identical block slots you plan to tile across one ciphertext when batching
        multiple same-(L,w) pairs. Only canonical tiled masks are emitted; H/V values remain
        one-block-per-entry because this file stores a single pair.

    Returns a dict that can be embedded into NPZ payload.
    In addition to H/V and valid/mv/mh, it now emits guard-aware runtime masks:
      - origin/top/left/interior
      - padv/padh   (where to inject pad after local rotate)
      - delta_i0 and local source-offset schedules for vin/hin/A/Brev (asc only)
    """
    if layout_w is None:
        layout_w = band_w
    if layout_w < band_w:
        raise ValueError(f"layout_w={layout_w} must be >= band_w={band_w}")
    if tile_blocks <= 0:
        raise ValueError("tile_blocks must be >= 1")

    ks_even, ks_odd = [], []
    H_even, V_even, valid_even = [], [], []
    H_odd, V_odd, valid_odd = [], [], []

    origin_even, top_even, left_even, interior_even = [], [], [], []
    origin_odd, top_odd, left_odd, interior_odd = [], [], [], []

    mv_even, mh_even, padv_even, padh_even = [], [], [], []
    mv_odd, mh_odd, padv_odd, padh_odd = [], [], [], []

    have_mv = mv_diags is not None
    have_mh = mh_diags is not None
    if have_mv != have_mh:
        raise ValueError("mv_diags and mh_diags must either both be provided or both be None")

    Z = 0.0 + 0.0j
    K = len(Hdiags)

    diag_len = np.array([len(Hdiags[k]) for k in range(K)], dtype=np.int32)
    delta_i0 = np.zeros((K,), dtype=np.int32)
    if K >= 2:
        for k in range(1, K):
            if i0s[k] >= 0 and i0s[k - 1] >= 0:
                delta_i0[k] = np.int32(i0s[k] - i0s[k - 1])

    runtime_ready = np.uint8(1 if order == "asc" else 0)
    if runtime_ready:
        vin_src_offset = delta_i0.astype(np.int32, copy=True)
        hin_src_offset = (delta_i0 - 1).astype(np.int32, copy=False)
        A_src_offset = delta_i0.astype(np.int32, copy=True)
        Brev_src_offset = (delta_i0 - 1).astype(np.int32, copy=False)
        max_abs_local_shift = np.int32(
            max(0,
                int(np.max(np.abs(vin_src_offset))),
                int(np.max(np.abs(hin_src_offset))),
                int(np.max(np.abs(A_src_offset))),
                int(np.max(np.abs(Brev_src_offset))))
        )
        guard_sufficient = np.uint8(1 if guard >= int(max_abs_local_shift) else 0)
    else:
        vin_src_offset = None
        hin_src_offset = None
        A_src_offset = None
        Brev_src_offset = None
        max_abs_local_shift = np.int32(-1)
        guard_sufficient = np.uint8(0)

    for k in range(K):
        Hk = Hdiags[k]
        Vk = Vdiags[k]
        n = int(len(Hk))
        i0 = int(i0s[k]) if k < len(i0s) else -1

        packedH, valid = _pack_to_block(Hk, layout_w=layout_w, guard=guard, order=order, pad_value=Z)
        packedV, _ = _pack_to_block(Vk, layout_w=layout_w, guard=guard, order=order, pad_value=Z)

        rt = _diag_runtime_masks(
            k=k,
            i0=i0,
            n=n,
            mv_raw=(mv_diags[k] if have_mv else None),
            mh_raw=(mh_diags[k] if have_mh else None),
        )

        packed_origin = _pack_mask_to_block(rt["origin"], layout_w=layout_w, guard=guard, order=order)
        packed_top = _pack_mask_to_block(rt["top"], layout_w=layout_w, guard=guard, order=order)
        packed_left = _pack_mask_to_block(rt["left"], layout_w=layout_w, guard=guard, order=order)
        packed_interior = _pack_mask_to_block(rt["interior"], layout_w=layout_w, guard=guard, order=order)

        if have_mv:
            packed_mv = _pack_mask_to_block(rt["mv"], layout_w=layout_w, guard=guard, order=order)
            packed_mh = _pack_mask_to_block(rt["mh"], layout_w=layout_w, guard=guard, order=order)
            packed_padv = _pack_mask_to_block(rt["padv"], layout_w=layout_w, guard=guard, order=order)
            packed_padh = _pack_mask_to_block(rt["padh"], layout_w=layout_w, guard=guard, order=order)

        if (k % 2) == 0:
            ks_even.append(k)
            H_even.append(packedH)
            V_even.append(packedV)
            valid_even.append(valid)
            origin_even.append(packed_origin)
            top_even.append(packed_top)
            left_even.append(packed_left)
            interior_even.append(packed_interior)
            if have_mv:
                mv_even.append(packed_mv)
                mh_even.append(packed_mh)
                padv_even.append(packed_padv)
                padh_even.append(packed_padh)
        else:
            ks_odd.append(k)
            H_odd.append(packedH)
            V_odd.append(packedV)
            valid_odd.append(valid)
            origin_odd.append(packed_origin)
            top_odd.append(packed_top)
            left_odd.append(packed_left)
            interior_odd.append(packed_interior)
            if have_mv:
                mv_odd.append(packed_mv)
                mh_odd.append(packed_mh)
                padv_odd.append(packed_padv)
                padh_odd.append(packed_padh)

    def _to_obj(lst):
        obj = np.empty((len(lst),), dtype=object)
        for i, x in enumerate(lst):
            obj[i] = x
        return obj

    block_meta = make_block_masks(layout_w=layout_w, band_w=band_w, guard=guard)

    payload = dict(
        ks_even=np.array(ks_even, dtype=np.int32),
        ks_odd=np.array(ks_odd, dtype=np.int32),
        Hsym_even=_to_obj(H_even),
        Vsym_even=_to_obj(V_even),
        valid_even=_to_obj(valid_even),
        Hsym_odd=_to_obj(H_odd),
        Vsym_odd=_to_obj(V_odd),
        valid_odd=_to_obj(valid_odd),
        origin_even=_to_obj(origin_even),
        top_even=_to_obj(top_even),
        left_even=_to_obj(left_even),
        interior_even=_to_obj(interior_even),
        origin_odd=_to_obj(origin_odd),
        top_odd=_to_obj(top_odd),
        left_odd=_to_obj(left_odd),
        interior_odd=_to_obj(interior_odd),
        diag_len=diag_len,
        delta_i0=delta_i0,
        pack_order=str(order),
        runtime_ready=runtime_ready,
        max_abs_local_shift=max_abs_local_shift,
        guard_sufficient=guard_sufficient,
        tile_blocks=np.int32(tile_blocks),
        **block_meta,
    )

    if runtime_ready:
        payload["vin_src_offset"] = vin_src_offset
        payload["hin_src_offset"] = hin_src_offset
        payload["A_src_offset"] = A_src_offset
        payload["Brev_src_offset"] = Brev_src_offset

    if tile_blocks > 1:
        payload["data_mask_tiled"] = tile_block_pattern(block_meta["data_mask_block"], tile_blocks)
        payload["active_mask_tiled"] = tile_block_pattern(block_meta["active_mask_block"], tile_blocks)
        payload["left_guard_tiled"] = tile_block_pattern(block_meta["left_guard_block"], tile_blocks)
        payload["right_guard_tiled"] = tile_block_pattern(block_meta["right_guard_block"], tile_blocks)
        # Guard-aware runtime masks repeated across same-(L,w) blocks.
        # Useful when one ciphertext batches multiple independent pairs.
        payload["valid_even_tiled"] = tile_object_blocks(payload["valid_even"], tile_blocks)
        payload["valid_odd_tiled"] = tile_object_blocks(payload["valid_odd"], tile_blocks)
        payload["origin_even_tiled"] = tile_object_blocks(payload["origin_even"], tile_blocks)
        payload["origin_odd_tiled"] = tile_object_blocks(payload["origin_odd"], tile_blocks)
        payload["top_even_tiled"] = tile_object_blocks(payload["top_even"], tile_blocks)
        payload["top_odd_tiled"] = tile_object_blocks(payload["top_odd"], tile_blocks)
        payload["left_even_tiled"] = tile_object_blocks(payload["left_even"], tile_blocks)
        payload["left_odd_tiled"] = tile_object_blocks(payload["left_odd"], tile_blocks)
        payload["interior_even_tiled"] = tile_object_blocks(payload["interior_even"], tile_blocks)
        payload["interior_odd_tiled"] = tile_object_blocks(payload["interior_odd"], tile_blocks)

    if have_mv and have_mh:
        payload["mv_even"] = _to_obj(mv_even)
        payload["mh_even"] = _to_obj(mh_even)
        payload["padv_even"] = _to_obj(padv_even)
        payload["padh_even"] = _to_obj(padh_even)
        payload["mv_odd"] = _to_obj(mv_odd)
        payload["mh_odd"] = _to_obj(mh_odd)
        payload["padv_odd"] = _to_obj(padv_odd)
        payload["padh_odd"] = _to_obj(padh_odd)
        if tile_blocks > 1:
            payload["mv_even_tiled"] = tile_object_blocks(payload["mv_even"], tile_blocks)
            payload["mh_even_tiled"] = tile_object_blocks(payload["mh_even"], tile_blocks)
            payload["padv_even_tiled"] = tile_object_blocks(payload["padv_even"], tile_blocks)
            payload["padh_even_tiled"] = tile_object_blocks(payload["padh_even"], tile_blocks)
            payload["mv_odd_tiled"] = tile_object_blocks(payload["mv_odd"], tile_blocks)
            payload["mh_odd_tiled"] = tile_object_blocks(payload["mh_odd"], tile_blocks)
            payload["padv_odd_tiled"] = tile_object_blocks(payload["padv_odd"], tile_blocks)
            payload["padh_odd_tiled"] = tile_object_blocks(payload["padh_odd"], tile_blocks)

    return payload



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
    # delta ~1 when match, ~0 when mismatch, but keep it real-valued (no thresholding)
    x = -delta
    y = vin_sym
    z = hin_sym
    hout = H_tilde(x, y, z)
    vout = V_tilde(x, y, z)
    return vout, hout


# ------------------------------------------------------------
# 5) Banded Myers DP in anti-diagonal order (symbol domain)
def banded_myers_symbol(a: str, b: str, w: int, missing: str = "zero", pad_alpha: float = 1.0, return_masks: bool = False):
    """Windowed anti-diagonal Myers in Cha symbol domain.

    This implementation uses *Option C* for missing-neighbor handling:
    it constructs the per-diagonal input vectors (vin/hin) and injects the chosen
    fill symbol directly at missing slots, avoiding any mask*ct "selection" step.

    missing:
      - 'zero' : fill missing vin/hin with 0-symbol
      - 'pad'  : fill missing vin/hin with (+1)-symbol g(+1)=i/2 scaled by pad_alpha (conservative boundary-like fill)
      - 'error': raise on missing vin/hin for interior cells (i>0, j>0)

    return_masks:
      If True, also return per-diagonal uint8 masks mv/mh indicating whether vin/hin
      existed (1) or was missing (0) for each stored cell on that diagonal.
      (Note: these masks are only length <= w per diagonal, not full DP-sized.)
    """
    L = len(a)
    assert len(b) == L

    i0s: list[int] = []
    Hdiags: list[np.ndarray] = []
    Vdiags: list[np.ndarray] = []

    mv_diags: list[np.ndarray] | None = [] if return_masks else None
    mh_diags: list[np.ndarray] | None = [] if return_masks else None

    Z = 0.0 + 0.0j
    P1 = g_map(1) * float(pad_alpha)  # +1 symbol (scaled)

    if missing not in ("zero", "pad", "error"):
        raise ValueError("missing must be one of {'zero','pad','error'}")

    fill = Z if missing == "zero" else P1  # used for missing='pad' and for provisional fill

    for k in range(0, 2 * L + 1):
        rng = band_diag_indices(L, w, k)
        if rng is None:
            i0s.append(-1)
            Hdiags.append(np.zeros((0,), dtype=np.complex128))
            Vdiags.append(np.zeros((0,), dtype=np.complex128))
            if return_masks:
                mv_diags.append(np.zeros((0,), dtype=np.uint8))  # type: ignore[union-attr]
                mh_diags.append(np.zeros((0,), dtype=np.uint8))  # type: ignore[union-attr]
            continue

        i0, i1 = rng
        n = i1 - i0 + 1

        Hcur = np.zeros((n,), dtype=np.complex128)
        Vcur = np.zeros((n,), dtype=np.complex128)

        # Prepare vin/hin vectors for this diagonal from the previous diagonal window.
        # For missing slots, inject the chosen fill symbol directly (Option C).
        vin_vec = np.full((n,), fill, dtype=np.complex128)
        hin_vec = np.full((n,), fill, dtype=np.complex128)

        if return_masks:
            mv_cur = np.zeros((n,), dtype=np.uint8)
            mh_cur = np.zeros((n,), dtype=np.uint8)
        else:
            mv_cur = mh_cur = None

        if k - 1 >= 0 and len(Vdiags) > 0:
            prev_i0 = i0s[k - 1]
            if prev_i0 >= 0:
                Vprev = Vdiags[k - 1]
                Hprev = Hdiags[k - 1]
                p_len = int(len(Vprev))

                for t in range(n):
                    i = i0 + t 
                    j = k - i

                    # Boundary cells do not use vin/hin; leave whatever is in vin_vec/hin_vec.
                    if i == 0 or j == 0:
                        continue

                    # vin uses previous diagonal at same i: (i, j-1) lies on k-1
                    u = i - prev_i0
                    if 0 <= u < p_len:
                        vin_vec[t] = Vprev[u]
                        if return_masks:
                            mv_cur[t] = 1
                    else:
                        if return_masks:
                            mv_cur[t] = 0

                    # hin uses previous diagonal at i-1: (i-1, j) lies on k-1
                    v = (i - 1) - prev_i0
                    if 0 <= v < p_len:
                        hin_vec[t] = Hprev[v]
                        if return_masks:
                            mh_cur[t] = 1
                    else:
                        if return_masks:
                            mh_cur[t] = 0

        # Now compute this diagonal.
        for t in range(n):
            i = i0 + t
            j = k - i

            if i == 0 and j == 0:
                Hcur[t] = Z
                Vcur[t] = Z
                continue

            if i == 0:
                # top row: H[0,j]=+1, V[0,j]=0
                Hcur[t] = g_map(1)
                Vcur[t] = Z
                continue

            if j == 0:
                # left col: V[i,0]=+1, H[i,0]=0
                Hcur[t] = Z
                Vcur[t] = g_map(1)
                continue

            # missing='error' means: for interior cells, both neighbors must exist in previous window
            if missing == "error":
                if return_masks:
                    if mv_cur[t] == 0 or mh_cur[t] == 0:
                        raise RuntimeError(
                            f"Missing neighbor at (i={i}, j={j}) on k={k}: mv={int(mv_cur[t])} mh={int(mh_cur[t])}"
                        )
                else:
                    prev_i0 = i0s[k - 1]
                    if prev_i0 < 0:
                        raise RuntimeError(f"Missing neighbor at (i={i}, j={j}) on k={k}: prev diag empty")
                    p_len = int(len(Vdiags[k - 1]))
                    p_lo = prev_i0
                    p_hi = prev_i0 + p_len - 1
                    if (i < p_lo or i > p_hi) or ((i - 1) < p_lo or (i - 1) > p_hi):
                        raise RuntimeError(
                            f"Missing neighbor at (i={i}, j={j}) on k={k}: i-range prev=[{p_lo},{p_hi}]"
                        )

            delta = eq_from_chars(a[i - 1], b[j - 1])
            vout, hout = myers_step_symbolmapped_sym(delta, vin_vec[t], hin_vec[t])
            Vcur[t] = vout
            Hcur[t] = hout

        i0s.append(i0)
        Hdiags.append(Hcur)
        Vdiags.append(Vcur)

        if return_masks:
            mv_diags.append(mv_cur)  # type: ignore[union-attr]
            mh_diags.append(mh_cur)  # type: ignore[union-attr]

    if return_masks:
        return i0s, Hdiags, Vdiags, mv_diags, mh_diags
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
def save_npz(path: Path, L: int, w: int, pid: str, i0s, Hsym, Vsym, emit_int: bool,
             packed_parity: dict | None = None):
    """Save NPZ.
    Always saves:
      - i0 (per k), Hsym/Vsym (object arrays of variable length diags)

    Optionally saves:
      - integer H/V and dist (plaintext-only debug)
      - packed even/odd fixed-length-w diags and masks (HEaaN-friendly)
    """
    Hsym_obj = np.empty((len(Hsym),), dtype=object)
    Vsym_obj = np.empty((len(Vsym),), dtype=object)
    for k in range(len(Hsym)):
        Hsym_obj[k] = Hsym[k]
        Vsym_obj[k] = Vsym[k]

    payload = dict(
        pid=pid,
        L=np.int32(L),
        w=np.int32(w),
        i0=np.array(i0s, dtype=np.int32),
        Hsym=Hsym_obj,
        Vsym=Vsym_obj,
    )

    if emit_int:
        H_int = quantize_diags(Hsym)
        V_int = quantize_diags(Vsym)

        H_obj = np.empty((len(H_int),), dtype=object)
        V_obj = np.empty((len(V_int),), dtype=object)
        for k in range(len(H_int)):
            H_obj[k] = H_int[k]
            V_obj[k] = V_int[k]

        dist = recover_dist_from_V_int(L, i0s, V_int)
        payload["H"] = H_obj
        payload["V"] = V_obj
        payload["dist"] = np.int32(dist if dist is not None else -1)

    if packed_parity is not None:
        # packed even/odd fields (object arrays)
        for kk, vv in packed_parity.items():
            payload[kk] = vv

    np.savez_compressed(path, **payload)


# ------------------------------------------------------------
# 8) CLI
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", type=Path, required=True, help="pair FASTA file (4-line blocks)")
    ap.add_argument("-k", "--pair", type=int, default=0, help="pair index")
    ap.add_argument("-w", "--band", type=int, default=7,
                    help="band width: number of cells per anti-diagonal (windowed around i≈j)")
    ap.add_argument("-p", "--prefix", type=str, default="out", help="output prefix")
    ap.add_argument("-o", "--outdir", type=Path, default=Path("out"), help="output dir")

    ap.add_argument("--missing", type=str, default="pad", choices=["zero", "pad", "error"],
                    help="how to handle missing neighbors due to band/window: zero-fill, pad-fill, or error")
    ap.add_argument("--pad_alpha", type=float, default=1.0,
                    help="scale for pad symbol when --missing pad (pad = alpha * g(+1) = alpha*(i/2))")

    ap.add_argument("--emit_packed", action="store_true",
                    help="also emit fixed-length-w packed even/odd diags + masks for HEaaN")
    ap.add_argument("--emit_masks", action="store_true",
                    help="deprecated: runtime masks are now emitted automatically with --emit_packed")

    ap.add_argument("--pack_order", type=str, default="asc", choices=["asc", "desc"],
                    help="slot order inside each packed diagonal (use 'asc' for HE runtime; 'desc' is mainly for drawing/debug)")
    ap.add_argument("--layout_w", type=int, default=None,
                    help="fixed HE layout width reserved for one packed diagonal block; defaults to --band")
    ap.add_argument("--guard", type=int, default=0,
                    help="number of guard slots on each side of the packed data region (HE block packing)")
    ap.add_argument("--tile_blocks", type=int, default=1,
                    help="optional runtime metadata for tiling one packed block across this many same-(L,w) pair blocks in one ciphertext")
    ap.add_argument("--no_emit_int", action="store_true",
                    help="do NOT emit integer H/V and dist in output npz (plaintext-only debug)")

    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    pid, a, b = get_pair(args.input, args.pair)
    if len(a) != len(b):
        raise SystemExit(f"Lengths differ: len(a)={len(a)} len(b)={len(b)}")

    L = len(a)
    layout_w = args.layout_w if args.layout_w is not None else args.band
    if layout_w < args.band:
        raise SystemExit(f"layout_w ({layout_w}) must be >= band ({args.band})")

    packed_parity = None
    if args.emit_packed:
        # Runtime packing now always emits mv/mh and derived padv/padh masks.
        # They are needed once guard slots are introduced; using (1-mv) directly would
        # incorrectly fill guards and inactive layout slots.
        i0s, Hsym, Vsym, mv, mh = banded_myers_symbol(
            a, b, args.band, missing=args.missing, pad_alpha=args.pad_alpha, return_masks=True
        )
        packed_parity = build_packed_even_odd(
            i0s, Hsym, Vsym,
            mv_diags=mv, mh_diags=mh,
            band_w=args.band, layout_w=layout_w, guard=args.guard, order=args.pack_order,
            tile_blocks=args.tile_blocks,
        )
    else:
        i0s, Hsym, Vsym = banded_myers_symbol(
            a, b, args.band, missing=args.missing, pad_alpha=args.pad_alpha, return_masks=False
        )
    out_path = args.outdir / f"{args.prefix}.myersHVsym.npz"
    save_npz(out_path, L, args.band, pid, i0s, Hsym, Vsym,
             emit_int=(not args.no_emit_int), packed_parity=packed_parity)

    if not args.no_emit_int:
        # quick sanity print (plaintext-only, after single g_inv)
        V_int = quantize_diags(Vsym)
        dist = recover_dist_from_V_int(L, i0s, V_int)
        print(f"[ok] pid={pid} L={L} w={args.band} dist={dist}")
    else:
        print(f"[ok] pid={pid} L={L} w={args.band} (symbol-only)")
    if args.emit_packed:
        runtime_ready = int(packed_parity.get("runtime_ready", 0)) if packed_parity is not None else 0
        guard_sufficient = int(packed_parity.get("guard_sufficient", 0)) if packed_parity is not None else 0
        print(f"[packed] order={args.pack_order} missing={args.missing} pad_alpha={args.pad_alpha} layout_w={layout_w} guard={args.guard} tile_blocks={args.tile_blocks}")
        print("         fields: ks_even/ks_odd, Hsym_even/Hsym_odd, Vsym_even/Vsym_odd, valid_*, origin/top/left/interior, mv/mh, padv/padh")
        print(f"         runtime_ready={runtime_ready} guard_sufficient={guard_sufficient}")
        if args.pack_order != "asc":
            print("         note: desc order disables the simple local rotate schedule; use --pack_order asc for HE runtime")

    print(f"wrote: {out_path}")


if __name__ == "__main__":
    main()