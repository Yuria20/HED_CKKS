#!/usr/bin/env python3
# levenshtein_fast.py
#
# Fast plaintext edit distance (Levenshtein) using Myers' bit-parallel algorithm.
# Works well for strings up to a few thousand chars in pure Python thanks to big-int bit ops.
#
# Usage:
#   python3 levenshtein_fast.py --a ACGT... --b ACGT...
#   python3 levenshtein_fast.py -i out/pairsrnd.fa -k 0
#
# Notes:
# - This computes the exact Levenshtein distance (unit costs: ins/del/sub = 1).
# - For DNA alphabets (ACGTN) it's extremely fast.
# - For very long strings (>> 10^5) Python may slow down; for your range (<=2048) it's fine.

import argparse
from pathlib import Path

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

def myers_levenshtein(a: str, b: str) -> int:
    """Exact Levenshtein distance using Myers bit-parallel algorithm (1999).
    Complexity ~ O(len(b) * ceil(len(a)/word_bits)), but with Python big-int it's fast.
    """
    # Ensure a is the pattern (columns), b is text (rows)
    # Prefer smaller pattern for speed
    if len(a) > len(b):
        a, b = b, a

    m = len(a)
    n = len(b)
    if m == 0:
        return n
    if n == 0:
        return m

    # Bitmask for m bits
    mask = (1 << m) - 1
    high_bit = 1 << (m - 1)

    # Peq: char -> bitmask of positions where pattern has that char
    Peq = {}
    for i, ch in enumerate(a):
        Peq[ch] = Peq.get(ch, 0) | (1 << i)

    VP = mask  # all 1s within m bits
    VN = 0
    D = m

    for ch in b:
        PM = Peq.get(ch, 0)
        # Core Myers update
        X = PM | VN
        # D0 uses integer addition; keep within m bits with &mask where needed
        D0 = (((X & VP) + VP) ^ VP) | X
        D0 &= mask

        HP = VN | ~(D0 | VP)
        HN = VP & D0

        # update distance using highest bit
        if HP & high_bit:
            D += 1
        elif HN & high_bit:
            D -= 1

        HP = (HP << 1) | 1
        HN = (HN << 1)

        VP = (HN | ~(D0 | HP)) & mask
        VN = (HP & D0) & mask

    return D

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", type=str, default=None, help="string a")
    ap.add_argument("--b", type=str, default=None, help="string b")
    ap.add_argument("-i", "--input", type=Path, default=None, help="pair FASTA file (4-line blocks)")
    ap.add_argument("-k", "--pair", type=int, default=0, help="pair index in FASTA")
    args = ap.parse_args()

    if args.input is not None:
        pid, a, b = get_pair(args.input, args.pair)
        dist = myers_levenshtein(a, b)
        print(f"pid={pid} len={len(a)} dist={dist}")
        return

    if args.a is None or args.b is None:
        raise SystemExit("Provide either --a/--b or -i/--input")

    print(myers_levenshtein(args.a, args.b))

if __name__ == "__main__":
    main()
