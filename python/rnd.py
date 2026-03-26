#!/usr/bin/env python3
# generate_dataset_py.py
import argparse, math, random
from pathlib import Path

ALPHABET = "ACGT"

def rand_iid(lo: int, hi: int) -> int:
    # C 코드 rand_iid(0, n)처럼 [lo, hi) 균등
    return random.randrange(lo, hi)

def sequence_generate_random(length: int) -> str:
    return "".join(ALPHABET[rand_iid(0, len(ALPHABET))] for _ in range(length))

def sequence_extract(seqlong: str, seqshort_length: int) -> tuple[str, int]:
    # seqlong에서 seqshort_length만큼 substring 추출 + offset 반환
    max_off = len(seqlong) - seqshort_length
    offset = rand_iid(0, max_off + 1)  # inclusive max_off
    return seqlong[offset:offset + seqshort_length], offset

def gen_mismatch(seq: list[str], errors_log: list[tuple]) -> None:
    pos = rand_iid(0, len(seq))
    base = seq[pos]
    repl = base
    while repl == base:
        repl = ALPHABET[rand_iid(0, len(ALPHABET))]
    seq[pos] = repl
    errors_log.append(("M", pos, base, repl))

def gen_deletion(seq: list[str], errors_log: list[tuple]) -> None:
    if not seq:
        return
    pos = rand_iid(0, len(seq))
    seq.pop(pos)
    errors_log.append(("D", pos))

def gen_insertion(seq: list[str], errors_log: list[tuple]) -> None:
    pos = rand_iid(0, len(seq) + 1)
    ch = ALPHABET[rand_iid(0, len(ALPHABET))]
    seq.insert(pos, ch)
    errors_log.append(("I", pos))

def sequence_generate_errors(seq: str, num_errors: int) -> tuple[str, list[tuple]]:
    s = list(seq)
    log = []
    for _ in range(num_errors):
        t = rand_iid(0, 3)  # 0,1,2 균등
        if t == 0:
            gen_mismatch(s, log)
        elif t == 1:
            gen_deletion(s, log)
        else:
            gen_insertion(s, log)
    return "".join(s), log

def gen_large_deletion(seq: list[str], deletion_length: int, errors_log: list[tuple]) -> None:
    if deletion_length <= 0 or deletion_length >= len(seq):
        return
    pos = rand_iid(0, len(seq))
    # pos부터 deletion_length만큼 “삭제 효과”
    # (끝 넘어가면 실제로는 가능한 만큼만 없어짐)
    end = min(pos + deletion_length, len(seq))
    del seq[pos:end]
    errors_log.append(("N", pos, deletion_length))

def sequence_generate_indels(seq: str, max_indels: int, deletion_length: int) -> tuple[str, list[tuple]]:
    s = list(seq)
    log = []
    indels = rand_iid(0, max_indels + 1)  # 0..max_indels
    for _ in range(indels):
        gen_large_deletion(s, deletion_length, log)
    return "".join(s), log

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", default=None)
    ap.add_argument("-n", "--num-patterns", type=int, default=1000, dest="num_reads")
    ap.add_argument("-l", "--length", type=int, default=100, help="pattern length")
    ap.add_argument("--length-diff", type=float, default=1.0)
    ap.add_argument("-e", "--error", type=float, default=0.04, dest="error_degree",
                    help="(<1) ratio of pattern length, (>=1) absolute #errors")
    ap.add_argument("--indels", default="0,0", help="NUM,LENGTH (large deletions)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    random.seed(args.seed)

    indels_num, indels_len = args.indels.split(",")
    indels_num = int(float(indels_num))
    indels_len = int(float(indels_len))

    pattern_length = args.length
    text_length = int(math.ceil(pattern_length * args.length_diff))
    seqlong_length = max(pattern_length, text_length)
    seqshort_length = min(pattern_length, text_length)

    if args.error_degree >= 1.0:
        num_errors = int(args.error_degree)
    else:
        num_errors = int(math.ceil(pattern_length * args.error_degree))

    out = open(args.output, "w") if args.output else None
    f = out if out else __import__("sys").stdout

    for i in range(args.num_reads):
        seqlong = sequence_generate_random(seqlong_length)
        seqshort, offset = sequence_extract(seqlong, seqshort_length)

        errors_log = []
        if num_errors > 0:
            seqshort, log = sequence_generate_errors(seqshort, num_errors)
            errors_log.extend(log)

        if indels_num > 0 and indels_len > 0:
            seqshort, log = sequence_generate_indels(seqshort, indels_num, indels_len)
            errors_log.extend(log)

        if args.debug:
            f.write(f"#DEBUG i={i} offset={offset} errors={errors_log}\n")

        # levenshtein_fast.py의 iter_pairs()는:
        # >pid\nREAD\n<ref\nREF\n
        f.write(f">r{i}\n{seqshort}\n<ref{i}\n{seqlong}\n")

    if out:
        out.close()

if __name__ == "__main__":
    main()