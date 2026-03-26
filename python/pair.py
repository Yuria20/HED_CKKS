import os
import pysam

def rc(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def n2a(seq: str) -> str:
    # N/n -> A 로 강제 치환
    return seq.replace("N", "A").replace("n", "A")

BAM = "out/alnq.bam"
REF = "out/chr22chr.fa"   # chr22만 쓰는 현재 파이프라인 기준
OUT = "out/pairs.fa"

MINL = int(os.environ.get("MINL", "150"))
MAXL = int(os.environ.get("MAXL", "2048"))

bam = pysam.AlignmentFile(BAM, "rb")
ref = pysam.FastaFile(REF)

n = 0
with open(OUT, "w") as f:
    for r in bam.fetch(until_eof=True):
        if r.query_sequence is None:
            continue

        read_seq = n2a(r.query_sequence.upper())
        L = len(read_seq)
        if L < MINL or L > MAXL:
            continue

        chrom = bam.get_reference_name(r.reference_id)
        start = r.reference_start
        ref_seq = n2a(ref.fetch(chrom, start, start + L).upper())

        if r.is_reverse:
            ref_seq = rc(ref_seq)

        # 안전장치: 혹시라도 남아있으면 마지막에 한번 더 치환
        read_seq = n2a(read_seq)
        ref_seq = n2a(ref_seq)

        f.write(f">{r.query_name}\n{read_seq}\n")
        f.write(f"<{chrom}:{start}-{start+L}\n{ref_seq}\n")
        n += 1

print("pairs:", n, "->", OUT)