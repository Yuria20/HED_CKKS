package main

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const AlphabetSize = 256

type DeltaTable struct {
	DeltaBlocks []*rlwe.Ciphertext //eq 결과 ciphertext들
	KOfSlot     [][]int //각 slot이 어느 anti-diagonal k=i+j 에 속하는지
	PosInDiag   [][]int //해당 anti-diagonal 안에서 몇 번째 원소인지
	Valid       [][]bool //padding이 아닌 실제 데이터인지
	M           int //문자열 길이
	Slots       int //block 당 slot 수
	PairCount   int //전체 비교쌍 개수 = m*m
}

type PairIndex struct {
	IA  int // a의 index (0-based)
	JB  int // b의 index (0-based)
	K   int // anti-diagonal index (1-based 좌표계 기준 k=i+j)
	Pos int // diagonal 내부 위치
}

//복소원 매핑
func ByteToRoot(x int) complex128 {
	return complex(
		math.Cos(2*math.Pi*(float64(x)/float64(AlphabetSize))),
		math.Sin(2*math.Pi*(float64(x)/float64(AlphabetSize))),
	)
}

func DiagLen(k, m int) int {
	if k <= m+1 {
		return k - 1
	}
	return 2*m - k + 1
}

// 전개 순서:
// A: a1 | a1,a2 | a1,a2,a3 | ... | a_{m-1},a_m | a_m
// B: b1 | b2,b1 | b3,b2,b1 | ... | b_m,b_{m-1} | b_m
func buildIndexes(m int) (
	pairs []PairIndex,
	err error,
) {
	if m == 0 {
		return nil, fmt.Errorf("empty input")
	}

	// 위쪽 삼각형: k = 2 .. m+1
	for k := 2; k <= m+1; k++ {
		diagLen := k - 1
		for t := 0; t < diagLen; t++ {
			i := t + 1
			j := diagLen - t

			pairs = append(pairs, PairIndex{
				IA:  i - 1,
				JB:  j - 1,
				K:   k,
				Pos: t,
			})
		}
	}

	// 아래쪽 삼각형: k = m+2 .. 2m
	for k := m + 2; k <= 2*m; k++ {
		diagLen := 2*m - k + 1
		start := m - diagLen + 1

		for t := 0; t < diagLen; t++ {
			i := start + t
			j := m - t

			pairs = append(pairs, PairIndex{
				IA:  i - 1,
				JB:  j - 1,
				K:   k,
				Pos: t,
			})
		}
	}

	return
}

func packBlocks(
	pairs []PairIndex,
	slots int,
) (
	blocks [][]PairIndex,
	kOfSlot [][]int,
	posInDiag [][]int,
	valid [][]bool,
	err error,
) {
	if slots <= 0 {
		return nil, nil, nil, nil, fmt.Errorf("invalid slots: %d", slots)
	}

	n := len(pairs)
	for start := 0; start < n; start += slots {
		end := start + slots
		if end > n {
			end = n
		}

		curPairs := append([]PairIndex{}, pairs[start:end]...)
		curK := make([]int, 0, slots)
		curPos := make([]int, 0, slots)
		curValid := make([]bool, 0, slots)

		for _, p := range curPairs {
			curK = append(curK, p.K)
			curPos = append(curPos, p.Pos)
			curValid = append(curValid, true)
		}

		for len(curPairs) < slots {
			curPairs = append(curPairs, PairIndex{IA: -1, JB: -1, K: 0, Pos: -1})
			curK = append(curK, 0)
			curPos = append(curPos, -1)
			curValid = append(curValid, false)
		}

		blocks = append(blocks, curPairs)
		kOfSlot = append(kOfSlot, curK)
		posInDiag = append(posInDiag, curPos)
		valid = append(valid, curValid)
	}

	return
}

// [a0, a1, ..., a_{m-1}] 암호화
func encryptString(values []int, cc Context) (*rlwe.Ciphertext, error) {
	slots := cc.params.MaxSlots()
	buf := make([]complex128, slots)

	for i := 0; i < len(values); i++ {
		buf[i] = ByteToRoot(values[i])
	}

	pt := ckks.NewPlaintext(*cc.params, cc.params.MaxLevel())
	if err := cc.encoder.Encode(buf, pt); err != nil {
		return nil, err
	}
	return cc.encryptor.EncryptNew(pt)
}

func encryptPacked(values []complex128, cc Context) (*rlwe.Ciphertext, error) {
	pt := ckks.NewPlaintext(*cc.params, cc.params.MaxLevel())
	if err := cc.encoder.Encode(values, pt); err != nil {
		return nil, err
	}
	return cc.encryptor.EncryptNew(pt)
}

func makeMask(slot, slots, level int, cc Context) (*rlwe.Plaintext, error) {
	if slot < 0 || slot >= slots {
		return nil, fmt.Errorf("invalid mask slot=%d for slots=%d", slot, slots)
	}

	mask := make([]complex128, cc.params.MaxSlots())
	mask[slot] = 1

	pt := ckks.NewPlaintext(*cc.params, level)
	if err := cc.encoder.Encode(mask, pt); err != nil {
		return nil, err
	}
	return pt, nil
}

func encryptZero(cc Context) (*rlwe.Ciphertext, error) {
	zeros := make([]complex128, cc.params.MaxSlots())
	return encryptPacked(zeros, cc)
}


// 원본 문자열 암호문에서 필요한 index만 뽑아 현재 block slot 배치로 재구성
// srcIndices[dst] = 원본 문자열에서 dst slot으로 가져올 index
func BlockFromEncrypted(
	ctOrig *rlwe.Ciphertext,
	srcIndices []int,
	slots int,
	cc Context,
) (*rlwe.Ciphertext, error) {

	var acc *rlwe.Ciphertext

	for dst, src := range srcIndices {
		if src < 0 {
			continue
		}

		rot := src - dst
		ctRot, err := cc.eval.RotateNew(ctOrig, rot)
		if err != nil {
			return nil, fmt.Errorf("rotate dst=%d src=%d: %w", dst, src, err)
		}

		maskPt, err := makeMask(dst, slots, ctRot.Level(), cc)
		if err != nil {
			return nil, fmt.Errorf("make mask dst=%d src=%d: %w", dst, src, err)
		}

		ctMasked, err := cc.eval.MulNew(ctRot, maskPt)
		if err != nil {
			return nil, fmt.Errorf("mul mask dst=%d src=%d: %w", dst, src, err)
		}

		if err := cc.eval.Rescale(ctMasked, ctMasked); err != nil {
			return nil, fmt.Errorf("rescale dst=%d src=%d: %w", dst, src, err)
		}

		if acc == nil {
			acc = ctMasked
		} else {
			if err := cc.eval.Add(acc, ctMasked, acc); err != nil {
				return nil, fmt.Errorf("add dst=%d src=%d: %w", dst, src, err)
			}
		}
	}

	if acc == nil {
		return encryptZero(cc)
	}

	return acc, nil
}

func ABlock(
	ctAOrig *rlwe.Ciphertext,
	block []PairIndex,
	slots int,
	cc Context,
) (*rlwe.Ciphertext, error) {
	src := make([]int, len(block))
	for i, p := range block {
		src[i] = p.IA
	}
	return BlockFromEncrypted(ctAOrig, src, slots, cc)
}

func BBlock(
	ctBOrig *rlwe.Ciphertext,
	block []PairIndex,
	slots int,
	cc Context,
) (*rlwe.Ciphertext, error) {
	src := make([]int, len(block))
	for i, p := range block {
		src[i] = p.JB
	}
	return BlockFromEncrypted(ctBOrig, src, slots, cc)
}

func PrecompDeltas(a, b []int, cc Context) (*DeltaTable, error) {
	if len(a) != len(b) {
		return nil, fmt.Errorf("different length")
	}

	m := len(a)
	slots := m * m
	if slots > cc.params.MaxSlots() {
		slots = cc.params.MaxSlots()
	}

	pairs, err := buildIndexes(m)
	if err != nil {
		return nil, err
	}

	blocks, kOfSlot, posInDiag, valid, err := packBlocks(pairs, slots)
	if err != nil {
		return nil, err
	}

	ctAOrig, err := encryptString(a, cc)
	if err != nil {
		return nil, fmt.Errorf("encrypt a: %w", err)
	}

	ctBOrig, err := encryptString(b, cc)
	if err != nil {
		return nil, fmt.Errorf("encrypt b: %w", err)
	}

	deltaBlocks := make([]*rlwe.Ciphertext, len(blocks))
	for blk := range blocks {
		ctAblk, err := ABlock(ctAOrig, blocks[blk], slots, cc)
		if err != nil {
			return nil, fmt.Errorf("A block %d: %w", blk, err)
		}

		ctBblk, err := BBlock(ctBOrig, blocks[blk], slots, cc)
		if err != nil {
			return nil, fmt.Errorf("B block %d: %w", blk, err)
		}

		// ===== debug: block 0의 gather 결과 확인 =====
		if blk == 0 {
			aVals := make([]complex128, cc.params.MaxSlots())
			bVals := make([]complex128, cc.params.MaxSlots())

			if err := cc.encoder.Decode(cc.decryptor.DecryptNew(ctAblk), aVals); err != nil {
				return nil, fmt.Errorf("decode A block %d: %w", blk, err)
			}
			if err := cc.encoder.Decode(cc.decryptor.DecryptNew(ctBblk), bVals); err != nil {
				return nil, fmt.Errorf("decode B block %d: %w", blk, err)
			}

			fmt.Println("=== gathered block 0 before EvalEquality ===")
			for s := 0; s < slots; s++ {
				if !valid[blk][s] {
					continue
				}
				fmt.Printf("slot=%d k=%d pos=%d A=%0.4f%+0.4fi B=%0.4f%+0.4fi\n",
					s,
					kOfSlot[blk][s],
					posInDiag[blk][s],
					real(aVals[s]), imag(aVals[s]),
					real(bVals[s]), imag(bVals[s]),
				)
			}
		}
		// ===== debug 끝 =====


		ctEq := EvalEquality(ctBblk, ctBblk, cc)
		deltaBlocks[blk] = ctEq
	}

	return &DeltaTable{
		DeltaBlocks: deltaBlocks,
		KOfSlot:     kOfSlot,
		PosInDiag:   posInDiag,
		Valid:       valid,
		M:           m,
		Slots:       slots,
		PairCount:   m * m,
	}, nil
}


