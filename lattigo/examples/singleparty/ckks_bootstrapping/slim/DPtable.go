package main

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

const AlphabetSize = 256

type DeltaTable struct {
	DeltaBlocks []*rlwe.Ciphertext // anti-diagonal별 eq 결과
	KOfSlot     [][]int            // 각 block의 slot이 속한 k
	PosInDiag   [][]int            // 각 block의 slot이 anti-diagonal 내부에서 몇 번째인지
	Valid       [][]bool           // 유효 slot 여부
	M           int                // 문자열 길이
	Slots       int                // 전체 사용 가능 slot 수
	PairCount   int                // 전체 비교쌍 개수 = N*m*m
}

// 복소원 매핑
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

func encryptPacked(values []complex128, cc Context) (*rlwe.Ciphertext, error) {
	pt := ckks.NewPlaintext(*cc.params, cc.max_lv)
	if err := cc.encoder.Encode(values, pt); err != nil {
		return nil, err
	}
	return cc.encryptor.EncryptNew(pt)
}

func encryptZero(cc Context) (*rlwe.Ciphertext, error) {
	zeros := make([]complex128, cc.params.MaxSlots())
	return encryptPacked(zeros, cc)
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

func chunkSize(m, n, maxSlots int) (int, error) {
	inputSlots := n * m
	if inputSlots == 0 {
		return 0, fmt.Errorf("empty input")
	}

	chunk := maxSlots / inputSlots
	if chunk < 1 {
		return 0, fmt.Errorf("invalid chunk")
	}
	return chunk, nil
}

// query 하나를 N개 block에 복제해서 sparse하게 pack
func packQSparse(q []int, n int, chunk int, cc Context) (*rlwe.Ciphertext, error) {
	m := len(q)
	buf := make([]complex128, cc.params.MaxSlots())

	for blk := 0; blk < n; blk++ {
		base := blk * m * chunk
		for i := 0; i < m; i++ {
			pos := base + i*chunk
			if pos >= len(buf) {
				return nil, fmt.Errorf("packQSparse overflow")
			}
			buf[pos] = ByteToRoot(q[i])
		}
	}

	return encryptPacked(buf, cc)
}

// target N개를 sparse하게 pack
func packTSparse(ts [][]int, chunk int, cc Context) (*rlwe.Ciphertext, error) {
	if len(ts) == 0 {
		return nil, fmt.Errorf("empty targets")
	}

	m := len(ts[0])
	buf := make([]complex128, cc.params.MaxSlots())

	for blk := 0; blk < len(ts); blk++ {
		if len(ts[blk]) != m {
			return nil, fmt.Errorf("different target length at %d", blk)
		}

		base := blk * m * chunk
		for j := 0; j < m; j++ {
			pos := base + j*chunk
			if pos >= len(buf) {
				return nil, fmt.Errorf("packTSparse overflow")
			}
			buf[pos] = ByteToRoot(ts[blk][j])
		}
	}

	return encryptPacked(buf, cc)
}

// 각 target block 내부에서 문자 순서를 뒤집음
// GBFV의 ReverseCharacters 대응
func revSparse(ct *rlwe.Ciphertext, m, n, chunk int, cc Context) (*rlwe.Ciphertext, error) {
	acc, err := encryptZero(cc)
	if err != nil {
		return nil, err
	}

	for blk := 0; blk < n; blk++ {
		base := blk * m * chunk

		for i := 0; i < m; i++ {
			src := base + i*chunk
			dst := base + (m-1-i)*chunk
			rot := src - dst

			ctRot, err := cc.eval.RotateNew(ct, rot)
			if err != nil {
				return nil, fmt.Errorf("rev rotate blk=%d i=%d: %w", blk, i, err)
			}

			maskPt, err := makeMask(dst, cc.params.MaxSlots(), ctRot.Level(), cc)
			if err != nil {
				return nil, fmt.Errorf("rev mask blk=%d i=%d: %w", blk, i, err)
			}

			tmp, err := cc.eval.MulNew(ctRot, maskPt)
			if err != nil {
				return nil, fmt.Errorf("rev mul blk=%d i=%d: %w", blk, i, err)
			}

			if err := cc.eval.Rescale(tmp, tmp); err != nil {
				return nil, fmt.Errorf("rev rescale blk=%d i=%d: %w", blk, i, err)
			}

			if err := cc.eval.Add(acc, tmp, acc); err != nil {
				return nil, fmt.Errorf("rev add blk=%d i=%d: %w", blk, i, err)
			}
		}
	}

	return acc, nil
}

// k번째 anti-diagonal에 필요한 문자들을 한 block으로 확장
// reverse=false : query 쪽 인덱스 사용
// reverse=true  : reversed target 쪽 인덱스 사용
func expandDiag(ct *rlwe.Ciphertext, m, n, chunk, k int, reverse bool, cc Context) (*rlwe.Ciphertext, []int, []int, error) {
	acc, err := encryptZero(cc)
	if err != nil {
		return nil, nil, nil, err
	}

	kOfSlot := []int{}
	posOfSlot := []int{}
	dst := 0

	for blk := 0; blk < n; blk++ {
		var pairs [][2]int

		if k <= m+1 {
			diagLen := k - 1
			for t := 0; t < diagLen; t++ {
				i := t
				j := diagLen - 1 - t
				pairs = append(pairs, [2]int{i, j})
			}
		} else {
			diagLen := 2*m - k + 1
			start := m - diagLen
			for t := 0; t < diagLen; t++ {
				i := start + t
				j := m - 1 - t
				pairs = append(pairs, [2]int{i, j})
			}
		}

		for t, p := range pairs {
			idx := p[0]
			if reverse {
				idx = p[1]
			}

			src := blk*m*chunk + idx*chunk
			rot := src - dst

			ctRot, err := cc.eval.RotateNew(ct, rot)
			if err != nil {
				return nil, nil, nil, fmt.Errorf("expand rotate blk=%d k=%d pos=%d: %w", blk, k, t, err)
			}

			maskPt, err := makeMask(dst, cc.params.MaxSlots(), ctRot.Level(), cc)
			if err != nil {
				return nil, nil, nil, fmt.Errorf("expand mask blk=%d k=%d pos=%d: %w", blk, k, t, err)
			}

			tmp, err := cc.eval.MulNew(ctRot, maskPt)
			if err != nil {
				return nil, nil, nil, fmt.Errorf("expand mul blk=%d k=%d pos=%d: %w", blk, k, t, err)
			}

			if err := cc.eval.Rescale(tmp, tmp); err != nil {
				return nil, nil, nil, fmt.Errorf("expand rescale blk=%d k=%d pos=%d: %w", blk, k, t, err)
			}

			if err := cc.eval.Add(acc, tmp, acc); err != nil {
				return nil, nil, nil, fmt.Errorf("expand add blk=%d k=%d pos=%d: %w", blk, k, t, err)
			}

			kOfSlot = append(kOfSlot, k)
			posOfSlot = append(posOfSlot, t)
			dst++
		}
	}

	return acc, kOfSlot, posOfSlot, nil
}

// query 1개 vs target N개를 sparse input + expand 방식으로 eq table 생성
func Delta1N(q []int, ts [][]int, cc Context) (*DeltaTable, error) {
	if len(ts) == 0 {
		return nil, fmt.Errorf("empty targets")
	}

	m := len(q)
	n := len(ts)

	for i := 0; i < n; i++ {
		if len(ts[i]) != m {
			return nil, fmt.Errorf("different length at target %d", i)
		}
	}

	chunk, err := chunkSize(m, n, cc.params.MaxSlots())
	if err != nil {
		return nil, err
	}

	ctQ, err := packQSparse(q, n, chunk, cc)
	if err != nil {
		return nil, fmt.Errorf("packQSparse: %w", err)
	}

	ctT0, err := packTSparse(ts, chunk, cc)
	if err != nil {
		return nil, fmt.Errorf("packTSparse: %w", err)
	}

	ctT, err := revSparse(ctT0, m, n, chunk, cc)
	if err != nil {
		return nil, fmt.Errorf("revSparse: %w", err)
	}

	var deltaBlocks []*rlwe.Ciphertext
	var kOfSlot [][]int
	var posInDiag [][]int
	var valid [][]bool

	for k := 2; k <= 2*m; k++ {
		ctA, ks, ps, err := expandDiag(ctQ, m, n, chunk, k, false, cc)
		if err != nil {
			return nil, fmt.Errorf("expandDiag query k=%d: %w", k, err)
		}

		ctB, _, _, err := expandDiag(ctT, m, n, chunk, k, true, cc)
		if err != nil {
			return nil, fmt.Errorf("expandDiag target k=%d: %w", k, err)
		}

		ctEq := EvalEquality(ctA, ctB, cc)

		deltaBlocks = append(deltaBlocks, ctEq)
		kOfSlot = append(kOfSlot, ks)
		posInDiag = append(posInDiag, ps)

		v := make([]bool, len(ks))
		for i := range v {
			v[i] = true
		}
		valid = append(valid, v)
	}

	return &DeltaTable{
		DeltaBlocks: deltaBlocks,
		KOfSlot:     kOfSlot,
		PosInDiag:   posInDiag,
		Valid:       valid,
		M:           m,
		Slots:       cc.params.MaxSlots(),
		PairCount:   n * m * m,
	}, nil
}