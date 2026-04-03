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

// 전개:
// A: a1 | a1,a2 | a1,a2,a3 | ... | a_{m-1},a_m | a_m
// B: b1 | b2,b1 | b3,b2,b1 | ... | b_m,b_{m-1} | b_m
func buildPaperStyleFlattened(a, b []int) (
	flatA []complex128,
	flatB []complex128,
	kList []int,
	posList []int,
	err error,
) {
	if len(a) != len(b) {
		return nil, nil, nil, nil, fmt.Errorf("a and b must have same length")
	}
	m := len(a)
	if m == 0 {
		return nil, nil, nil, nil, fmt.Errorf("empty input")
	}

	// 위쪽 삼각형: k = 2 .. m+1
	for k := 2; k <= m+1; k++ {
		diagLen := k - 1
		for t := 0; t < diagLen; t++ {
			// A: a1, a2, ..., a_{diagLen}
			i := t + 1
			// B: b_{diagLen}, b_{diagLen-1}, ..., b1
			j := diagLen - t

			flatA = append(flatA, ByteToRoot(a[i-1]))
			flatB = append(flatB, ByteToRoot(b[j-1]))
			kList = append(kList, k)
			posList = append(posList, t)
		}
	}

	// 아래쪽 삼각형: k = m+2 .. 2m
	for k := m + 2; k <= 2*m; k++ {
		diagLen := 2*m - k + 1
		start := m - diagLen + 1 // A 시작 인덱스(1-based)

		for t := 0; t < diagLen; t++ {
			// A: a_start, a_{start+1}, ..., a_m
			i := start + t
			// B: b_m, b_{m-1}, ..., b_start
			j := m - t

			flatA = append(flatA, ByteToRoot(a[i-1]))
			flatB = append(flatB, ByteToRoot(b[j-1]))
			kList = append(kList, k)
			posList = append(posList, t)
		}
	}

	return
}

func packBlocks(
	flatA, flatB []complex128,
	kList, posList []int,
	slots int,
) (
	aBlocks [][]complex128,
	bBlocks [][]complex128,
	kOfSlot [][]int,
	posInDiag [][]int,
	valid [][]bool,
	err error,
) {
	if len(flatA) != len(flatB) || len(flatA) != len(kList) || len(flatA) != len(posList) {
		return nil, nil, nil, nil, nil, fmt.Errorf("inconsistent lengths")
	}

	n := len(flatA)
	for start := 0; start < n; start += slots {
		end := start + slots
		if end > n {
			end = n
		}

		curA := append([]complex128{}, flatA[start:end]...)
		curB := append([]complex128{}, flatB[start:end]...)
		curK := append([]int{}, kList[start:end]...)
		curPos := append([]int{}, posList[start:end]...)
		curValid := make([]bool, end-start)
		for i := range curValid {
			curValid[i] = true
		}

		for len(curA) < slots {
			curA = append(curA, complex(1, 0))
			curB = append(curB, complex(1, 0))
			curK = append(curK, 0)
			curPos = append(curPos, -1)
			curValid = append(curValid, false)
		}

		aBlocks = append(aBlocks, curA)
		bBlocks = append(bBlocks, curB)
		kOfSlot = append(kOfSlot, curK)
		posInDiag = append(posInDiag, curPos)
		valid = append(valid, curValid)
	}

	return
}

func encryptPacked(values []complex128, cc Context) (*rlwe.Ciphertext, error) {
	pt := ckks.NewPlaintext(*cc.params, cc.max_lv)
	if err := cc.encoder.Encode(values, pt); err != nil {
		return nil, err
	}
	return cc.encryptor.EncryptNew(pt)
}

func PrecompDeltas(a, b []int, cc Context) (*DeltaTable, error) {
	if len(a) != len(b) {
		return nil, fmt.Errorf("different length")
	}

	m := len(a)
	slots := cc.params.MaxSlots()

	flatA, flatB, kList, posList, err := buildPaperStyleFlattened(a, b)
	if err != nil {
		return nil, err
	}

	aBlocks, bBlocks, kOfSlot, posInDiag, valid, err :=
		packBlocks(flatA, flatB, kList, posList, slots)
	if err != nil {
		return nil, err
	}

	deltaBlocks := make([]*rlwe.Ciphertext, len(aBlocks))
	for blk := range aBlocks {
		ctA, err := encryptPacked(aBlocks[blk], cc)
		if err != nil {
			return nil, fmt.Errorf("encrypt ctA block %d: %w", blk, err)
		}
		ctB, err := encryptPacked(bBlocks[blk], cc)
		if err != nil {
			return nil, fmt.Errorf("encrypt ctB block %d: %w", blk, err)
		}

		ctEq := EvalEquality(ctA, ctB, cc)
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


