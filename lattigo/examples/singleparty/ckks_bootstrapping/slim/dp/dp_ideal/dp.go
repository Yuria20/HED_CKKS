package dp_ideal

import (
	he "banded_dp/he"
	"banded_dp/util"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

func PreComputeDeltas(cc he.Context, hermite int, ctxt1 *rlwe.Ciphertext, ctxt2 *rlwe.Ciphertext) []*rlwe.Ciphertext {
	sz := cc.C_size
	str_len := cc.M
	b_sz := cc.Batch_size

	in_lv := min(ctxt1.Level(), ctxt2.Level())
	out_lv := in_lv - (int(math.Ceil(math.Log2(float64(sz)))) + 1)

	values1 := he.InverseInterLeave(he.DEC(cc, ctxt1), cc.M, cc.Slot_size)
	values2 := he.InverseInterLeave(he.DEC(cc, ctxt2), cc.M, cc.Slot_size)
	for i := 0; i < cc.Slot_size; i++ {
		values1[i] = complex(float64(he.RootToByte(values1[i], sz)), 0)
		values2[i] = complex(float64(he.RootToByte(values2[i], sz)), 0)
	}

	A := make([]int, str_len)
	for i := 0; i < len(A); i++ {
		A[i] = int(math.Round(real(values1[i])))
	}
	B_list := make([][]int, b_sz)
	for i := 0; i < len(B_list); i++ {
		B_list[i] = make([]int, str_len)
		for j := 0; j < str_len; j++ {
			B_list[i][j] = int(math.Round(real(values2[i*str_len+j])))
		}
	}

	GenTable := func(str1 []int, str2 []int) [][]complex128 {
		E := make([][]complex128, str_len)
		for i := 0; i < str_len; i++ {
			E[i] = make([]complex128, str_len)
			for j := 0; j < str_len; j++ {
				diff := str2[i] - str1[j]
				if diff == 0 {
					E[i][j] = -1
				} else {
					E[i][j] = 0
				}

			}
		}

		return ToCyclicAntiDiag2D(E)
	}

	E_list := make([][][]complex128, b_sz)
	for i := 0; i < b_sz; i++ {
		E_list[i] = GenTable(A, util.ReverseString(B_list)[i])

	}
	fmt.Println(A, B_list[0])
	fmt.Println(E_list[0])

	rtn_list := make([]*rlwe.Ciphertext, str_len)
	for i := 0; i < str_len; i++ {
		values_want := make([]complex128, cc.Params.MaxSlots())
		for j := 0; j < b_sz; j++ {
			for k := 0; k < str_len; k++ {
				values_want[j*str_len+k] = E_list[j][i][k]
			}
		}
		rtn_list[i], _ = he.ENC(cc, he.InterLeave(values_want, cc.M, cc.Slot_size), out_lv)
	}

	return rtn_list
}

// ToCyclicAntiDiag2D: m x m 행렬을 cyclic anti-diagonal 기준 2차원 배열로 변환
// result[d][j] = matrix[i][j],  d = (i+j) mod m
func ToCyclicAntiDiag2D[T any](matrix [][]T) [][]T {
	m := len(matrix)
	result := make([][]T, m)
	for d := range result {
		result[d] = make([]T, m)
	}

	for i := 0; i < m; i++ {
		for j := 0; j < m; j++ {
			d := (i + j) % m
			result[d][j] = matrix[i][j]
		}
	}
	return result
}

// FromCyclicAntiDiag2D: cyclic anti-diagonal 2D 배열 -> 원본 m x m 행렬 복원
// matrix[i][j] = result[d][j],  i = (d-j) mod m
func FromCyclicAntiDiag2D[T any](result [][]T) [][]T {
	m := len(result)
	matrix := make([][]T, m)
	for i := range matrix {
		matrix[i] = make([]T, m)
	}

	for d := 0; d < m; d++ {
		for j := 0; j < m; j++ {
			i := ((d-j)%m + m) % m // Go의 %는 음수를 반환할 수 있어 보정
			matrix[i][j] = result[d][j]
		}
	}
	return matrix
}

// levenshteinDistance 는 길이가 같은 두 int 슬라이스 사이의
// 레벤슈타인 편집 거리(삽입/삭제/치환 최소 횟수)를 계산합니다.
func LevenshteinDistance(a, b []int) int {
	m := len(a) // len(a) == len(b) 라고 가정

	prev := make([]int, m+1)
	curr := make([]int, m+1)

	for j := 0; j <= m; j++ {
		prev[j] = j
	}

	for i := 1; i <= m; i++ {
		curr[0] = i
		for j := 1; j <= m; j++ {
			if a[i-1] == b[j-1] {
				curr[j] = prev[j-1]
			} else {
				curr[j] = 1 + minInt(prev[j], minInt(curr[j-1], prev[j-1]))
			}
		}
		prev, curr = curr, prev
	}

	return prev[m]
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}
