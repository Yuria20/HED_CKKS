package main

import (
	"math"
	"math/cmplx"
)

// Description
// Given input n points that x_k = exp(2pi * k/t) (0<= k < n) and output points that y_k = LUT(k) where LUT is look up table,
// this function outputs coefficients of f(x) such that f(k)=y_k for each 0<= k <n.
// the hermiteOrder means that level of analytic stability.
// Note : the degree of f is k*n -1
func HermiteInterpolation(n int, hermiteOrder int, LUT func(int) int) []complex128 {
	k := hermiteOrder
	// n개의 점, 각 점마다 (k+1)개의 제약 조건 (0차~k차 도함수)
	// 총 제약 조건의 개수(N)이자, 계산될 다항식 계수의 개수
	numConstraints := n * (k + 1)

	// A는 (N x N) 행렬, b는 (N x 1) 벡터
	A := make([][]complex128, numConstraints)
	b := make([]complex128, numConstraints)

	// 기본 n차 단위근 z = e^(2*pi*i / n)
	z := cmplx.Exp(2 * math.Pi * 1i / complex(float64(n), 0))

	row := 0 // A 행렬과 b 벡터의 현재 행 인덱스
	for i := 0; i < n; i++ {
		// i번째 보간점 z_i = z^i
		zi := cmplx.Pow(z, complex(float64(i), 0))

		// 0차 도함수부터 k차 도함수까지 제약 조건 설정
		for m := 0; m <= k; m++ {
			A[row] = make([]complex128, numConstraints)

			// P(x) = c_0 + c_1*x + ... + c_{N-1}*x^{N-1}
			// m차 도함수는 P^(m)(x) = sum_{j=m}^{N-1} c_j * (j!/(j-m)!) * x^(j-m)
			for j := 0; j < numConstraints; j++ {
				if j < m {
					A[row][j] = 0
				} else {
					// 계수 부분: j! / (j-m)!
					prod := complex(1.0, 0)
					for p := 0; p < m; p++ {
						prod *= complex(float64(j-p), 0)
					}
					// 변수 부분: (z_i)^(j-m)
					zi_pow := cmplx.Pow(zi, complex(float64(j-m), 0))
					A[row][j] = prod * zi_pow
				}
			}

			// b 벡터 (제약 조건의 결과값) 설정
			if m == 0 {
				b[row] = complex(float64(LUT(i)), 0)
			} else {
				b[row] = 0
			}

			row++ // 다음 행으로 이동
		}
	}

	// 선형 시스템 A * list = b 를 풀어 계수(list)를 구합니다.
	// list의 길이는 numConstraints 입니다.
	list := solveLinearSystem(A, b)

	// 계산된 계수 슬라이스를 그대로 반환합니다.
	return list
}

func solveLinearSystem(A [][]complex128, b []complex128) []complex128 {
	n := len(b)
	for i := 0; i < n; i++ {
		maxRow := i
		for k := i + 1; k < n; k++ {
			if cmplx.Abs(A[k][i]) > cmplx.Abs(A[maxRow][i]) {
				maxRow = k
			}
		}
		A[i], A[maxRow] = A[maxRow], A[i]
		b[i], b[maxRow] = b[maxRow], b[i]

		pivot := A[i][i]
		for j := i; j < n; j++ {
			A[i][j] /= pivot
		}
		b[i] /= pivot

		for k := i + 1; k < n; k++ {
			factor := A[k][i]
			for j := i; j < n; j++ {
				A[k][j] -= factor * A[i][j]
			}
			b[k] -= factor * b[i]
		}
	}

	x := make([]complex128, n)
	for i := n - 1; i >= 0; i-- {
		x[i] = b[i]
		for j := i + 1; j < n; j++ {
			x[i] -= A[i][j] * x[j]
		}
	}

	return x
}
