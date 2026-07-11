package util

import (
	"math/rand"
)

// [0,sz) 안에서 랜덤 정수 sampling
func SampleUint8(sz int) int {
	return rand.Intn(sz)
}

// 길이 m인 string을 sampling 함.
// 각 문자는 [0,sz) 사이에서 sampling 됨.
func SampleStringInts(m int, sz int) []int {
	s := make([]int, m)
	for i := 0; i < m; i++ {
		s[i] = SampleUint8(sz)
		//s[i] = i //SampleUint8(sz)
	}
	return s
}

// 길이 m인 string을 N개 만큼 sampling 함.
// 각 문자는 [0,sz) 사이에서 sampling 됨.
func SampleStringIntMany(N int, m int, sz int) [][]int {
	valuesArr := make([][]int, N)
	for i := 0; i < N; i++ {
		valuesArr[i] = SampleStringInts(m, sz)
	}
	return valuesArr
}

func ReverseString(vals [][]int) [][]int {
	N := len(vals)
	m := len(vals[0])

	valuesArr := make([][]int, N)
	for i := 0; i < N; i++ {
		valuesArr[i] = make([]int, m)
		for j := 1; j <= m; j++ {
			valuesArr[i][j-1] = vals[i][m-j]
		}
	}
	return valuesArr
}

// ???
func RandomSetMatch(query []int, n int, alphabetSize int) [][]int {
	m := len(query)
	targets := make([][]int, n)

	for i := 0; i < n; i++ {
		targets[i] = SampleStringInts(m, alphabetSize)
	}

	if n > 0 {
		match := make([]int, m)
		copy(match, query)
		targets[0] = match
	}

	return targets
}
