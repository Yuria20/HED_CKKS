package util

import (
	dp "banded_dp/dp/dp_sim"
	he "banded_dp/he"
	"fmt"
)

// String 묶음을 한번에 출력함.
// tag = "uint8" 이면 숫자로 출력.
// tag = "str" 이면 string으로 출력.
func PrintString(strpt [][]int, tag string) {

	// strpt : N x m
	m := len(strpt[0])
	N := len(strpt)
	for i := 0; i < 2; i++ {
		fmt.Printf("str[%d] : ", i)
		for j := 0; j < m; j++ {
			if tag == "uint8" {
				fmt.Print(strpt[i][j])
				fmt.Print(" ")
			} else if tag == "str" {
				fmt.Print(string(strpt[i][j]))
			}
		}
		fmt.Println()
	}
	fmt.Println("......")
	for i := N - 1; i < N; i++ {
		fmt.Printf("str[%d] : ", i)
		for j := 0; j < m; j++ {
			if tag == "uint8" {
				fmt.Print(strpt[i][j])
				fmt.Print(" ")
			} else if tag == "str" {
				fmt.Print(string(strpt[i][j]))
			}
		}
		fmt.Println()
	}
}

// StrCipher 객체를 복호화하고 출력함.
// tag = "uint8" 이면 숫자로 출력.
// tag = "str" 이면 string으로 출력.
func StrDecryptAndPrint(cc he.Context, strct dp.StrCipher, tag string) {
	PrintString(strct.StrDecrypt(cc), tag)
}
