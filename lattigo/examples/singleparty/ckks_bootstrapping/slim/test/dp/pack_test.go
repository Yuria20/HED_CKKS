package main

import (
	he "banded_dp/he"
	util "banded_dp/util"
	"flag"
	"fmt"
	"testing"
)

func Test_pack(t *testing.T) {
	flag.Parse()

	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 13

	var cc he.Context
	fmt.Print("Generate Context... ")
	cc.GenContext(LogN, 10, 4, 1)
	fmt.Println("Done!")

	str_len := 8
	char_size := 127

	fmt.Println("Packing Test... ")
	var pack_ct he.StrCipher
	pack_ct.Length = str_len

	strArr_want := util.SampleStringIntMany(cc.Params.MaxSlots(), str_len, char_size)
	util.PrintString(strArr_want, "uint8")

	pack_ct.StrEncrypt(cc, char_size, cc.Exp_lv, strArr_want)

	util.StrDecryptAndPrint(cc, pack_ct, "uint8")

	strArr_test := pack_ct.StrDecrypt(cc)

	equality_test(strArr_want, strArr_test)
}

func equality_test(want, test [][]int) {

	if len(want) != len(test) {
		fmt.Println("error!")
		return
	}

	if len(want[0]) != len(test[0]) {
		fmt.Println("error!")
		return
	}

	N := len(want)
	m := len(want[0])

	for i := 0; i < N; i++ {
		for j := 0; j < m; j++ {
			if want[i][j] != test[i][j] {
				fmt.Printf("want[%d][%d] = %d and test[%d][%d] = %d\n", i, j, want[i][j], i, j, test[i][j])
				return
			}
		}
	}

	fmt.Println("pass!")
}
