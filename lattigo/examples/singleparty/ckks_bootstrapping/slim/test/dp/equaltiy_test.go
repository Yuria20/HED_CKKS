// Package main implements an example showcasing slim for bootstrapping for fixed-point approximate
// arithmetic over the reals/complexes numbers.
// This re-ordering of the bootstrapping steps was first proposed for the BFV/BGV schemes by Chen and Han
// in Homomorphic Lower Digits Removal and Improved FHE Bootstrapping (https://eprint.iacr.org/2018/067).
// It was also used by Kim and Guyot in Optimized Privacy-Preserving CNN Inference With Fully Homomorphic
// Encryption (https://ieeexplore.ieee.org/document/10089847) to efficiently perform the convolution in
// the coefficient domain.
//
// This example assumes that the user is already familiar with the bootstrapping and its different steps.
// See the basic example `lattigo/single_party/applications/reals_bootstrapping/basics` for an introduction into the
// bootstrapping.
//
// The usual order of the bootstrapping operations is:
//
// 0) User defined circuit in the slots domain
// 1) ScaleDown: Scale the ciphertext to q0/|m|
// 2) ModUp: Raise modulus from q0 to qL
// 3) CoeffsToSlots: Homomorphic encoding
// 4) EvalMod: Homomorphic modular reduction
// 5) SlotsToCoeffs (and go back to 0): Homomorphic Decoding
//
// This example instantiates a custom order of the circuit evaluating:
//
// 0) User defined circuit in the slots domain
// 1) SlotsToCoeffs: Homomorphic Decoding
// 2) User defined circuit in the coeffs domain
// 3) ScaleDown: Scale the ciphertext to q0/|m|
// 4) ModUp: Raise modulus from q0 to qL
// 5) CoeffsToSlots: Homomorphic encoding
// 6) EvalMod (and to back to 0): Homomorphic modular reduction
//
// Use the flag -short to run the examples fast but with insecure parameters.
package main

import (
	"banded_dp/debug"
	dp "banded_dp/dp"
	he "banded_dp/he"
	util "banded_dp/util"
	"flag"
	"fmt"
	"math"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func Test_equality(t *testing.T) {

	flag.Parse()

	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 16

	if *flagShort {
		LogN -= 3
	}

	var cc he.Context
	fmt.Print("Generate Context... ")
	cc.GenContext(LogN, 5, 4, 1)
	fmt.Println("Done!")

	str_len := 8
	cc.C_size = 128

	fmt.Println("Packing Test... ")
	var pack_ct he.StrCipher
	pack_ct.Length = str_len

	strArr := util.SampleStringIntMany(cc.Params.MaxSlots(), str_len, cc.C_size)
	util.PrintString(strArr, "uint8")

	pack_ct.StrEncrypt(cc, cc.Exp_lv, strArr)
	fmt.Println(pack_ct.Cipher[0].Level())
	eq_time := time.Now()
	ct, _ := dp.EvalEQ_exp(cc, 0, pack_ct.Cipher[0], pack_ct.Cipher[1])
	fmt.Println(ct.Level())
	fmt.Println("eq time : ", time.Since(eq_time))
	debug.PrintNoise(cc.Params, ct, he.DEC(cc, ct), cc.Decryptor, cc.Encoder)

	util.StrDecryptAndPrint(cc, pack_ct, "uint8")
}

func printDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %10.14f %10.14f %10.14f %10.14f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %10.14f %10.14f %10.14f %10.14f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])

	precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	fmt.Println(precStats.String())
	fmt.Println()

	return
}
