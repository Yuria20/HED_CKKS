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
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func main() {

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
	cc.C_size = 127

	fmt.Println("Packing Test... ")
	var pack_ct he.StrCipher
	pack_ct.Length = str_len

	strArr := util.SampleStringIntMany(cc.Params.MaxSlots(), str_len, cc.C_size)
	util.PrintString(strArr, "uint8")

	pack_ct.StrEncrypt(cc, cc.Exp_lv, strArr)
	fmt.Println(pack_ct.Cipher[0].Level())
	eq_time := time.Now()
	E := dp.GenEqTable(cc, pack_ct, pack_ct)
	fmt.Println("eq time : ", time.Since(eq_time))
	debug.PrintDebug(cc.Params, E.Deltas[0][0], he.DEC(cc, E.Deltas[0][0]), cc.Decryptor, cc.Encoder)

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

// func PrintDelta(dt *DeltaTable, ct *rlwe.Ciphertext, encoder *ckks.Encoder, decryptor *rlwe.Decryptor) error {
// 	slots := ct.Slots()
// 	if !ct.IsBatched {
// 		slots *= 2
// 	}

// 	values := make([]complex128, slots)
// 	if err := encoder.Decode(decryptor.DecryptNew(ct), values); err != nil {
// 		return err
// 	}

// 	for s := 0; s < slots; s++ {
// 		if !dt.Valid[0][s] {
// 			continue
// 		}
// 		fmt.Printf("slot=%d k=%d pos=%d eq≈%.4f\n",
// 			s,
// 			dt.KOfSlot[0][s],
// 			dt.PosInDiag[0][s],
// 			real(values[s]),
// 		)
// 	}
// 	return nil
// }

// func needRotsExpand(maxSlots, m, n int) []int {
// 	chunk := maxSlots / (n * m)
// 	if chunk < 1 {
// 		return nil
// 	}

// 	seen := make(map[int]bool)
// 	rots := make([]int, 0)

// 	// revSparse에서 필요한 회전
// 	for blk := 0; blk < n; blk++ {
// 		base := blk * m * chunk
// 		for i := 0; i < m; i++ {
// 			src := base + i*chunk
// 			dst := base + (m-1-i)*chunk
// 			r := src - dst
// 			if r != 0 && !seen[r] {
// 				seen[r] = true
// 				rots = append(rots, r)
// 			}
// 		}
// 	}

// 	// expandDiag에서 필요한 회전
// 	for k := 2; k <= 2*m; k++ {
// 		dst := 0
// 		diagLen := DiagLen(k, m)

// 		for blk := 0; blk < n; blk++ {
// 			for t := 0; t < diagLen; t++ {
// 				var i, j int
// 				if k <= m+1 {
// 					i = t
// 					j = diagLen - 1 - t
// 				} else {
// 					start := m - diagLen
// 					i = start + t
// 					j = m - 1 - t
// 				}

// 				srcQ := blk*m*chunk + i*chunk
// 				srcT := blk*m*chunk + j*chunk

// 				rq := srcQ - dst
// 				rt := srcT - dst

// 				if rq != 0 && !seen[rq] {
// 					seen[rq] = true
// 					rots = append(rots, rq)
// 				}
// 				if rt != 0 && !seen[rt] {
// 					seen[rt] = true
// 					rots = append(rots, rt)
// 				}
// 				dst++
// 			}
// 		}
// 	}

// 	return rots
// }

// func ijOf(k, pos, m int) (int, int) {
// 	if k <= m+1 {
// 		i := pos
// 		j := (k - 2) - pos
// 		return i, j
// 	}
// 	diagLen := 2*m - k + 1
// 	start := m - diagLen
// 	i := start + pos
// 	j := m - 1 - pos
// 	return i, j
// }
