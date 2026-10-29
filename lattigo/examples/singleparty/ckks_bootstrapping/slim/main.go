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
	"banded_dp/dp"
	"banded_dp/dp/dp_ideal"
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

	for i := 8; i >= 3; i-- {

		cc.M = (1 << i)
		fmt.Printf("===================== str_len = %d ===============================\n", cc.M)
		str_len := cc.M
		// [cpx2int] + (eval depth) * (# eval) + (eq depth)
		fmt.Print("Generate Context... \n")
		cc.GenContext(LogN, 1+2*4+10, 0, 0, str_len)
		fmt.Println("Done!")

		cc.C_size = 128
		cc.Batch_size = cc.Slot_size / cc.M

		for iter := 0; iter < 1; iter++ {
			var ct_A he.StrCipher
			ct_A.Length = str_len
			str_A := util.SampleStringIntMany(cc.Batch_size, str_len, cc.C_size)
			//str_A[0] = []int{109, 98, 77, 44, 77, 25, 58, 9}
			util.PrintString(str_A, "uint8")
			ct_A.StrEncrypt(cc, cc.Exp_lv, str_A)

			// Many 해당하는 string을 미리 reverse 시킨 후 암호화 합니다.
			var ct_B he.StrCipher
			ct_B.Length = str_len
			str_B := util.SampleStringIntMany(cc.Batch_size, str_len, cc.C_size)
			//str_B[0] = []int{91, 45, 103, 78, 72, 77, 25, 0}
			rev_str_B := util.ReverseString(str_B)
			//fmt.Println(rev_str_B[0])
			util.PrintString(str_B, "uint8")
			ct_B.StrEncrypt(cc, cc.Exp_lv, rev_str_B)

			var st dp.DPState
			dp_start := time.Now()
			st.GenDPstate(cc, ct_A, ct_B)
			st.ComputeDpTable(cc)
			dp_time := time.Since(dp_start).Seconds()
			fmt.Println("Latency : ", dp_time-st.Pre_time)
			fmt.Println("AMortized : ", (dp_time-st.Pre_time)/float64(cc.Batch_size))
			fmt.Println("Precomputation time : ", st.Pre_time)
			fmt.Println("Precomputation time (amortized) : ", st.Pre_time/float64(cc.Batch_size))
			fmt.Println("PCD time : ", st.Eq_time)
			fmt.Println("PCD time (amortized) : ", st.Eq_time/float64(cc.Batch_size))
			fmt.Println("bootstrapping time : ", st.Btp_time)
			fmt.Println("bootstrapping time (amortized) : ", st.Btp_time/float64(cc.Batch_size))
			fmt.Println("Others time : ", dp_time-st.Pre_time-st.Eq_time-st.Btp_time)
			fmt.Println("Others time (amortized) : ", (dp_time-st.Pre_time-st.Eq_time-st.Btp_time)/float64(cc.Batch_size))
			res_ideal := make([]int, cc.Batch_size)
			for k := 0; k < cc.Batch_size; k++ {
				res_ideal[k] = dp_ideal.LevenshteinDistance(str_A[k], str_B[k])
			}

			res_sim := st.Res_sim
			fmt.Println(res_ideal)
			fmt.Println(res_sim)

			for k := 0; k < cc.Batch_size; k++ {
				if res_ideal[k] != res_sim[k] {
					fmt.Println("k = ", k)
					fmt.Println(str_A[k])
					fmt.Println(str_B[k])

					break
				}
			}
			fmt.Println("Pass!")
		}
		fmt.Println()
		fmt.Println()
	}
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
