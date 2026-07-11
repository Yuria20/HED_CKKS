package dp

import (
	"banded_dp/debug"
	dp_ideal "banded_dp/dp/dp_ideal"
	dp_sim "banded_dp/dp/dp_sim"
	he "banded_dp/he"
	"fmt"
	"math"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type DPState struct {
	CurrAntiDiagIdx int
	E               []*rlwe.Ciphertext
	H               []*rlwe.Ciphertext
	V               []*rlwe.Ciphertext
	dist            *rlwe.Ciphertext

	Res_sim []int
}

func (state *DPState) GenDPstate(cc he.Context, strct1 he.StrCipher, strct2 he.StrCipher) {
	state.CurrAntiDiagIdx = 0

	PCD_start := time.Now()
	state.GenEqTable(cc, strct1, strct2)
	PCD_time := time.Since(PCD_start)
	fmt.Print("PCD time (latency)", float64(PCD_time.Seconds()))
	fmt.Println()
	fmt.Print("PCD time (amortized)", float64(PCD_time.Seconds())/float64(cc.Batch_size))
	fmt.Println()

	// Debug
	for i := 0; i < cc.M; i++ {
		fmt.Println(state.E[i].Level())
		//state.E[i], _ = he.BootstrapIdeal(cc, state.E[i])
	}

	state.H = make([]*rlwe.Ciphertext, cc.M)
	state.V = make([]*rlwe.Ciphertext, cc.M)
}

func (state *DPState) GenEqTable(cc he.Context, strct1 he.StrCipher, strct2 he.StrCipher) {

	// debug option! : 0 (ciphertext), 1 (plaintext), 2 (validation)
	debug_opt := debug.Debug.EQ
	hermite := 0
	ctxt1 := strct1.Cipher
	ctxt2 := strct2.Cipher

	if debug_opt == 0 {
		state.E = dp_sim.PreComputeDeltas(cc, hermite, ctxt1, ctxt2)
	} else if debug_opt == 1 {
		state.E = dp_ideal.PreComputeDeltas(cc, hermite, ctxt1, ctxt2)
	} else {
		state.E = dp_sim.PreComputeDeltas(cc, hermite, ctxt1, ctxt2)
		temp := dp_ideal.PreComputeDeltas(cc, hermite, ctxt1, ctxt2)

		for i := 0; i < cc.M; i++ {
			fmt.Println("index : ", i)
			debug.PrintDebug(cc, cc.Params, state.E[i], he.DEC(cc, temp[i]), cc.Decryptor, cc.Encoder)
		}
		//return ct, err
	}
}

func (state *DPState) ComputeDpTable(cc he.Context) {

	//opt := 0 // once
	opt := 1 // folding

	if opt == 0 {
		// forward => m-1, m th ciphertext
		state.DpOnce(cc)
		return
	} else if opt == 1 {
		state.DpFold(cc)
	}
}

func Reversed[T any](arr []T) []T {
	res := make([]T, len(arr))

	for i := range arr {
		res[len(arr)-1-i] = arr[i]
	}

	return res
}

func (state *DPState) DpFold(cc he.Context) {

	eval := cc.Eval

	// ones // forward(d=0), backwoard(d=1)
	plainList := make([][]*rlwe.Plaintext, 2)
	for d := 0; d < 2; d++ {
		plainList[d] = make([]*rlwe.Plaintext, cc.M)
		for i := 0; i < cc.M; i++ {
			values := make([]complex128, cc.Slot_size)
			for j := 0; j < cc.Batch_size; j++ {
				temp := make([]complex128, cc.M)
				for k := 0; k < cc.M; k++ {
					if k == i {
						temp[k] = complex(0.0, 0.5)
					}
				}

				if d == 1 {
					temp = Reversed(temp)
				}

				for k := 0; k < cc.M; k++ {
					values[j*cc.M+k] = temp[k]
				}

				if j == 0 {
					//fmt.Printf("d=%d, ctIdx=%d\n", d, i)
					//fmt.Println(temp)
				}
			}

			values = he.InterLeave(values, cc.M, cc.Slot_size)
			plainList[d][i] = ckks.NewPlaintext(cc.Params, cc.Max_lv)
			if err := cc.Encoder.Encode(values, plainList[d][i]); err != nil { // values: []complex128 또는 []float64
				panic(err)
			}
		}
	}

	// mask
	maskList := make([][]*rlwe.Plaintext, 2)
	for d := 0; d < 2; d++ {
		maskList[d] = make([]*rlwe.Plaintext, cc.M)
		for i := 0; i < cc.M; i++ {
			values := make([]complex128, cc.Slot_size)
			for j := 0; j < cc.Batch_size; j++ {
				temp := make([]complex128, cc.M)
				for k := 0; k < cc.M; k++ {
					if k <= i {
						if i == 0 {
							temp[k] = complex(0.5, 0.0)
						} else {
							temp[k] = complex(1.0, 0.0)
						}
					}
				}

				if d == 1 {
					temp = Reversed(temp)
				}

				for k := 0; k < cc.M; k++ {
					values[j*cc.M+k] = temp[k]
				}

				if j == 0 {
					//fmt.Printf("d=%d, ctIdx=%d\n", d, i)
					//fmt.Println(temp)
				}
			}

			values = he.InterLeave(values, cc.M, cc.Slot_size)
			maskList[d][i] = ckks.NewPlaintext(cc.Params, cc.Params.MaxLevel())
			if err := cc.Encoder.Encode(values, maskList[d][i]); err != nil { // values: []complex128 또는 []float64
				panic(err)
			}
		}
	}

	// 보정
	corrList := make([][]*rlwe.Plaintext, 2)
	for d := 0; d < 2; d++ {
		corrList[d] = make([]*rlwe.Plaintext, cc.M)
		for i := 0; i < cc.M; i++ {
			values := make([]complex128, cc.Slot_size)
			for j := 0; j < cc.Batch_size; j++ {

				temp := make([]complex128, cc.M)

				for k := 0; k < cc.M; k++ {
					if i < k {
						temp[k] = complex(0.0, 0.5)
					}
				}

				if d == 1 {
					temp = Reversed(temp)
				}

				for k := 0; k < cc.M; k++ {
					values[j*cc.M+k] = temp[k]
				}

				if j == 0 {
					//fmt.Printf("d=%d, ctIdx=%d\n", d, i)
					//fmt.Println(temp)
				}
			}

			values = he.InterLeave(values, cc.M, cc.Slot_size)
			corrList[d][i] = ckks.NewPlaintext(cc.Params, cc.Params.MaxLevel())
			if err := cc.Encoder.Encode(values, corrList[d][i]); err != nil { // values: []complex128 또는 []float64
				panic(err)
			}
		}
	}

	// up init
	H_up, err := eval.MulRelinNew(state.E[0], maskList[0][0])
	if err != nil {
		panic(err)
	}
	if err = eval.Rescale(H_up, H_up); err != nil {
		panic(err)
	}
	eval.Mul(H_up, -1, H_up)
	V_up := H_up.CopyNew()

	// down init
	H_down, err := eval.MulRelinNew(state.E[cc.M-2], maskList[1][0])
	if err != nil {
		panic(err)
	}
	if err = eval.Rescale(H_down, H_down); err != nil {
		panic(err)
	}
	eval.Mul(H_down, -1, H_down)
	V_down := H_down.CopyNew()

	var up_prev, up_mid *rlwe.Ciphertext
	up_mid = V_up.CopyNew()
	eval.Mul(up_mid, -2, up_mid)

	var down_prev, down_mid *rlwe.Ciphertext
	down_mid = V_down.CopyNew()
	eval.Mul(down_mid, -2, down_mid)

	//debug.PrintDebug(cc, cc.Params, up_mid, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)

	for i := 1; i < cc.M-1; i++ {
		// H[k], 1
		eval.Add(H_up, plainList[0][i], H_up)

		// 1, V[k]
		err := eval.Rotate(V_up, -1*cc.Slot_size/cc.M, V_up)
		if err != nil {
			panic(err)
		}
		eval.Add(V_up, plainList[0][0], V_up)

		// delta , 0, 0, 0
		E, err := eval.MulRelinNew(state.E[i], maskList[0][i])
		if err != nil {
			panic(err)
		}
		if err = eval.Rescale(E, E); err != nil {
			panic(err)
		}

		H_up, V_up = state.EvalHV(cc, E, H_up, V_up)
		// 보정
		eval.Sub(H_up, corrList[0][i], H_up)
		eval.Sub(V_up, corrList[0][i], V_up)

		H_up, _ = he.BootstrapIdeal(cc, H_up)
		V_up, _ = he.BootstrapIdeal(cc, V_up)

		//fmt.Printf("UP : %d 번째 Table Update! \n", i)
		//debug.PrintDebug(cc, cc.Params, H_up, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//debug.PrintDebug(cc, cc.Params, V_up, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)

		// complex => int
		temp := V_up.CopyNew()
		eval.Add(V_up, temp, temp)
		eval.MulRelin(temp, temp, temp)
		eval.Rescale(temp, temp)
		eval.Mul(temp, -1, temp)
		eval.Add(up_mid, temp, up_mid)

		//debug.PrintDebug(cc, cc.Params, up_mid, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)

		if i == cc.M-3 {
			up_prev = up_mid.CopyNew()
		}

		// down start
		eval.Add(H_down, plainList[1][i], H_down)
		// 1, V[k]
		err = eval.Rotate(V_down, cc.Slot_size/cc.M, V_down)
		if err != nil {
			panic(err)
		}
		eval.Add(V_down, plainList[1][0], V_down)

		// delta , 0, 0, 0
		var E_down *rlwe.Ciphertext
		if i != cc.M-1 {
			E_down, err = eval.MulRelinNew(state.E[cc.M-2-i], maskList[1][i])
		} else {
			E_down, err = eval.MulRelinNew(state.E[i], maskList[1][i])
		}
		if err != nil {
			panic(err)
		}
		if err = eval.Rescale(E_down, E_down); err != nil {
			panic(err)
		}

		H_down, V_down = state.EvalHV(cc, E_down, H_down, V_down)
		// 보정
		eval.Sub(H_down, corrList[1][i], H_down)
		eval.Sub(V_down, corrList[1][i], V_down)

		H_down, _ = he.BootstrapIdeal(cc, H_down)
		V_down, _ = he.BootstrapIdeal(cc, V_down)

		//fmt.Printf("DOWN : %d 번째 Table Update! \n", i)
		//debug.PrintDebug(cc, cc.Params, H_down, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//debug.PrintDebug(cc, cc.Params, V_down, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)

		// complex => int
		temp = V_down.CopyNew()
		eval.Add(V_down, temp, temp)
		eval.MulRelin(temp, temp, temp)
		eval.Rescale(temp, temp)
		eval.Mul(temp, -1, temp)
		eval.Add(down_mid, temp, down_mid)

		if i == cc.M-3 {
			down_prev = down_mid.CopyNew()
		}
	}

	up_prev_values := make([]complex128, cc.Slot_size)
	down_prev_values := make([]complex128, cc.Slot_size)
	for j := 0; j < cc.Batch_size; j++ {

		temp := make([]complex128, cc.M)

		for k := 0; k < cc.M-2; k++ {
			temp[k] = complex(float64(k+1), 0)
		}

		for k := 0; k < cc.M; k++ {
			up_prev_values[j*cc.M+k] = temp[k]
		}

		temp = Reversed(temp)

		for k := 0; k < cc.M; k++ {
			down_prev_values[j*cc.M+k] = temp[k]
		}
	}

	up_mid_values := make([]complex128, cc.Slot_size)
	down_mid_values := make([]complex128, cc.Slot_size)
	for j := 0; j < cc.Batch_size; j++ {

		temp := make([]complex128, cc.M)

		for k := 0; k < cc.M-1; k++ {
			temp[k] = complex(float64(k+1), 0)
		}

		for k := 0; k < cc.M; k++ {
			up_mid_values[j*cc.M+k] = temp[k]
		}

		temp = Reversed(temp)

		for k := 0; k < cc.M; k++ {
			down_mid_values[j*cc.M+k] = temp[k]
		}
	}

	up_prev_values = he.InterLeave(up_prev_values, cc.M, cc.Slot_size)
	down_prev_values = he.InterLeave(down_prev_values, cc.M, cc.Slot_size)
	up_mid_values = he.InterLeave(up_mid_values, cc.M, cc.Slot_size)
	down_mid_values = he.InterLeave(down_mid_values, cc.M, cc.Slot_size)

	eval.Add(up_prev, up_prev_values, up_prev)
	eval.Add(down_prev, down_prev_values, down_prev)
	eval.Add(up_mid, up_mid_values, up_mid)
	eval.Add(down_mid, down_mid_values, down_mid)

	// fmt.Println(up_prev.Level(), down_prev.Level())
	//debug.PrintDebug(cc, cc.Params, up_prev, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
	//debug.PrintDebug(cc, cc.Params, down_prev, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)

	// debug.PrintDebug(cc, cc.Params, up_mid, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
	// debug.PrintDebug(cc, cc.Params, down_mid, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)

	min_prev := state.Min_prev(cc, up_prev, down_prev)
	min_mid := state.Min_mid(cc, up_mid, down_mid)
	min_final := state.Min_final(cc, min_prev, min_mid)

	_, res_sim := debug.PrintDebugSpecial(cc, cc.Params, min_final, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
	state.Res_sim = res_sim
	return
}

func (state *DPState) Min_prev(cc he.Context, up_prev *rlwe.Ciphertext, down_prev *rlwe.Ciphertext) (min_prev *rlwe.Ciphertext) {
	values1 := he.InverseInterLeave(he.DEC(cc, up_prev), cc.M, cc.Slot_size)
	values2 := he.InverseInterLeave(he.DEC(cc, down_prev), cc.M, cc.Slot_size)
	eq := he.InverseInterLeave(he.DEC(cc, state.E[cc.M-1]), cc.M, cc.Slot_size)

	values := make([]complex128, cc.Slot_size)
	fmt.Println(cc.Batch_size)
	for j := 0; j < cc.Batch_size; j++ {
		min_prev := 1000000.0
		for k := 0; k < cc.M-2; k++ {
			if j == 0 {
				//fmt.Println(eq[j*cc.M+k])
				//fmt.Println(real(values1[j*cc.M+k]), real(values2[j*cc.M+k+2]), (real(eq[j*cc.M+k+1]) + 1.0))
				//fmt.Println(real(values1[j*cc.M+k]) + real(values2[j*cc.M+k+2]) + (real(eq[j*cc.M+k]) + 1.0))
			}
			min_prev = math.Min(min_prev, real(values1[j*cc.M+k])+real(values2[j*cc.M+k+2])+(real(eq[j*cc.M+k+1])+1.0))
		}
		values[cc.M*j] = complex(min_prev, 0)
		//fmt.Println(min_prev)
	}

	min_prev, _ = he.ENC(cc, he.InterLeave(values, cc.M, cc.Slot_size), cc.Max_lv)
	return min_prev
}

func (state *DPState) Min_mid(cc he.Context, up_mid *rlwe.Ciphertext, down_mid *rlwe.Ciphertext) (min_mid *rlwe.Ciphertext) {
	values1 := he.InverseInterLeave(he.DEC(cc, up_mid), cc.M, cc.Slot_size)
	values2 := he.InverseInterLeave(he.DEC(cc, down_mid), cc.M, cc.Slot_size)

	values := make([]complex128, cc.Slot_size)
	fmt.Println(cc.Batch_size)
	for j := 0; j < cc.Batch_size; j++ {
		min_mid := 1000000.0
		for k := 0; k < cc.M-1; k++ {
			if j == 0 {
				//fmt.Println(eq[j*cc.M+k])
				//fmt.Println(real(values1[j*cc.M+k]), real(values2[j*cc.M+k+1]))
				//fmt.Println(real(values1[j*cc.M+k]) + real(values2[j*cc.M+k+1]))
			}
			min_mid = math.Min(min_mid, real(values1[j*cc.M+k])+real(values2[j*cc.M+k+1]))
		}
		values[cc.M*j] = complex(min_mid, 0)
		//fmt.Println(min_mid)
	}

	min_mid, _ = he.ENC(cc, he.InterLeave(values, cc.M, cc.Slot_size), cc.Max_lv)
	return min_mid
}

func (state *DPState) Min_final(cc he.Context, min_prev *rlwe.Ciphertext, min_mid *rlwe.Ciphertext) (min_final *rlwe.Ciphertext) {
	values1 := he.InverseInterLeave(he.DEC(cc, min_prev), cc.M, cc.Slot_size)
	values2 := he.InverseInterLeave(he.DEC(cc, min_mid), cc.M, cc.Slot_size)

	values := make([]complex128, cc.Slot_size)
	fmt.Println(cc.Batch_size)
	for j := 0; j < cc.Batch_size; j++ {
		min_mid := 1000000.0
		for k := 0; k < 1; k++ {
			if j == 0 {
				//fmt.Println(eq[j*cc.M+k])
				//fmt.Println(real(values1[j*cc.M+k]) + real(values2[j*cc.M+cc.M-1-k]))
				//fmt.Println(real(values1[j*cc.M+k]), real(values2[j*cc.M+k]))
			}

			min_mid = math.Min(real(values1[j*cc.M+k]), real(values2[j*cc.M+k]))
		}
		values[cc.M*j+cc.M-1] = complex(min_mid, 0)
		//fmt.Println(min_mid)
	}

	min_final, _ = he.ENC(cc, he.InterLeave(values, cc.M, cc.Slot_size), cc.Max_lv)
	return min_final
}

func (state *DPState) DpOnce(cc he.Context) {

	eval := cc.Eval

	// ones
	plainList := make([]*rlwe.Plaintext, cc.M)
	for i := 0; i < cc.M; i++ {
		values := make([]complex128, cc.Slot_size)
		for j := 0; j < cc.Slot_size; j++ {
			if j%cc.M == i {
				values[j] = complex(0.0, 0.5)
			}
		}

		values = he.InterLeave(values, cc.M, cc.Slot_size)
		plainList[i] = ckks.NewPlaintext(cc.Params, cc.Max_lv)
		if err := cc.Encoder.Encode(values, plainList[i]); err != nil { // values: []complex128 또는 []float64
			panic(err)
		}
	}

	// mask
	maskList := make([]*rlwe.Plaintext, cc.M)
	for i := 0; i < cc.M; i++ {
		values := make([]complex128, cc.Slot_size)
		for j := 0; j < cc.Batch_size; j++ {
			for k := 0; k < cc.M; k++ {
				if k <= i {
					if i == 0 {
						values[j*cc.M+k] = complex(0.5, 0.0)
					} else {
						values[j*cc.M+k] = complex(1.0, 0.0)
					}
				}
			}
		}

		values = he.InterLeave(values, cc.M, cc.Slot_size)
		maskList[i] = ckks.NewPlaintext(cc.Params, cc.Params.MaxLevel())
		if err := cc.Encoder.Encode(values, maskList[i]); err != nil { // values: []complex128 또는 []float64
			panic(err)
		}
	}

	// 보정
	corrList := make([]*rlwe.Plaintext, cc.M)
	for i := 0; i < cc.M; i++ {
		values := make([]complex128, cc.Slot_size)
		for j := 0; j < cc.Batch_size; j++ {
			for k := 0; k < cc.M; k++ {
				if i < k {
					values[j*cc.M+k] = complex(0.0, 0.5)
				}
			}
		}

		values = he.InterLeave(values, cc.M, cc.Slot_size)
		corrList[i] = ckks.NewPlaintext(cc.Params, cc.Params.MaxLevel())
		if err := cc.Encoder.Encode(values, corrList[i]); err != nil { // values: []complex128 또는 []float64
			panic(err)
		}
	}

	// forward start
	H, err := eval.MulRelinNew(state.E[0], maskList[0])
	if err != nil {
		panic(err)
	}
	if err = eval.Rescale(H, H); err != nil {
		panic(err)
	}
	eval.Mul(H, -1, H)

	V := H.CopyNew()

	for i := 1; i < cc.M; i++ {
		// H[k], 1
		eval.Add(H, plainList[i], H)

		// 1, V[k]
		err := eval.Rotate(V, -1*cc.Slot_size/cc.M, V)
		if err != nil {
			panic(err)
		}
		eval.Add(V, plainList[0], V)

		// delta , 0, 0, 0
		E, err := eval.MulRelinNew(state.E[i], maskList[i])
		if err != nil {
			panic(err)
		}
		if err = eval.Rescale(E, E); err != nil {
			panic(err)
		}
		//debug.PrintDebug(cc, cc.Params, E, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)

		//debug.PrintDebug(cc, cc.Params, H, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//debug.PrintDebug(cc, cc.Params, V, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//fmt.Println("입력 H,V level ", H.Level(), V.Level())

		H, V = state.EvalHV(cc, E, H, V)
		// 보정
		eval.Sub(H, corrList[i], H)
		eval.Sub(V, corrList[i], V)

		//debug.PrintDebug(cc, cc.Params, H, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//debug.PrintDebug(cc, cc.Params, V, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//fmt.Println(E.Level(), H.Level(), V.Level())
		H, _ = he.BootstrapIdeal(cc, H)
		V, _ = he.BootstrapIdeal(cc, V)
		//fmt.Println(E.Level(), H.Level(), V.Level())
		//fmt.Printf("%d 번째 Table Update! \n", i)
		//debug.PrintDebug(cc, cc.Params, H, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//debug.PrintDebug(cc, cc.Params, V, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
	}

	///// 아래부터 debugging
	// mask
	maskList2 := make([]*rlwe.Plaintext, cc.M)
	for i := 0; i < cc.M; i++ {
		values := make([]complex128, cc.Slot_size)
		for j := 0; j < cc.Batch_size; j++ {
			for k := 0; k < cc.M; k++ {
				if k > i {
					values[j*cc.M+k] = complex(1.0, 0.0)
				}
			}
		}

		values = he.InterLeave(values, cc.M, cc.Slot_size)
		maskList2[i] = ckks.NewPlaintext(cc.Params, cc.Params.MaxLevel())
		if err := cc.Encoder.Encode(values, maskList2[i]); err != nil { // values: []complex128 또는 []float64
			panic(err)
		}
	}

	//acc := V.CopyNew()
	acc := V.CopyNew()
	eval.Add(V, acc, acc)
	eval.MulRelin(acc, acc, acc)
	eval.Rescale(acc, acc)
	eval.Mul(acc, -1, acc)
	eval.Add(acc, cc.M, acc)

	for i := cc.M; i < 2*cc.M-1; i++ {
		// H[k], 1
		//eval.Add(H, plainList[i], H)

		// 1, V[k]
		err := eval.Rotate(V, -1*cc.Slot_size/cc.M, V)
		if err != nil {
			panic(err)
		}
		//eval.Add(V, plainList[0], V)

		// delta , 0, 0, 0
		E, err := eval.MulRelinNew(state.E[i-cc.M], maskList2[i-cc.M])
		if err != nil {
			panic(err)
		}
		if err = eval.Rescale(E, E); err != nil {
			panic(err)
		}
		//debug.PrintDebug(cc, cc.Params, E, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)

		//
		//debug.PrintDebug(cc, cc.Params, H, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//debug.PrintDebug(cc, cc.Params, V, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//fmt.Println("입력 H,V level ", H.Level(), V.Level())

		H, V = state.EvalHV(cc, E, H, V)

		temp := V.CopyNew()
		eval.Add(V, temp, temp)
		eval.MulRelin(temp, temp, temp)
		eval.Rescale(temp, temp)
		eval.Mul(temp, -1, temp)
		eval.Add(acc, temp, acc)

		//debug.PrintDebug(cc, cc.Params, H, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//debug.PrintDebug(cc, cc.Params, V, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//fmt.Println(E.Level(), H.Level(), V.Level())
		H, _ = he.BootstrapIdeal(cc, H)
		V, _ = he.BootstrapIdeal(cc, V)
		//fmt.Println(E.Level(), H.Level(), V.Level())
		//fmt.Printf("%d 번째 Table Update! \n", i)
		//debug.PrintDebug(cc, cc.Params, H, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
		//debug.PrintDebug(cc, cc.Params, V, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
	}

	_, res_sim := debug.PrintDebugSpecial(cc, cc.Params, acc, make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
	state.Res_sim = res_sim

	return
}

func (state *DPState) EvalHV(cc he.Context, delta *rlwe.Ciphertext, H_prev *rlwe.Ciphertext, V_prev *rlwe.Ciphertext) (H_curr *rlwe.Ciphertext, V_curr *rlwe.Ciphertext) {
	eval := cc.Eval
	var err error

	// conj
	H_conj, _ := eval.ConjugateNew(H_prev)
	V_conj, _ := eval.ConjugateNew(V_prev)

	// double Re
	A, _ := eval.AddNew(H_prev, H_conj)
	B, _ := eval.AddNew(V_prev, V_conj)

	// P = A * B
	P, _ := eval.MulRelinNew(A, B)
	if err = eval.Rescale(P, P); err != nil {
		panic(err)
	}

	// R = A + B - P
	R, _ := eval.AddNew(A, B)
	eval.Sub(R, P, R)

	// G_H
	eval.Mul(V_conj, complex(0, 1), V_conj)
	eval.Sub(V_conj, complex(0, 0.5), V_conj)
	G_H, _ := eval.AddNew(V_prev, V_conj)

	// G_V
	eval.Mul(H_conj, complex(0, 1), H_conj)
	eval.Sub(H_conj, complex(0, 0.5), H_conj)
	G_V, _ := eval.AddNew(H_prev, H_conj)

	E_H, _ := eval.MulRelinNew(G_H, delta)
	if err = eval.Rescale(E_H, E_H); err != nil {
		panic(err)
	}
	E_V, _ := eval.MulRelinNew(G_V, delta)
	if err = eval.Rescale(E_V, E_V); err != nil {
		panic(err)
	}

	S_H, _ := eval.AddNew(G_H, E_H)
	eval.MulRelin(S_H, R, S_H)
	if err = eval.Rescale(S_H, S_H); err != nil {
		panic(err)
	}
	S_V, _ := eval.AddNew(G_V, E_V)
	eval.MulRelin(S_V, R, S_V)
	if err = eval.Rescale(S_V, S_V); err != nil {
		panic(err)
	}

	H_curr, _ = eval.SubNew(V_prev, complex(0, 0.5))
	eval.Add(H_curr, E_H, H_curr)
	eval.Sub(S_H, H_curr, H_curr)

	V_curr, _ = eval.SubNew(H_prev, complex(0, 0.5))
	eval.Add(V_curr, E_V, V_curr)
	eval.Sub(S_V, V_curr, V_curr)

	return H_curr, V_curr
}

func EvalEQ_exp(cc he.Context, hermite int, ctxt1 *rlwe.Ciphertext, ctxt2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	// debug option! : 0 (ciphertext), 1 (plaintext), 2 (validation)
	debug_opt := debug.Debug.EQ

	if debug_opt == 0 {
		return dp_sim.EvalEQ_exp2(cc, hermite, ctxt1, ctxt2)
	} else if debug_opt == 1 {
		return dp_ideal.EvalEQ_exp(cc, hermite, ctxt1, ctxt2)
	} else {
		ct, err := dp_sim.EvalEQ_exp2(cc, hermite, ctxt1, ctxt2)
		temp, err := dp_ideal.EvalEQ_exp(cc, hermite, ctxt1, ctxt2)
		if err != nil {
			panic(err)
		}
		debug.PrintDebug(cc, cc.Params, ct, he.DEC(cc, temp), cc.Decryptor, cc.Encoder)

		return ct, err
	}
}

type DPstate_mid struct {
	CurrAntiDiagIdx int
	m               int
	n               int
	E               []*rlwe.Ciphertext // m ciphertext
	H               []*rlwe.Ciphertext // m ciphertext
	V               []*rlwe.Ciphertext // m ciphertext
	dist            *rlwe.Ciphertext
}
