package dp_sim

import (
	he "banded_dp/he"
	"fmt"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

func PreComputeDeltas2(cc he.Context, hermite int, ctxt1 *rlwe.Ciphertext, ctxt2 *rlwe.Ciphertext) []*rlwe.Ciphertext {
	//sz := cc.C_size
	str_len := cc.M
	//b_sz := cc.Batch_size
	eval := cc.Eval
	var err error

	bit_start := time.Now()
	bit_ctxt1_List, _ := EvalBitDecomp(cc, 0, ctxt1)
	bit_ctxt2_List, _ := EvalBitDecomp(cc, 0, ctxt2)
	fmt.Println(bit_ctxt1_List[0].Level(), bit_ctxt2_List[0].Level())
	bit_time := time.Since(bit_start)
	fmt.Println(bit_time)

	for i := 0; i < len(bit_ctxt1_List); i++ {
		bit_ctxt1_List[i], _ = he.BootstrapIdeal(cc, bit_ctxt1_List[i])
		bit_ctxt2_List[i], _ = he.BootstrapIdeal(cc, bit_ctxt2_List[i])
	}
	fmt.Println(bit_ctxt1_List[0].Level(), bit_ctxt2_List[0].Level())
	// Note : str.A가 column에 있고 str.B가 row에 있을 예정임.
	//roted_list := make([]*rlwe.Ciphertext, str_len)
	//temp_list := make([]*rlwe.Ciphertext, str_len)
	//temp_list[0] = ctxt2

	shift_idx := []int{}
	for i := 0; i < str_len; i++ {
		shift_idx = append(shift_idx, i*cc.Slot_size/cc.M)
	}
	fmt.Println(shift_idx)

	rot_start := time.Now()

	for i := 0; i < len(bit_ctxt1_List); i++ {
		temp, err := eval.ConjugateNew(bit_ctxt2_List[i])
		if err != nil {
			panic(err)
		}
		eval.Add(bit_ctxt1_List[i], temp, temp)

		temp_list, err := eval.Evaluator.RotateHoistedNew(temp, shift_idx)
		if err != nil {
			fmt.Println(err)
			panic(err)
		}

		for k := 0; k < cc.M; k++ {
			imag, _ := eval.ConjugateNew(temp_list[(k)*cc.Slot_size/cc.M])
			real, _ := eval.AddNew(temp_list[(k)*cc.Slot_size/cc.M], imag)
			eval.Sub(temp_list[(k)*cc.Slot_size/cc.M], imag, imag)
			eval.MulRelin(real, 0.5, real)
			eval.Rescale(real, real)
			eval.MulRelin(imag, 0.5, imag)
			eval.Rescale(imag, imag)
		}
	}
	rot_time := time.Since(rot_start)

	for i := 1; i <= str_len; i++ {
		//roted_list[i-1] = temp_list[(str_len-i)*cc.Slot_size/cc.M]
		//debug.PrintDebug(cc, cc.Params, roted_list[i-1], make([]complex128, cc.Params.MaxSlots()), cc.Decryptor, cc.Encoder)
	}

	rtn_list := make([]*rlwe.Ciphertext, str_len)
	rot_start = time.Now()
	for i := 0; i < str_len; i++ {
		rot_start := time.Now()
		rtn_list[i], err = EvalEQ_exp3(cc, 0, bit_ctxt1_List, bit_ctxt2_List)
		fmt.Println("eq time : ", time.Since(rot_start))

		if err != nil {
			panic(nil)
		}
	}
	eq_time := time.Since(rot_start)

	fmt.Println("bit time : ", bit_time)
	fmt.Println("rot time : ", rot_time)
	fmt.Println("eq time : ", eq_time)

	return rtn_list
}

func PreComputeDeltas(cc he.Context, hermite int, ctxt1 *rlwe.Ciphertext, ctxt2 *rlwe.Ciphertext) []*rlwe.Ciphertext {
	if cc.M != 256 {
		//sz := cc.C_size
		str_len := cc.M
		//b_sz := cc.Batch_size
		eval := cc.Eval
		var err error

		// Note : str.A가 column에 있고 str.B가 row에 있을 예정임.
		roted_list := make([]*rlwe.Ciphertext, str_len)

		shift_idx := []int{}
		for i := 0; i < str_len; i++ {
			shift_idx = append(shift_idx, i*cc.Slot_size/cc.M)
		}
		//fmt.Println(shift_idx)

		//rot_start := time.Now()
		temp_list, err := eval.Evaluator.RotateHoistedNew(ctxt2, shift_idx)
		if err != nil {
			panic(err)
		}
		//rot_time := time.Since(rot_start)

		for i := 1; i <= str_len; i++ {
			roted_list[i-1] = temp_list[(str_len-i)*cc.Slot_size/cc.M]
		}

		rtn_list := make([]*rlwe.Ciphertext, str_len)
		//	rot_start = time.Now()
		// conjugate : exp(y) ==> exp(-y)
		err = eval.Conjugate(ctxt1, ctxt1)
		if err != nil {
			panic(err)
		}
		for i := 0; i < str_len; i++ {
			//temp := eval.DropLevelNew(ctxt1, (i%4)*2)
			temp := eval.DropLevelNew(ctxt1, 0)
			//roted_list[i] = eval.DropLevelNew(roted_list[i], (i%4)*2)
			roted_list[i] = eval.DropLevelNew(roted_list[i], 0)
			if i == 0 {
				fmt.Println("levels : ", temp.Level(), roted_list[i].Level())
			}

			rtn_list[i], err = EvalEQ_exp2(cc, 0, temp, roted_list[i])
			if i == 0 {
				fmt.Println(rtn_list[i].Level())
			}

			//fmt.Printf("equlaity index : %d && level : %d\n", i, rtn_list[i].Level())
			if err != nil {
				panic(nil)
			}
		}
		//eq_time := time.Since(rot_start)

		//fmt.Println("rot time : ", rot_time)
		//fmt.Println("eq time : ", eq_time)

		return rtn_list
	} else {
		//sz := cc.C_size
		str_len := cc.M
		//b_sz := cc.Batch_size
		eval := cc.Eval
		var err error

		// Note : str.A가 column에 있고 str.B가 row에 있을 예정임.
		roted_list := make([]*rlwe.Ciphertext, str_len)

		shift_idx := []int{}
		for i := 0; i < str_len/2; i++ {
			shift_idx = append(shift_idx, i*cc.Slot_size/cc.M)
		}
		//fmt.Println(shift_idx)

		//rot_start := time.Now()
		temp_list, err := eval.Evaluator.RotateHoistedNew(ctxt2, shift_idx)
		if err != nil {
			panic(err)
		}

		ctxt2_tmp, _ := eval.RotateNew(ctxt2, (str_len/2)*cc.Slot_size/cc.M)
		temp_list2, err := eval.Evaluator.RotateHoistedNew(ctxt2_tmp, shift_idx)
		if err != nil {
			panic(err)
		}

		//rot_time := time.Since(rot_start)

		for i := 1; i <= str_len; i++ {
			if i <= str_len/2 {
				roted_list[i-1] = temp_list2[(str_len/2-i)*cc.Slot_size/cc.M]
			} else {
				roted_list[i-1] = temp_list[(str_len-i)*cc.Slot_size/cc.M]
			}

		}

		rtn_list := make([]*rlwe.Ciphertext, str_len)
		//rot_start = time.Now()
		// conjugate : exp(y) ==> exp(-y)
		err = eval.Conjugate(ctxt1, ctxt1)
		if err != nil {
			panic(err)
		}
		for i := 0; i < str_len; i++ {
			//temp := eval.DropLevelNew(ctxt1, (i%4)*2)
			temp := eval.DropLevelNew(ctxt1, 0)
			//roted_list[i] = eval.DropLevelNew(roted_list[i], (i%4)*2)
			roted_list[i] = eval.DropLevelNew(roted_list[i], 0)
			//if i == 0 {
			//	fmt.Println("levels : ", temp.Level(), roted_list[i].Level())
			//}

			rtn_list[i], err = EvalEQ_exp2(cc, 0, temp, roted_list[i])
			if i == 0 {
				//fmt.Println(rtn_list[i].Level())
			}

			//fmt.Printf("equlaity index : %d && level : %d\n", i, rtn_list[i].Level())
			if err != nil {
				panic(nil)
			}
		}
		//eq_time := time.Since(rot_start)

		//	fmt.Println("rot time : ", rot_time)
		//fmt.Println("eq time : ", eq_time)

		return rtn_list
	}
}
