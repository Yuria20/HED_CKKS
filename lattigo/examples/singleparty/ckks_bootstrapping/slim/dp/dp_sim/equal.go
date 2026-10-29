package dp_sim

import (
	he "banded_dp/he"
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Note : 현재는 cipher-cipher equality 상태로 작성되어 있습니다. 추후에 수정 예정입니다.
// Input : 각 암호문은 N/2개의 7bit integer가 각각 exp(2pi*i* x/128), exp(2pi*i* y/128) form으로 저장되어 있습니다.
// Output : 출력 암호문은 x==y 이면 1을 그렇지 않으면 0이 저장됩니다.
// description : 1. Compute exp(2pi*i*(x-y)/128) : ctxt2를 conjugation 한 뒤, ctxt1과 곱합니다.
//  2. Hermite Interpoation : equality Look-up Table을 평가합니다.
func EvalEQ_exp(cc he.Context, hermite int, ciphertext1 *rlwe.Ciphertext, ciphertext2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	eval := cc.Eval
	params := cc.Params
	poly_eval := cc.PolyEval
	sz := cc.C_size

	eq := func(x int) int {
		if x == 0 {
			return -1
		} else {
			return 0
		}
	}

	var err error
	var result *rlwe.Ciphertext

	// ciphertext2 복사본
	ct2Conj := ciphertext2.CopyNew()

	// conjugate : exp(y) ==> exp(-y)
	err = eval.Conjugate(ct2Conj, ct2Conj)
	if err != nil {
		panic(err)
	}

	// mult : exp(x) * exp(-y) ==> exp(x-y)
	result, err = eval.MulRelinNew(ciphertext1, ct2Conj)
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(result, result)
	if err != nil {
		panic(err)
	}

	var f_coeffs bignum.Polynomial

	// hermite order :
	f_coeffs = bignum.NewPolynomial(0, he.HermiteInterpolation(sz, hermite, eq), nil)

	var poly_f polynomial.Polynomial
	poly_f = polynomial.NewPolynomial(f_coeffs)

	result, err = poly_eval.Evaluate(result, poly_f, params.DefaultScale())
	if err != nil {
		panic(err)
	}

	return result, err
}

// Note : 현재는 cipher-cipher equality 상태로 작성되어 있습니다. 추후에 수정 예정입니다.
// Input : 각 암호문은 N/2개의 7bit integer가 각각 exp(2pi*i* x/128), exp(2pi*i* y/128) form으로 저장되어 있습니다.
// Output : 출력 암호문은 x==y 이면 1을 그렇지 않으면 0이 저장됩니다.
// description : 1. Compute exp(2pi*i*(x-y)/128) : ctxt2를 conjugation 한 뒤, ctxt1과 곱합니다.
//  2. Hermite Interpoation : equality Look-up Table을 평가합니다.
func EvalEQ_exp2(cc he.Context, hermite int, ciphertext1 *rlwe.Ciphertext, ciphertext2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	eval := cc.Eval
	sz := cc.C_size
	log_sz := int(math.Log2(float64(sz)))

	var err error
	var result *rlwe.Ciphertext

	// mult : exp(x) * exp(-y) ==> exp(x-y)
	result, err = eval.MulRelinNew(ciphertext1, ciphertext2)
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(result, result)
	if err != nil {
		panic(err)
	}

	powers := make([]*rlwe.Ciphertext, log_sz)
	powers[0] = result
	for i := 1; i < log_sz; i++ {
		powers[i], err = eval.MulRelinNew(powers[i-1], powers[i-1])
		if err != nil {
			panic(err)
		}
		err = eval.Rescale(powers[i], powers[i])
		if err != nil {
			panic(err)
		}
	}

	result, err = eval.AddNew(powers[0], 1)
	if err != nil {
		panic(err)
	}
	for i := 1; i < log_sz; i++ {
		temp, err := eval.AddNew(powers[i], 1)
		if err != nil {
			panic(err)
		}
		result, err = eval.MulRelinNew(result, temp)
		if err != nil {
			panic(err)
		}
		err = eval.Rescale(result, result)
		if err != nil {
			panic(err)
		}
	}

	result, err = eval.MulRelinNew(result, -1.0/float64(sz))
	if err != nil {
		panic(err)
	}
	err = eval.Rescale(result, result)
	if err != nil {
		panic(err)
	}

	return result, err
}

func EvalEQ_exp3(cc he.Context, hermite int, ctList []*rlwe.Ciphertext, ctList2 []*rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	eval := cc.Eval
	// params := cc.Params
	// poly_eval := cc.PolyEval
	//sz := cc.C_size
	//log_sz := int(math.Log2(float64(sz)))

	// eq := func(x int) int {
	// 	if x == 0 {
	// 		return 1
	// 	} else {
	// 		return 0
	// 	}
	// }

	var err error
	var result *rlwe.Ciphertext

	ciphertextList := make([]*rlwe.Ciphertext, len(ctList))
	ciphertextList2 := make([]*rlwe.Ciphertext, len(ctList))

	//cp_start := time.Now()
	for i := 0; i < len(ctList); i++ {
		ciphertextList[i] = ctList[i].CopyNew()
		ciphertextList2[i] = ctList2[i].CopyNew()
	}
	//fmt.Println(time.Since(cp_start))

	for i := 0; i < len(ciphertextList); i++ {
		//fmt.Println(i)
		eval.MulRelin(ciphertextList[i], ciphertextList2[i], ciphertextList[i])
		if err != nil {
			panic(err)
		}

		err = eval.Rescale(ciphertextList[i], ciphertextList[i])
		if err != nil {
			panic(err)
		}

	}

	for i := 0; i < len(ciphertextList); i++ {
		//fmt.Println(ciphertext1[i].Level())
	}

	eval.MulRelin(ciphertextList[0], ciphertextList[1], ciphertextList[0])
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(ciphertextList[0], ciphertextList[0])
	if err != nil {
		panic(err)
	}

	eval.MulRelin(ciphertextList[2], ciphertextList[3], ciphertextList[1])
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(ciphertextList[1], ciphertextList[1])
	if err != nil {
		panic(err)
	}

	eval.MulRelin(ciphertextList[4], ciphertextList[5], ciphertextList[2])
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(ciphertextList[2], ciphertextList[2])
	if err != nil {
		panic(err)
	}

	for i := 0; i < len(ciphertextList); i++ {
		//fmt.Println(ciphertext1[i].Level())
	}

	eval.MulRelin(ciphertextList[0], ciphertextList[1], ciphertextList[0])
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(ciphertextList[0], ciphertextList[0])
	if err != nil {
		panic(err)
	}

	eval.MulRelin(ciphertextList[2], ciphertextList[6], ciphertextList[1])
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(ciphertextList[1], ciphertextList[1])
	if err != nil {
		panic(err)
	}

	result, _ = eval.MulRelinNew(ciphertextList[0], ciphertextList[1])
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(result, result)
	if err != nil {
		panic(err)
	}

	fmt.Println(result.Level())
	return result, err
}

func EvalBitDecomp(cc he.Context, hermite int, ciphertext1 *rlwe.Ciphertext) ([]*rlwe.Ciphertext, error) {

	eval := cc.Eval
	params := cc.Params
	poly_eval := cc.PolyEval
	sz := cc.C_size

	bit01 := func(x int) (int, int) {
		return (x >> 0) & 1, (x >> 1) & 1
	}
	bit23 := func(x int) (int, int) {
		return (x >> 2) & 1, (x >> 3) & 1
	}
	bit45 := func(x int) (int, int) {
		return (x >> 4) & 1, (x >> 5) & 1
	}
	bit6 := func(x int) int {
		return (x >> 6) & 1
	}

	var err error
	var bitPack []*rlwe.Ciphertext

	// hermite order :
	eval_poly_list := []bignum.Polynomial{
		bignum.NewPolynomial(0, he.HermiteInterpolation2(sz, hermite, bit01), nil),
		bignum.NewPolynomial(0, he.HermiteInterpolation2(sz, hermite, bit23), nil),
		bignum.NewPolynomial(0, he.HermiteInterpolation2(sz, hermite, bit45), nil),
		bignum.NewPolynomial(0, he.HermiteInterpolation(sz, hermite, bit6), nil),
	}

	poly_list := make([]polynomial.Polynomial, 4)
	for i := 0; i < 4; i++ {
		if poly_list[i] = polynomial.NewPolynomial(eval_poly_list[i]); err != nil {
			panic(err)
		}
	}

	p_slice := make([]interface{}, 4)
	for i := 0; i < len(p_slice); i++ {
		p_slice[i] = poly_list[i]
	}

	bitPack, err = poly_eval.EvaluateMultiPoly(ciphertext1, p_slice, params.DefaultScale())
	if err != nil {
		panic(err)
	}

	rtn := make([]*rlwe.Ciphertext, 7)

	for i := 0; i < 3; i++ {
		conj, _ := eval.ConjugateNew(bitPack[i])
		rtn[2*i], _ = eval.AddNew(bitPack[i], conj)
		rtn[2*i+1], _ = eval.SubNew(bitPack[i], conj)
		eval.Mul(rtn[2*i+1], complex(0, -1), rtn[2*i+1])
	}
	rtn[6] = bitPack[3]

	for i := 0; i < 7; i++ {
		//debug.PrintDebug(cc, cc.Params, rtn[i], make([]complex128, cc.Slot_size), cc.Decryptor, cc.Encoder)
	}

	return rtn, err
}
