package main

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Note : 현재는 cipher-cipher equality 상태로 작성되어 있습니다. 추후에 수정 예정입니다.
// Input : 각 암호문은 N/2개의 8bit integer가 각각 exp(2pi*i* x/256), exp(2pi*i* y/256) form으로 저장되어 있습니다.
// Output : 출력 암호문은 x==y 이면 1을 그렇지 않으면 0이 저장됩니다.
// description : 1. Compute exp(2pi*i*(x-y)/256) : ctxt2를 conjugation 한 뒤, ctxt1과 곱합니다.
//  2. Hermite Interpoation : equality Look-up Table을 평가합니다.
func EvalEquality(ciphertext1 *rlwe.Ciphertext, ciphertext2 *rlwe.Ciphertext, cc Context) *rlwe.Ciphertext {

	eval := cc.eval
	params := cc.params
	poly_eval := cc.poly_eval

	eq := func(x int) int {
		if x == 0 {
			return 1
		} else {
			return 0
		}
	}

	var err error
	var result *rlwe.Ciphertext

	// conjugate : exp(y) ==> exp(-y)
	err = eval.Conjugate(ciphertext2, ciphertext2)
	if err != nil {
		fmt.Println(err)
	}

	// mult : exp(x) * exp(-y) ==> exp(x-y)
	result, _ = eval.MulRelinNew(ciphertext1, ciphertext2)
	eval.Rescale(result, result)

	var f_coeffs bignum.Polynomial

	// hermite order : 1
	f_coeffs = bignum.NewPolynomial(0, HermiteInterpolation(256, 1, eq), nil)
	var poly_f polynomial.Polynomial
	poly_f = polynomial.NewPolynomial(f_coeffs)

	// eval poly : exp(x-y) ==> eq(x-y)
	result, err = poly_eval.Evaluate(result, poly_f, params.DefaultScale())
	if err != nil {
		panic(err)
	}

	return result
}
