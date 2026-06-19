package dp_sim

import (
	he "banded_dp/he"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

// Note : 현재는 cipher-cipher equality 상태로 작성되어 있습니다. 추후에 수정 예정입니다.
// Input : 각 암호문은 N/2개의 8bit integer가 각각 exp(2pi*i* x/256), exp(2pi*i* y/256) form으로 저장되어 있습니다.
// Output : 출력 암호문은 x==y 이면 1을 그렇지 않으면 0이 저장됩니다.
// description : 1. Compute exp(2pi*i*(x-y)/256) : ctxt2를 conjugation 한 뒤, ctxt1과 곱합니다.
//  2. Hermite Interpoation : equality Look-up Table을 평가합니다.
func EvalEquality(ciphertext1 *rlwe.Ciphertext, ciphertext2 *rlwe.Ciphertext, cc he.Context) *rlwe.Ciphertext {

	eval := cc.Eval
	params := cc.Params
	poly_eval := cc.PolyEval

	eq := func(x int) int {
		if x == 0 {
			return 1
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

	// fmt.Printf("before mul: ct1 level=%d scale=2^%.2f | ct2 level=%d scale=2^%.2f\n",
	// 	ciphertext1.Level(), math.Log2(ciphertext1.Scale.Float64()),
	// 	ct2Conj.Level(), math.Log2(ct2Conj.Scale.Float64()),
	// )

	// mult : exp(x) * exp(-y) ==> exp(x-y)
	result, err = eval.MulRelinNew(ciphertext1, ct2Conj)
	if err != nil {
		panic(err)
	}

	err = eval.Rescale(result, result)
	if err != nil {
		panic(err)
	}

	// fmt.Printf("after rescale: result level=%d scale=2^%.2f logQ=%d\n",
	// 	result.Level(),
	// 	math.Log2(result.Scale.Float64()),
	// 	params.LogQLvl(result.Level()),
	// )

	var f_coeffs bignum.Polynomial

	// hermite order : 1
	f_coeffs = bignum.NewPolynomial(0, he.HermiteInterpolation(256, 1, eq), nil)

	var poly_f polynomial.Polynomial
	poly_f = polynomial.NewPolynomial(f_coeffs)

	// // 여기서 bootstrap
	// result, err = cc.btp.Bootstrap(result)
	// if err != nil {
	// 	panic(err)
	// }

	// fmt.Printf("after bootstrap: result level=%d scale=2^%.2f logQ=%d\n",
	// 	result.Level(),
	// 	math.Log2(result.Scale.Float64()),
	// 	params.LogQLvl(result.Level()),
	// )

	// 그 다음 polynomial evaluation
	result, err = poly_eval.Evaluate(result, poly_f, params.DefaultScale())
	if err != nil {
		panic(err)
	}

	return result
}
