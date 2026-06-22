package he

import (
	"math"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
)

func must(err error) {
	if err != nil {
		panic(err)
	}
}

// DiBootstrapMany with Multiple Look-up Table Evaluation
// 이 함수는 하나의 암호문에 대해서 다수의 Look-Up Table을 평가하여 출력하는 함수입니다.
// 세부사항은 DiBootstrap 함수를 참고하세요.
func DiBootstrapMany(cc Context, ciphertext *rlwe.Ciphertext, LookUpTable func(int) int) (ctOut *rlwe.Ciphertext, err error) {

	params := cc.Params
	eval := cc.Eval
	base := 1 << cc.LogB
	hermite_order := cc.H_btp
	polyeval := cc.PolyEval

	if ciphertext, err = eval.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	// Step 3: scale to q/|m|
	if ciphertext, _, err = eval.ScaleDown(ciphertext); err != nil {
		panic(err)
	}

	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = eval.ModUp(ciphertext); err != nil {
		panic(err)
	}

	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.
	var real, imag *rlwe.Ciphertext
	if real, imag, err = eval.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}

	// Step 6 : EvalMod (Homomorphic modular reduction)

	// compute real part

	if imag, err = eval.EvalSinAndScale(real, 2*math.Pi); err != nil {
		panic(err)
	}

	if real, err = eval.EvalCosAndScale(real, 2*math.Pi); err != nil {
		panic(err)
	}

	// Recombines the real and imaginary part
	if err = eval.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if ciphertext, err = eval.Evaluator.AddNew(real, imag); err != nil {
		panic(err)
	}

	eval_poly := bignum.NewPolynomial(0, HermiteInterpolation(base, hermite_order, LookUpTable), nil)
	poly := polynomial.NewPolynomial(eval_poly)

	ctOut, err = polyeval.Evaluate(ciphertext, poly, params.DefaultScale())

	return ctOut, err
}
