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

func BootstrapIdeal(cc Context, ciphertext *rlwe.Ciphertext) (ct *rlwe.Ciphertext, err error) {
	return ENC(cc, DEC(cc, ciphertext), cc.Max_lv-9)
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

func DiBootstrapSpecial(cc Context, ciphertextList []*rlwe.Ciphertext) (ctOutList []*rlwe.Ciphertext) {

	params := cc.Params
	eval := cc.Eval
	base := 81
	//hermite_order := cc.H_btp
	polyeval := cc.PolyEval
	var err error

	ciphertext, _ := eval.AddNew(ciphertextList[0], 1)
	pow := 3
	for i := 1; i < len(ciphertextList); i++ {
		eval.Add(ciphertextList[i], 1, ciphertextList[i])
		eval.Mul(ciphertextList[i], pow, ciphertextList[i])
		eval.Add(ciphertext, ciphertextList[i], ciphertext)
		pow *= 3
	}

	// if ciphertext.Level() > 4 {
	// 	eval.DropLevel(ciphertext, ciphertext.Level()-4)
	// }

	//t := time.Now()
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

	//fmt.Println("S2C+ModUp : ", time.Since(t))

	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.

	//t = time.Now()
	var real, imag *rlwe.Ciphertext
	if real, imag, err = eval.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}
	//fmt.Println("C2S : ", time.Since(t))

	// Step 6 : EvalMod (Homomorphic modular reduction)

	// compute real part

	//t = time.Now()
	if imag, err = eval.EvalSinAndScale(real, 2*math.Pi); err != nil {
		panic(err)
	}

	if real, err = eval.EvalCosAndScale(real, 2*math.Pi); err != nil {
		panic(err)
	}

	//fmt.Println("EvalExp : ", time.Since(t))

	// Recombines the real and imaginary part
	if err = eval.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if ciphertext, err = eval.Evaluator.AddNew(real, imag); err != nil {
		panic(err)
	}

	f1 := func(x int) complex128 {
		re := x % 3
		im := ((x - (x % 3)) / 3) % 3

		re -= 1
		im -= 1

		return complex(float64(re)/8.0, float64(im)/8.0)
	}
	f2 := func(x int) complex128 {
		x = (x - (x % 3)) / 3
		x = (x - (x % 3)) / 3

		re := x % 3
		im := ((x - (x % 3)) / 3) % 3

		re -= 1
		im -= 1

		return complex(float64(re)/8.0, float64(im)/8.0)
	}
	// f3 := func(x int) complex128 {
	// 	x = (x - (x % 3)) / 3
	// 	x = (x - (x % 3)) / 3
	// 	if x%3 == 1 {
	// 		return complex(0, 0)
	// 	} else if x%3 == 2 {
	// 		return complex(0, 0.5)
	// 	} else {
	// 		return complex(0.5, 0)
	// 	}
	// }
	// f4 := func(x int) complex128 {
	// 	x = (x - (x % 3)) / 3
	// 	x = (x - (x % 3)) / 3
	// 	x = (x - (x % 3)) / 3
	// 	if x%3 == 1 {
	// 		return complex(0, 0)
	// 	} else if x%3 == 2 {
	// 		return complex(0, 0.5)
	// 	} else {
	// 		return complex(0.5, 0)
	// 	}
	// }

	eval_poly_list := []bignum.Polynomial{
		bignum.NewPolynomial(0, HermiteInterpolation3(base, 2, f1), nil),
		bignum.NewPolynomial(0, HermiteInterpolation3(base, 2, f2), nil),
		//bignum.NewPolynomial(0, HermiteInterpolation3(base, 1, f3), nil),
		//bignum.NewPolynomial(0, HermiteInterpolation3(base, 1, f4), nil),
	}

	poly_list := make([]polynomial.Polynomial, len(eval_poly_list))
	for i := 0; i < len(eval_poly_list); i++ {
		if poly_list[i] = polynomial.NewPolynomial(eval_poly_list[i]); err != nil {
			panic(err)
		}
	}

	p_slice := make([]interface{}, len(poly_list))
	for i := 0; i < len(p_slice); i++ {
		p_slice[i] = poly_list[i]
	}

	//fmt.Println(ciphertext.Level())
	//eval_t := time.Now()
	ctOutList, err = polyeval.EvaluateMultiPoly(ciphertext, p_slice, params.DefaultScale())
	if err != nil {
		panic(err)
	}
	//fmt.Println(time.Since(eval_t))
	//fmt.Println(ctOutList[0].Level())

	int2cpxCoeff := bignum.NewPolynomial(0, []complex128{0, complex(-1, 1), complex(4, 4)}, nil)
	int2cpxPoly := polynomial.NewPolynomial(int2cpxCoeff)

	extract := func(ct *rlwe.Ciphertext) (*rlwe.Ciphertext, *rlwe.Ciphertext) {
		ct_conj, _ := eval.ConjugateNew(ct)
		out1, _ := eval.AddNew(ct, ct_conj)
		out2, _ := eval.SubNew(ct, ct_conj)
		eval.Mul(out2, complex(0, -1), out2)
		out1, _ = polyeval.Evaluate(out1, int2cpxPoly, params.DefaultScale())
		out2, _ = polyeval.Evaluate(out2, int2cpxPoly, params.DefaultScale())

		return out1, out2
	}

	out1, out2 := extract(ctOutList[0])
	out3, out4 := extract(ctOutList[1])
	return []*rlwe.Ciphertext{out1, out2, out3, out4}
}

func DiBootstrapSpecial3(cc Context, ciphertext *rlwe.Ciphertext) (rtn *rlwe.Ciphertext) {

	params := cc.Params
	eval := cc.Eval
	base := cc.M
	//hermite_order := cc.H_btp
	polyeval := cc.PolyEval
	var err error

	eval.Mul(ciphertext, float64(81)/float64(base), ciphertext)
	eval.Rescale(ciphertext, ciphertext)

	// if ciphertext.Level() > 4 {
	// 	eval.DropLevel(ciphertext, ciphertext.Level()-4)
	// }

	//t := time.Now()
	if ciphertext, err = eval.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	// Step 3: scale to q/|m|
	//if ciphertext, _, err = eval.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = eval.ModUp(ciphertext); err != nil {
		panic(err)
	}

	//fmt.Println("S2C+ModUp : ", time.Since(t))

	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.

	//	t = time.Now()
	var real, imag *rlwe.Ciphertext
	if real, imag, err = eval.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}
	//fmt.Println("C2S : ", time.Since(t))

	// Step 6 : EvalMod (Homomorphic modular reduction)

	// compute real part

	//t = time.Now()
	if imag, err = eval.EvalSinAndScale(real, 2*math.Pi); err != nil {
		panic(err)
	}

	if real, err = eval.EvalCosAndScale(real, 2*math.Pi); err != nil {
		panic(err)
	}
	//fmt.Println("EvalExp : ", time.Since(t))

	// Recombines the real and imaginary part
	if err = eval.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if ciphertext, err = eval.Evaluator.AddNew(real, imag); err != nil {
		panic(err)
	}

	f1 := func(x int) complex128 {
		return complex(float64(x), 0)
	}

	coeff := bignum.NewPolynomial(0, HermiteInterpolation3(base, 1, f1), nil)
	poly := polynomial.NewPolynomial(coeff)

	//fmt.Println(ciphertext.Level())
	//eval_t := time.Now()
	rtn, err = polyeval.Evaluate(ciphertext, poly, params.DefaultScale())
	if err != nil {
		panic(err)
	}
	//fmt.Println(time.Since(eval_t))

	return rtn
}

func DiBootstrapSpecial2(cc Context, ciphertext, ciphertext2 *rlwe.Ciphertext) (rtn *rlwe.Ciphertext) {

	params := cc.Params
	eval := cc.Eval
	base := cc.M
	//hermite_order := cc.H_btp
	polyeval := cc.PolyEval
	var err error

	eval.Sub(ciphertext, ciphertext2, ciphertext)
	eval.Mul(ciphertext, float64(81)/float64(2*base), ciphertext)
	eval.Rescale(ciphertext, ciphertext)

	// if ciphertext.Level() > 4 {
	// 	eval.DropLevel(ciphertext, ciphertext.Level()-4)
	// }

	//t := time.Now()
	if ciphertext, err = eval.SlotsToCoeffs(ciphertext, nil); err != nil {
		panic(err)
	}

	// Step 3: scale to q/|m|
	//if ciphertext, _, err = eval.ScaleDown(ciphertext); err != nil {
	//	panic(err)
	//}

	// Step 4 : Extend the basis from q to Q
	if ciphertext, err = eval.ModUp(ciphertext); err != nil {
		panic(err)
	}

	//fmt.Println("S2C+ModUp : ", time.Since(t))

	// Step 5 : CoeffsToSlots (Homomorphic encoding)
	// Note: expects the result to be given in bit-reversed order
	// Also, we need the homomorphic encoding to split the real and
	// imaginary parts into two pure real ciphertexts, because the
	// homomorphic modular reduction is only defined on the reals.
	// The `imag` ciphertext can be ignored if the original input
	// is purely real.

	//t = time.Now()
	var real, imag *rlwe.Ciphertext
	if real, imag, err = eval.CoeffsToSlots(ciphertext); err != nil {
		panic(err)
	}
	//fmt.Println("C2S : ", time.Since(t))

	// Step 6 : EvalMod (Homomorphic modular reduction)

	// compute real part

	//t = time.Now()
	if imag, err = eval.EvalSinAndScale(real, 2*math.Pi); err != nil {
		panic(err)
	}

	if real, err = eval.EvalCosAndScale(real, 2*math.Pi); err != nil {
		panic(err)
	}
	//fmt.Println("EvalExp : ", time.Since(t))

	// Recombines the real and imaginary part
	if err = eval.Evaluator.Mul(imag, 1i, imag); err != nil {
		panic(err)
	}

	if ciphertext, err = eval.Evaluator.AddNew(real, imag); err != nil {
		panic(err)
	}

	f1 := func(x int) complex128 {
		if 0 <= x && x < cc.M {
			return complex(1, 0)
		} else {
			return complex(0, 0)
		}
	}

	coeff := bignum.NewPolynomial(0, HermiteInterpolation3(2*base, 1, f1), nil)
	poly := polynomial.NewPolynomial(coeff)

	//fmt.Println(ciphertext.Level())
	//eval_t := time.Now()
	rtn, err = polyeval.Evaluate(ciphertext, poly, params.DefaultScale())
	if err != nil {
		panic(err)
	}
	//fmt.Println(time.Since(eval_t))

	return rtn
}
