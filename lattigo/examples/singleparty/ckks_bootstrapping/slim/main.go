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
	"flag"
	"fmt"
	"math"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/bignum"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
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

	//============================
	//=== 1) SCHEME PARAMETERS ===
	//============================

	// In this example, for a pratical purpose, the residual parameters and bootstrapping
	// parameters are the same. But in practice the residual parameters would not contain the
	// moduli for the CoeffsToSlots and EvalMod steps.
	// With LogN=16, LogQP=1221 and H=192, these parameters achieve well over 128-bit of security.
	// For the purpose of the example, only one prime is allocated to the circuit in the slots domain
	// and no prime is allocated to the circuit in the coeffs domain.

	// scaling factor
	LogDefaultScale := 45

	// base modulus
	q0 := []int{55}                      // 3) ScaleDown & 4) ModUp
	qiSlotsToCoeffs := []int{39, 39, 39} // 1) SlotsToCoeffs

	// len(qiCircuitsSlots) = available depth
	qiCircuitSlots := []int{45, 45, 45, 45, 45, 45, 45, 45, 45, 45} // 0) Circuit in the slot domain
	qiEvalMod := []int{60, 60, 60, 60, 60, 60, 60, 60}              // 6) EvalMod
	qiCoeffsToSlots := []int{56, 56, 56, 56}                        // 5) CoeffsToSlots

	LogQ := append(q0, qiSlotsToCoeffs...)
	LogQ = append(LogQ, qiCircuitSlots...)
	max_lv := len(LogQ) - 1
	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,                      // Log2 of the ring degree
		LogQ:            LogQ,                      // Log2 of the ciphertext modulus
		LogP:            []int{61, 61, 61, 61, 61}, // Log2 of the key-switch auxiliary prime moduli
		LogDefaultScale: LogDefaultScale,           // Log2 of the scale
		Xs:              ring.Ternary{H: 192},
	})

	if err != nil {
		panic(err)
	}

	//====================================
	//=== 2) BOOTSTRAPPING PARAMETERS ===
	//====================================

	// CoeffsToSlots parameters (homomorphic encoding)
	CoeffsToSlotsParameters := dft.MatrixLiteral{
		Type:         dft.HomomorphicEncode,
		Format:       dft.RepackImagAsReal, // Returns the real and imaginary part into separate ciphertexts
		LogSlots:     params.LogMaxSlots(),
		LevelQ:       params.MaxLevelQ(),
		LevelP:       params.MaxLevelP(),
		LogBSGSRatio: 1,
		Levels:       []int{1, 1, 1, 1}, //qiCoeffsToSlots
	}

	// Parameters of the homomorphic modular reduction x mod 1
	Mod1ParametersLiteral := mod1.ParametersLiteral{
		LevelQ:          params.MaxLevel() - CoeffsToSlotsParameters.Depth(true),
		LogScale:        60,               // Matches qiEvalMod
		Mod1Type:        mod1.CosDiscrete, // Multi-interval Chebyshev interpolation
		Mod1Degree:      30,               // Depth 5
		DoubleAngle:     3,                // Depth 3
		K:               16,               // With EphemeralSecretWeight = 32 and 2^{15} slots, ensures < 2^{-138.7} failure probability
		LogMessageRatio: 10,               // q/|m| = 2^10
		Mod1InvDegree:   0,                // Depth 0
	}

	// SlotsToCoeffs parameters (homomorphic decoding)
	SlotsToCoeffsParameters := dft.MatrixLiteral{
		Type:         dft.HomomorphicDecode,
		LogSlots:     params.LogMaxSlots(),
		LogBSGSRatio: 1,
		LevelP:       params.MaxLevelP(),
		Levels:       []int{1, 1, 1}, // qiSlotsToCoeffs
	}

	SlotsToCoeffsParameters.LevelQ = len(SlotsToCoeffsParameters.Levels)

	// Custom bootstrapping.Parameters.
	// All fields are public and can be manually instantiated.
	btpParams := bootstrapping.Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParameters,
		Mod1ParametersLiteral:   Mod1ParametersLiteral,
		CoeffsToSlotsParameters: CoeffsToSlotsParameters,
		EphemeralSecretWeight:   32, // > 128bit secure for LogN=16 and LogQP = 115.
		CircuitOrder:            bootstrapping.DecodeThenModUp,
	}

	// We pring some information about the bootstrapping parameters (which are identical to the residual parameters in this example).
	// We can notably check that the LogQP of the bootstrapping parameters is smaller than 1550, which ensures
	// 128-bit of security as explained above.
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())

	//===========================
	//=== 3) KEYGEN & ENCRYPT ===
	//===========================

	// Now that both the residual and bootstrapping parameters are instantiated, we can
	// instantiate the usual necessary object to encode, encrypt and decrypt.

	// Scheme context and keys
	kgen := rlwe.NewKeyGenerator(params)

	sk, pk := kgen.GenKeyPairNew()

	encoder := ckks.NewEncoder(params)
	decryptor := rlwe.NewDecryptor(params, sk)
	encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println()
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}
	fmt.Println("Done")

	//========================
	//=== 4) BOOTSTRAPPING ===
	//========================

	// Instantiates the evaluator

	// ckks evaluator 에는 add, mult, rotation과 같은 함수들이 있음.
	// 시간나면 member 함수 목록들을 보면서 모르는것들을 공부해나가면 좋음.
	eval := ckks.NewEvaluator(params, evk)

	// polynomial evaluation은 polynomial evaluator로 진행함.
	// 미리 계산된 계수만 있으면, 최적의 방법으로 다항식을 계산해줌.
	poly_eval := polynomial.NewEvaluator(params, eval)

	// Generate a random plaintext with values uniformly distributed in [-1, 1] for the real and imaginary part.
	valuesWant := make([]complex128, params.MaxSlots())
	debug_values := make([]int, params.MaxSlots())
	for i := range valuesWant {
		x := sampling.RandInt(big.NewInt(256))
		debug_values[i] = int(x.Int64())
		// exp(2pi*i * x / 256) = cos(2pi * x / 256) + i * sin(2pi* x / 256)
		valuesWant[i] = complex(math.Cos(2*math.Pi*(float64(x.Int64())/256.0)), math.Sin(2*math.Pi*(float64(x.Int64())/256.0)))
	}

	// We encrypt at level 0
	plaintext := ckks.NewPlaintext(params, max_lv)
	if err := encoder.Encode(valuesWant, plaintext); err != nil {
		panic(err)
	}

	// Encrypt
	ciphertext, err := encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}

	valuesWant2 := make([]complex128, params.MaxSlots())
	debug_values2 := make([]int, params.MaxSlots())
	for i := range valuesWant {
		x := sampling.RandInt(big.NewInt(256))
		debug_values2[i] = int(x.Int64())
		// exp(2pi*i * x / 256) = cos(2pi * x / 256) + i * sin(2pi* x / 256)
		valuesWant2[i] = complex(math.Cos(2*math.Pi*(float64(x.Int64())/256.0)), math.Sin(2*math.Pi*(float64(x.Int64())/256.0)))
	}

	// We encrypt at level 0
	plaintext2 := ckks.NewPlaintext(params, max_lv)
	if err := encoder.Encode(valuesWant2, plaintext2); err != nil {
		panic(err)
	}

	// Encrypt
	ciphertext2, err := encryptor.EncryptNew(plaintext2)
	if err != nil {
		panic(err)
	}

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	valuesTest := printDebug(params, ciphertext, valuesWant, decryptor, encoder)

	// polynomial evaluation

	// Target : f(x) = x^3 + 1
	// 계수를 먼저 bignum.Polynomial에 초기화 해줘야함.
	// 여기서는 degree, depth 같은 메타정보와, chebyshev basis로 변환해주는 역할을 함.

	// x == 3 ==> 1
	// x != 3 ==> 0
	eq := func(x int) int {
		if x == 0 {
			return 1
		} else {
			return 0
		}
	}

	err = eval.Conjugate(ciphertext2, ciphertext2)
	if err != nil {
		fmt.Println(err)
	}

	var result *rlwe.Ciphertext
	result, _ = eval.MulRelinNew(ciphertext, ciphertext2)
	eval.Rescale(result, result)

	var f_coeffs bignum.Polynomial
	f_coeffs = bignum.NewPolynomial(0, HermiteInterpolation(256, 1, eq), nil)

	// 그 다음에는, 이걸 polnomial 객체로 wrapping 시켜서, 계수정보를 plaintext로 encoding하는 작업을 함으로써
	// 원 버튼 딸깍용으로 만들어주는 역할.
	var poly_f polynomial.Polynomial
	poly_f = polynomial.NewPolynomial(f_coeffs)

	result, err = poly_eval.Evaluate(result, poly_f, params.DefaultScale())
	if err != nil {
		panic(err)
	}

	fmt.Println("Done")

	for i := 0; i < len(debug_values); i++ {
		if debug_values[i] == debug_values2[i] {
			valuesTest[i] = 1
		} else {
			valuesTest[i] = 0
		}
	}

	//==================
	//=== 5) DECRYPT ===
	//==================

	fmt.Println("test")

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	printDebug(params, result, valuesTest, decryptor, encoder)
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
