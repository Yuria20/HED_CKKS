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
	"math/rand"
	"time"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

type Context struct {
	params    *ckks.Parameters
	encoder   *ckks.Encoder
	encryptor *rlwe.Encryptor
	decryptor *rlwe.Decryptor
	eval      *ckks.Evaluator
	btp       *bootstrapping.Evaluator
	poly_eval *polynomial.Evaluator
	max_lv    int
}

func main() {

	flag.Parse()
	rand.Seed(time.Now().UnixNano())

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
	qiCircuitSlots := []int{45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45} // 0) Circuit in the slot domain
	qiEvalMod := []int{60, 60, 60, 60, 60, 60, 60, 60}                                              // 6) EvalMod
	qiCoeffsToSlots := []int{56, 56, 56, 56}                                                        // 5) CoeffsToSlots

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

	lengths := []int{8}
	N := 4 // target 개수
	maxM := 0
	for _, m := range lengths {
		if m > maxM {
			maxM = m
		}
	}

	fmt.Println()
	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}

	rots := needRotsExpand(params.MaxSlots(), maxM, N)

	galEls := make([]uint64, 0, len(rots))
	for _, r := range rots {
		galEls = append(galEls, params.GaloisElement(r))
	}

	gk := kgen.GenGaloisKeysNew(galEls, sk)
	for i, key := range gk {
		evk.GaloisKeys[galEls[i]] = key
	}

	eval := ckks.NewEvaluator(params, evk)

	var btp *bootstrapping.Evaluator

	if btp, err = bootstrapping.NewEvaluator(btpParams, evk); err != nil {
		panic(err)
	}

	poly_eval := polynomial.NewEvaluator(params, eval)

	var cc Context
	cc.params = &params
	cc.encoder = encoder
	cc.encryptor = encryptor
	cc.decryptor = decryptor
	cc.eval = eval
	cc.btp = btp
	cc.poly_eval = poly_eval
	cc.max_lv = max_lv

	for _, m := range lengths {
		fmt.Printf("\n==============================\n")
		fmt.Printf("packed 1:N experiment | m=%d | N=%d\n", m, N)
		fmt.Printf("==============================\n")

		query := RandomStringInts(m, AlphabetSize)
		targets := RandomSetMatch(query, N, AlphabetSize)

		fmt.Println("query =", query)
		for i := 0; i < len(targets); i++ {
			fmt.Printf("target[%d] = %v\n", i, targets[i])
		}

		dt, err := Delta1N(query, targets, cc)
		if err != nil {
			panic(err)
		}

		fmt.Println("m =", dt.M)
		fmt.Println("pair count =", dt.PairCount)
		fmt.Println("num blocks =", len(dt.DeltaBlocks))


		for blk, ct := range dt.DeltaBlocks {
			fmt.Printf("\n=== Delta block %d ===\n", blk)

			slots := ct.Slots()
			if !ct.IsBatched {
				slots *= 2
			}

			values := make([]complex128, slots)
			if err := encoder.Decode(decryptor.DecryptNew(ct), values); err != nil {
				panic(err)
			}

			limit := len(dt.Valid[blk])
			if len(values) < limit {
				limit = len(values)
			}
			if limit == 0 {
				continue
			}

			k := dt.KOfSlot[blk][0]
			diagLen := DiagLen(k, dt.M)
			numTargets := limit / diagLen

			fmt.Printf("k=%d diagLen=%d targets=%d\n", k, diagLen, numTargets)

			for t := 0; t < numTargets; t++ {
				fmt.Printf("target[%d]\n", t)
				start := t * diagLen
				end := start + diagLen

				for s := start; s < end; s++ {
					if !dt.Valid[blk][s] {
						continue
					}
					pos := dt.PosInDiag[blk][s]
					i, j := ijOf(k, pos, dt.M)
					fmt.Printf("  (i=%d,j=%d) eq≈%.4f\n", i, j, real(values[s]))
				}
			}
		}
	}

	//========================
	//=== 4) BOOTSTRAPPING ===
	//========================

	// Instantiates the evaluator

	// ckks evaluator 에는 add, mult, rotation과 같은 함수들이 있음.
	// 시간나면 member 함수 목록들을 보면서 모르는것들을 공부해나가면 좋음.


	/*
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
		//fmt.Println("Precision of values vs. ciphertext")
		//valuesTest := printDebug(params, ciphertext, valuesWant, decryptor, encoder)
		valuesTest := make([]complex128, params.MaxSlots())

		// equality test
		ciphertext = EvalEquality(ciphertext, ciphertext2, cc)

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
		valuesEQtest := printDebug(params, ciphertext, valuesTest, decryptor, encoder)

		eq_acc := 0.0
		for i := 0; i < len(valuesEQtest); i++ {
			if math.Round(real(valuesEQtest[i])) == math.Round(real(valuesTest[i])) {
				eq_acc += 1.0
			}
		}
		fmt.Println("equality accuracy : ", eq_acc/float64(len(valuesEQtest)))*/
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

func RandomStringInts(m int, alphabetSize int) []int {
	s := make([]int, m)
	for i := 0; i < m; i++ {
		s[i] = rand.Intn(alphabetSize)
	}
	return s
}

func RandomSet(n int, m int, alphabetSize int) [][]int {
	targets := make([][]int, n)
	for i := 0; i < n; i++ {
		targets[i] = RandomStringInts(m, alphabetSize)
	}
	return targets
}

func RandomSetMatch(query []int, n int, alphabetSize int) [][]int {
	m := len(query)
	targets := make([][]int, n)

	for i := 0; i < n; i++ {
		targets[i] = RandomStringInts(m, alphabetSize)
	}

	if n > 0 {
		match := make([]int, m)
		copy(match, query)
		targets[0] = match
	}

	return targets
}

func PrintDelta(dt *DeltaTable, ct *rlwe.Ciphertext, encoder *ckks.Encoder, decryptor *rlwe.Decryptor) error {
	slots := ct.Slots()
	if !ct.IsBatched {
		slots *= 2
	}

	values := make([]complex128, slots)
	if err := encoder.Decode(decryptor.DecryptNew(ct), values); err != nil {
		return err
	}

	for s := 0; s < slots; s++ {
		if !dt.Valid[0][s] {
			continue
		}
		fmt.Printf("slot=%d k=%d pos=%d eq≈%.4f\n",
			s,
			dt.KOfSlot[0][s],
			dt.PosInDiag[0][s],
			real(values[s]),
		)
	}
	return nil
}

func needRotsExpand(maxSlots, m, n int) []int {
	chunk := maxSlots / (n * m)
	if chunk < 1 {
		return nil
	}

	seen := make(map[int]bool)
	rots := make([]int, 0)

	// revSparse에서 필요한 회전
	for blk := 0; blk < n; blk++ {
		base := blk * m * chunk
		for i := 0; i < m; i++ {
			src := base + i*chunk
			dst := base + (m-1-i)*chunk
			r := src - dst
			if r != 0 && !seen[r] {
				seen[r] = true
				rots = append(rots, r)
			}
		}
	}

	// expandDiag에서 필요한 회전
	for k := 2; k <= 2*m; k++ {
		dst := 0
		diagLen := DiagLen(k, m)

		for blk := 0; blk < n; blk++ {
			for t := 0; t < diagLen; t++ {
				var i, j int
				if k <= m+1 {
					i = t
					j = diagLen - 1 - t
				} else {
					start := m - diagLen
					i = start + t
					j = m - 1 - t
				}

				srcQ := blk*m*chunk + i*chunk
				srcT := blk*m*chunk + j*chunk

				rq := srcQ - dst
				rt := srcT - dst

				if rq != 0 && !seen[rq] {
					seen[rq] = true
					rots = append(rots, rq)
				}
				if rt != 0 && !seen[rt] {
					seen[rt] = true
					rots = append(rots, rt)
				}
				dst++
			}
		}
	}

	return rots
}

func ijOf(k, pos, m int) (int, int) {
	if k <= m+1 {
		i := pos
		j := (k - 2) - pos
		return i, j
	}
	diagLen := 2*m - k + 1
	start := m - diagLen
	i := start + pos
	j := m - 1 - pos
	return i, j
}