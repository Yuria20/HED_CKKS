package he

import (
	"fmt"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type Context struct {
	Params    *ckks.Parameters
	Encoder   *ckks.Encoder
	Encryptor *rlwe.Encryptor
	Decryptor *rlwe.Decryptor
	Eval      *bootstrapping.Evaluator
	Poly_eval *polynomial.Evaluator
	Max_lv    int
	logB      int
	H         int
}

func (cc Context) GenContext(logN int) {
	LogN := logN

	LogDefaultScale := 45

	q0 := []int{50}                                                         // ScaleDown & ModUp
	qiSlotsToCoeffs := []int{39, 39, 39}                                    // SlotsToCoeffs
	qiCircuitSlots := []int{45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45} // Circuit in slots
	qiLookUpTable := []int{50, 50, 50, 50, 50, 50}
	qiEvalMod := []int{60, 60, 60, 60, 60, 60, 60, 60} // EvalMod
	qiCoeffsToSlots := []int{56, 56, 56}               // CoeffsToSlots

	LogQ := append([]int{}, q0...)
	LogQ = append(LogQ, qiSlotsToCoeffs...)
	LogQ = append(LogQ, qiCircuitSlots...)

	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	var err error
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,
		LogQ:            LogQ,
		LogP:            []int{61, 61, 61, 61, 61},
		LogDefaultScale: LogDefaultScale,
		Xs:              ring.Ternary{H: 192},
	})
	if err != nil {
		panic(err)
	}

	// === Bootstrapping sub-params ===
	CoeffsToSlotsParams := dft.MatrixLiteral{
		Type:         dft.HomomorphicEncode,
		Format:       dft.RepackImagAsReal,
		LogSlots:     params.LogMaxSlots(),
		LevelQ:       params.MaxLevelQ(),
		LevelP:       params.MaxLevelP(),
		LogBSGSRatio: 1,
		Levels:       []int{1, 1, 1},
	}

	Mod1Params := mod1.ParametersLiteral{
		LevelQ:          params.MaxLevel() - CoeffsToSlotsParams.Depth(true),
		LogScale:        60,
		Mod1Type:        mod1.CosDiscrete,
		Mod1Degree:      30,
		DoubleAngle:     3,
		K:               16,
		LogMessageRatio: 10,
		Mod1InvDegree:   0,
	}

	SlotsToCoeffsParams := dft.MatrixLiteral{
		Type:         dft.HomomorphicDecode,
		LogSlots:     params.LogMaxSlots(),
		LogBSGSRatio: 1,
		LevelP:       params.MaxLevelP(),
		Levels:       []int{1, 1, 1},
	}
	SlotsToCoeffsParams.LevelQ = len(SlotsToCoeffsParams.Levels)

	btpParams := bootstrapping.Parameters{
		ResidualParameters:      params,
		BootstrappingParameters: params,
		SlotsToCoeffsParameters: SlotsToCoeffsParams,
		Mod1ParametersLiteral:   Mod1Params,
		CoeffsToSlotsParameters: CoeffsToSlotsParams,
		EphemeralSecretWeight:   32,
		CircuitOrder:            bootstrapping.DecodeThenModUp,
	}
	fmt.Printf("Bootstrapping parameters: logN=%d, logSlots=%d, H(%d; %d), sigma=%f, logQP=%f, levels=%d, scale=2^%d\n",
		btpParams.BootstrappingParameters.LogN(),
		btpParams.BootstrappingParameters.LogMaxSlots(),
		btpParams.BootstrappingParameters.XsHammingWeight(),
		btpParams.EphemeralSecretWeight,
		btpParams.BootstrappingParameters.Xe(),
		btpParams.BootstrappingParameters.LogQP(),
		btpParams.BootstrappingParameters.QCount(),
		btpParams.BootstrappingParameters.LogDefaultScale())

	//fmt.Printf("Maximum Level : %d, Minimum Level : %d\n", MaxLv)
	// === Keys + basic objects ===
	kgen := rlwe.NewKeyGenerator(params)
	sk, pk := kgen.GenKeyPairNew()

	Encoder := ckks.NewEncoder(params)
	Decryptor := rlwe.NewDecryptor(params, sk)
	Encryptor := rlwe.NewEncryptor(params, pk)

	fmt.Println("Generating bootstrapping evaluation keys...")
	evk, _, err := btpParams.GenEvaluationKeys(sk)
	if err != nil {
		panic(err)
	}

	Eval, err := bootstrapping.NewEvaluator(btpParams, evk)
	if err != nil {
		panic(err)
	}

	cc.Params = &params
	cc.Encoder = Encoder
	cc.Encryptor = Encryptor
	cc.Decryptor = Decryptor
	cc.Max_lv = len(LogQ) - 1
	cc.Eval = Eval
}
