package he

import (
	"fmt"
	"math/big"

	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/dft"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/mod1"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/ring"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type Context struct {
	Params    ckks.Parameters
	Encoder   *ckks.Encoder
	Encryptor *rlwe.Encryptor
	Decryptor *rlwe.Decryptor
	Eval      *bootstrapping.Evaluator
	PolyEval  *polynomial.Evaluator
	Min_lv    int
	Avail_lv  int
	Exp_lv    int
	Max_lv    int
	LogB      int
	H_btp     int
	C_size    int
}

func (cc *Context) GenContext(logN int, Avail_lv int, LogB int, H int) {
	LogN := logN

	LogDefaultScale := 50
	sz := LogDefaultScale

	q0 := []int{sz}                      // ScaleDown & ModUp
	qiSlotsToCoeffs := []int{39, 39, 39} // SlotsToCoeffs
	qiCircuitSlots := []int{}            // Circuit in slots
	for i := 0; i < Avail_lv; i++ {
		qiCircuitSlots = append(qiCircuitSlots, []int{sz}...)
	}
	qiLookUpTable := []int{}
	for i := 0; i < LogB+H; i++ {
		qiLookUpTable = append(qiLookUpTable, []int{sz}...)
	}
	qiEvalMod := []int{sz, sz, sz, sz, sz, sz, sz, sz} // EvalMod
	qiCoeffsToSlots := []int{42, 42, 42}               // CoeffsToSlots

	/// print
	arrays := []struct {
		name string
		vals []int
	}{
		{"q0", q0},
		{"S2C", qiSlotsToCoeffs},
		{"Circuit", qiCircuitSlots},
		{"LUT", qiLookUpTable},
		{"EvalExp", qiEvalMod},
		{"C2S", qiCoeffsToSlots},
	}

	for _, a := range arrays {
		fmt.Printf("%s : %v\n", a.name, a.vals)
	}

	LogQ := append([]int{}, q0...)
	LogQ = append(LogQ, qiSlotsToCoeffs...)
	min_lv := len(LogQ) - 1
	LogQ = append(LogQ, qiCircuitSlots...)
	max_lv := len(LogQ) - 1
	LogQ = append(LogQ, qiLookUpTable...)
	exp_lv := len(LogQ) - 1
	LogQ = append(LogQ, qiEvalMod...)
	LogQ = append(LogQ, qiCoeffsToSlots...)

	var err error
	params, err := ckks.NewParametersFromLiteral(ckks.ParametersLiteral{
		LogN:            LogN,
		LogQ:            LogQ,
		LogP:            []int{sz + 1, sz + 1, sz + 1, sz + 1, sz + 1},
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
		LogScale:        sz,
		Mod1Type:        mod1.CosDiscrete,
		Mod1Degree:      30,
		DoubleAngle:     3,
		K:               16,
		LogMessageRatio: 0,
		Mod1InvDegree:   0,
	}

	B := 1 << LogB
	S2C_scale := *big.NewFloat(float64(1.0) / float64(B))
	SlotsToCoeffsParams := dft.MatrixLiteral{
		Type:         dft.HomomorphicDecode,
		LogSlots:     params.LogMaxSlots(),
		LogBSGSRatio: 1,
		LevelP:       params.MaxLevelP(),
		Levels:       []int{1, 1, 1},
		Scaling:      &S2C_scale,
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

	polyEval := polynomial.NewEvaluator(params, Eval.Evaluator)

	cc.Params = params
	cc.Encoder = Encoder
	cc.Encryptor = Encryptor
	cc.Decryptor = Decryptor
	cc.Min_lv = min_lv
	cc.Avail_lv = Avail_lv
	cc.Exp_lv = exp_lv
	cc.Max_lv = max_lv
	cc.Eval = Eval
	cc.LogB = LogB
	cc.H_btp = H
	cc.PolyEval = polyEval
}
