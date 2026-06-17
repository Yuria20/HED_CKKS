package he

import (
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/bootstrapping"
	"github.com/tuneinsight/lattigo/v6/circuits/ckks/polynomial"
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type Context struct {
	Params    *ckks.Parameters
	Encoder   *ckks.Encoder
	Encryptor *rlwe.Encryptor
	Decryptor *rlwe.Decryptor
	Eval      *ckks.Evaluator
	Btp       *bootstrapping.Evaluator
	Poly_eval *polynomial.Evaluator
	Max_lv    int
}
