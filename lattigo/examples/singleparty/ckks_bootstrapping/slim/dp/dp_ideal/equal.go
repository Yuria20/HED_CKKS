package dp_ideal

import (
	he "banded_dp/he"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

func EvalEQ_exp(cc he.Context, char_size int, hermite int, ctxt1 *rlwe.Ciphertext, ctxt2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {

	in_lv := min(ctxt1.Level(), ctxt2.Level())
	out_lv := in_lv - (int(math.Ceil(math.Log2(float64(char_size)))) + 1)

	values1 := he.DEC(cc, ctxt1)
	values2 := he.DEC(cc, ctxt2)

	values_want := make([]complex128, cc.Params.MaxSlots())

	for i := 0; i < len(values_want); i++ {
		// round
		val1_int := he.RootToByte(values1[i], char_size)
		val2_int := he.RootToByte(values2[i], char_size)

		// eq
		val_int := val1_int - val2_int
		if val_int == 0 {
			values_want[i] = complex(1, 0)
		} else {
			values_want[i] = complex(0, 0)
		}
	}

	return he.ENC(cc, values_want, out_lv)
}
