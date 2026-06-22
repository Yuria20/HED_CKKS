package dp

import (
	"banded_dp/debug"
	dp_ideal "banded_dp/dp/dp_ideal"
	dp_sim "banded_dp/dp/dp_sim"
	he "banded_dp/he"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

func EvalEquality(cc he.Context, char_size int, hermite int, ctxt1 *rlwe.Ciphertext, ctxt2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	// debug option! : 0 (plaintext), 1 (ciphertext), 2 (validation)
	debug_opt := debug.Debug.EQ

	if debug_opt == 0 {
		return dp_ideal.EvalEQ_exp(cc, char_size, hermite, ctxt1, ctxt2)
	} else if debug_opt == 1 {
		return dp_sim.EvalEQ_exp(cc, char_size, hermite, ctxt1, ctxt2)
	} else {
		ct, err := dp_sim.EvalEQ_exp(cc, char_size, hermite, ctxt1, ctxt2)
		temp, err := dp_ideal.EvalEQ_exp(cc, char_size, hermite, ctxt1, ctxt2)
		if err != nil {
			panic(err)
		}
		debug.PrintDebug(cc.Params, ct, he.DEC(cc, temp), cc.Decryptor, cc.Encoder)

		return ct, err
	}
}

func ComputeDelta() {

}

func ComputeH() {

}

func ComputeV() {

}
