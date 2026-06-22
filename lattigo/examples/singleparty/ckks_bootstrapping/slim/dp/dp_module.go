package dp

import (
	"banded_dp/debug"
	dp_ideal "banded_dp/dp/dp_ideal"
	dp_sim "banded_dp/dp/dp_sim"
	he "banded_dp/he"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

type EqTable struct {
	Deltas          [][]*rlwe.Ciphertext
	CurrAntiDiagIdx int
}

func GenEqTable(cc he.Context, strct1 he.StrCipher, strct2 he.StrCipher) EqTable {

	var E EqTable
	E.Deltas = make([][]*rlwe.Ciphertext, strct1.Length)
	E.CurrAntiDiagIdx = 0

	for i := 0; i < strct1.Length; i++ {
		E.Deltas[i] = make([]*rlwe.Ciphertext, strct2.Length)
		for j := 0; j < strct2.Length; j++ {
			E.Deltas[i][j], _ = EvalEQ_exp(cc, 0, strct1.Cipher[i], strct2.Cipher[j])
		}
	}

	return E
}

func EvalEQ_exp(cc he.Context, hermite int, ctxt1 *rlwe.Ciphertext, ctxt2 *rlwe.Ciphertext) (*rlwe.Ciphertext, error) {
	// debug option! : 0 (ciphertext), 1 (plaintext), 2 (validation)
	debug_opt := debug.Debug.EQ

	if debug_opt == 0 {
		return dp_sim.EvalEQ_exp(cc, hermite, ctxt1, ctxt2)
	} else if debug_opt == 1 {
		return dp_ideal.EvalEQ_exp(cc, hermite, ctxt1, ctxt2)
	} else {
		ct, err := dp_sim.EvalEQ_exp(cc, hermite, ctxt1, ctxt2)
		temp, err := dp_ideal.EvalEQ_exp(cc, hermite, ctxt1, ctxt2)
		if err != nil {
			panic(err)
		}
		debug.PrintDebug(cc.Params, ct, he.DEC(cc, temp), cc.Decryptor, cc.Encoder)

		return ct, err
	}
}
