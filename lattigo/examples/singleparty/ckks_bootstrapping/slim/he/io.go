package he

import (
	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

func ENC(cc Context, values []complex128, level int) (*rlwe.Ciphertext, error) {

	pt := ckks.NewPlaintext(cc.Params, level)
	if err := cc.Encoder.Encode(values, pt); err != nil {
		return nil, err
	}
	return cc.Encryptor.EncryptNew(pt)
}

func DEC(cc Context, ctxt *rlwe.Ciphertext) []complex128 {
	var err error
	values := make([]complex128, cc.Params.MaxSlots())
	plaintext := cc.Decryptor.DecryptNew(ctxt)
	if err = cc.Encoder.Decode(plaintext, values); err != nil {
		return nil
	}

	return InverseInterLeave(values, cc.M, cc.Slot_size)
}
