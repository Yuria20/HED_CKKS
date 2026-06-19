package dp_sim

import (
	he "banded_dp/he"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type StrCipher struct {
	Length int
	Cipher []*rlwe.Ciphertext
}

func (strct *StrCipher) StrEncrypt(cc he.Context, strpt [][]int) {

	params := cc.Params
	encoder := cc.Encoder
	encryptor := cc.Encryptor
	max_lv := cc.Max_lv

	// strpt : N x m
	m := len(strpt[0])
	N := len(strpt)

	var err error
	strct.Cipher = make([]*rlwe.Ciphertext, m)

	for i := 0; i < m; i++ {
		values := make([]complex128, params.MaxSlots())
		for j := 0; j < N; j++ {
			values[j] = complex(float64(strpt[j][i]), 0)
		}

		plaintext := ckks.NewPlaintext(params, max_lv)
		if err := encoder.Encode(values, plaintext); err != nil {
			panic(err)
		}

		strct.Cipher[i], err = encryptor.EncryptNew(plaintext)
		if err != nil {
			panic(err)
		}

	}

	return
}

func (strct *StrCipher) StrDecrypt(cc he.Context) [][]int {

	params := cc.Params
	encoder := cc.Encoder
	decryptor := cc.Decryptor

	var strpt [][]int
	strpt = make([][]int, params.MaxSlots())
	for i := 0; i < params.MaxSlots(); i++ {
		strpt[i] = make([]int, strct.Length)
	}
	// strpt : N x m
	m := len(strpt[0])
	N := len(strpt)

	for i := 0; i < m; i++ {
		values := make([]complex128, params.MaxSlots())
		plaintext := decryptor.DecryptNew(strct.Cipher[i])
		encoder.Decode(plaintext, values)

		for j := 0; j < N; j++ {
			strpt[j][i] = int(math.Round(real(values[j])))
		}
	}

	return strpt
}
