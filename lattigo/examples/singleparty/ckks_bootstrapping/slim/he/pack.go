package he

import (
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

type StrCipher struct {
	Length int
	Cipher []*rlwe.Ciphertext
}

// int => root of unity
func ByteToRoot(x int, mod int) complex128 {
	return complex(
		math.Cos(2*math.Pi*(float64(x)/float64(mod))),
		math.Sin(2*math.Pi*(float64(x)/float64(mod))),
	)
}

// root of unity => int
func RootToByte(z complex128, mod int) int {
	theta := math.Atan2(imag(z), real(z)) // range: (-pi, pi]

	if theta < 0 {
		theta += 2 * math.Pi
	}

	x := int(math.Round(theta / (2 * math.Pi) * float64(mod)))

	// round 때문에 x == mod가 나올 수 있음
	x = ((x % mod) + mod) % mod

	return x
}

func (strct *StrCipher) StrEncrypt(cc Context, sz int, level int, strpt [][]int) {

	params := cc.Params

	// strpt : N x m
	m := len(strpt[0])
	N := len(strpt)

	var err error
	strct.Cipher = make([]*rlwe.Ciphertext, m)

	for i := 0; i < m; i++ {
		values := make([]complex128, params.MaxSlots())
		for j := 0; j < N; j++ {
			values[j] = ByteToRoot(strpt[j][i], sz)
		}

		strct.Cipher[i], err = ENC(cc, values, level)
		if err != nil {
			panic(err)
		}

	}

	return
}

func (strct *StrCipher) StrDecrypt(cc Context) [][]int {

	params := cc.Params

	var strpt [][]int
	strpt = make([][]int, params.MaxSlots())
	for i := 0; i < params.MaxSlots(); i++ {
		strpt[i] = make([]int, strct.Length)
	}
	// strpt : N x m
	m := len(strpt[0])
	N := len(strpt)

	for i := 0; i < m; i++ {
		values := DEC(cc, strct.Cipher[i])

		for j := 0; j < N; j++ {
			strpt[j][i] = int(math.Round(real(values[j])))
		}
	}

	return strpt
}
