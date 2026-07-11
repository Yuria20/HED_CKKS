package he

import (
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
)

type StrCipher struct {
	Length int
	Cipher *rlwe.Ciphertext
}

// int => root of unity
func ByteToRoot(x int, mod int) complex128 {
	return complex(
		math.Cos(2*math.Pi*(float64(x)/float64(mod))),
		math.Sin(2*math.Pi*(float64(x)/float64(mod))),
	)
}

// root ofe unity => int
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

func (strct *StrCipher) StrEncrypt(cc Context, level int, strpt [][]int) {

	params := cc.Params
	sz := cc.C_size

	// strpt : N x m
	str_len := cc.M
	b_sz := cc.Batch_size
	var err error

	values := make([]complex128, params.MaxSlots())
	for i := 0; i < b_sz; i++ {
		for j := 0; j < str_len; j++ {
			//values[i*str_len+j] = complex(float64(strpt[i][j]), 0.0) //ByteToRoot(strpt[i][j], sz)
			values[i*str_len+j] = ByteToRoot(strpt[i][j], sz)
		}
	}
	values = InterLeave(values, str_len, cc.Slot_size)
	// interleaving

	strct.Cipher, err = ENC(cc, values, level)
	if err != nil {
		panic(err)
	}

	return
}

func (strct *StrCipher) StrDecrypt(cc Context) [][]int {

	//params := cc.Params

	// strpt : N x m
	str_len := cc.M
	b_sz := cc.Batch_size

	var strpt [][]int
	strpt = make([][]int, b_sz)
	for i := 0; i < b_sz; i++ {
		strpt[i] = make([]int, str_len)
	}

	values := DEC(cc, strct.Cipher)
	values = InverseInterLeave(values, str_len, cc.Slot_size)

	for i := 0; i < b_sz; i++ {
		for j := 0; j < str_len; j++ {
			strpt[i][j] = int(math.Round(real(values[i*str_len+j])))
		}
	}

	return strpt
}

func InterLeave(vec []complex128, n int, N int) []complex128 {
	batch := N / n
	twistedVec := make([]complex128, N)

	for i := 0; i < n; i++ {
		for j := 0; j < batch; j++ {
			twistedVec[i*batch+j] = vec[j*n+i]
		}
	}

	return twistedVec
}

func InverseInterLeave(vec []complex128, n int, N int) []complex128 {
	batch := N / n
	InvtwistedVec := make([]complex128, N)

	for i := 0; i < n; i++ {
		for j := 0; j < batch; j++ {
			InvtwistedVec[j*n+i] = vec[i*batch+j]
		}
	}

	return InvtwistedVec
}
