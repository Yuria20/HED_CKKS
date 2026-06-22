package debug

import (
	"fmt"
	"math"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
)

type DebugBox struct {
	EQ int
}

// global
var Debug = DebugBox{
	// 0 : cipher, 1 : plain, 2 : validation
	EQ: 2,
}

func PrintDebug(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	max_L2_err := 100.0
	avg_L2_err := 0.0

	for i := 0; i < params.MaxSlots(); i++ {
		// compute real err
		real_round := real(valuesWant[i])
		real := real(valuesTest[i])
		real_err := math.Abs(real_round - real)

		// compute imag err
		imag_round := imag(valuesWant[i])
		imag := imag(valuesTest[i])
		imag_err := math.Abs(imag_round - imag)

		L2_err := math.Sqrt(real_err*real_err + imag_err*imag_err)
		log_L2_err := -1 * math.Log2(L2_err)

		avg_L2_err += log_L2_err
		max_L2_err = min(max_L2_err, log_L2_err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %10.14f %10.14f %10.14f %10.14f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %10.14f %10.14f %10.14f %10.14f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println("Average Noise (L2) : ", -1*avg_L2_err/float64(len(valuesTest)))
	fmt.Println("Maximum Noise (L2) : ", -1*max_L2_err)
	//precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	//fmt.Println(precStats.String())
	fmt.Println()

	return
}

func PrintNoise(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	max_L2_err := 100.0
	avg_L2_err := 0.0

	for i := 0; i < params.MaxSlots(); i++ {
		// compute real err
		real_round := math.Round(real(valuesTest[i]))
		real := real(valuesTest[i])
		real_err := math.Abs(real_round - real)

		// compute imag err
		imag_round := math.Round(imag(valuesTest[i]))
		imag := imag(valuesTest[i])
		imag_err := math.Abs(imag_round - imag)

		L2_err := math.Sqrt(real_err*real_err + imag_err*imag_err)
		log_L2_err := -1 * math.Log2(L2_err)

		avg_L2_err += log_L2_err
		max_L2_err = min(max_L2_err, log_L2_err)
	}

	fmt.Println()
	fmt.Printf("Level: %d (logQ = %d)\n", ciphertext.Level(), params.LogQLvl(ciphertext.Level()))

	fmt.Printf("Scale: 2^%f\n", math.Log2(ciphertext.Scale.Float64()))
	fmt.Printf("ValuesTest: %10.14f %10.14f %10.14f %10.14f...\n", valuesTest[0], valuesTest[1], valuesTest[2], valuesTest[3])
	fmt.Printf("ValuesWant: %10.14f %10.14f %10.14f %10.14f...\n", valuesWant[0], valuesWant[1], valuesWant[2], valuesWant[3])
	fmt.Println("Average Noise (L2) : ", -1*avg_L2_err/float64(len(valuesTest)))
	fmt.Println("Maximum Noise (L2) : ", -1*max_L2_err)
	//precStats := ckks.GetPrecisionStats(params, encoder, nil, valuesWant, valuesTest, 0, false)

	//fmt.Println(precStats.String())
	fmt.Println()

	return
}
