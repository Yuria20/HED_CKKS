package he

import (
	"flag"
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"github.com/tuneinsight/lattigo/v6/core/rlwe"
	"github.com/tuneinsight/lattigo/v6/schemes/ckks"
	"github.com/tuneinsight/lattigo/v6/utils/sampling"
)

var flagShort = flag.Bool("short", false, "run the example with a smaller and insecure ring degree.")

func Test_diBtp(t *testing.T) {

	flag.Parse()

	// Default LogN, which with the following defined parameters
	// provides a security of 128-bit.
	LogN := 16

	var cc Context
	cc.GenContext(LogN, 1, 4, 1)

	params := cc.Params
	encoder := cc.Encoder
	encryptor := cc.Encryptor
	decryptor := cc.Decryptor
	B := int(1 << cc.LogB)

	// Generate a random plaintext with values uniformly distributed in [-1, 1] for the real and imaginary part.
	valuesWant := make([]complex128, params.MaxSlots())
	for i := range valuesWant {
		valuesWant[i] = complex(float64(rand.Intn(B)), 0)
	}

	// We encrypt at level 0
	plaintext := ckks.NewPlaintext(params, cc.Min_lv+1)
	if err := encoder.Encode(valuesWant, plaintext); err != nil {
		panic(err)
	}

	// Encrypt
	ciphertext, err := encryptor.EncryptNew(plaintext)
	if err != nil {
		panic(err)
	}

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of values vs. ciphertext")
	valuesTest := printDebug_with_Noise(params, ciphertext, valuesWant, decryptor, encoder)

	fmt.Println("Bootstrapping...")

	start := time.Now()
	ciphertext, err = DiBootstrapMany(cc, ciphertext, func(x int) int { return x })
	must(err)
	btp_time := time.Since(start)

	fmt.Printf("Done... ")
	fmt.Println(btp_time)

	//==================
	//=== 5) DECRYPT ===
	//==================

	// Decrypt, print and compare with the plaintext values
	fmt.Println()
	fmt.Println("Precision of ciphertext vs. Bootstrap(ciphertext)")
	PrintDebug(params, ciphertext, valuesTest, decryptor, encoder)
}

func printDebug_with_Noise(params ckks.Parameters, ciphertext *rlwe.Ciphertext, valuesWant []complex128, decryptor *rlwe.Decryptor, encoder *ckks.Encoder) (valuesTest []complex128) {

	slots := ciphertext.Slots()

	if !ciphertext.IsBatched {
		slots *= 2
	}

	valuesTest = make([]complex128, slots)

	if err := encoder.Decode(decryptor.DecryptNew(ciphertext), valuesTest); err != nil {
		panic(err)
	}

	// insert noise
	for i := range valuesWant {
		valuesTest[i] = valuesTest[i] + sampling.RandComplex128(-1e-4, 1e-4)
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
