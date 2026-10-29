# banded_HED

## 환경

- Go version: `go1.23.3 linux/amd64`
- Project folder: `./lattigo/examples/singleparty/ckks_bootstrapping/slim`

가끔 내부 module을 수정해야 할 때가 있어서, lattigo 폴더를 통째로 넣어놨습니다.
해당 폴더에서 아래 명령어를 실행하면 코드가 돌아갑니다.

```bash
go run .
```

---

## 1. 사소한 변경사항

1. 가독성을 위해 문자 공간의 크기를 나타내던 `char_size`를 Context 내부에 할당하고,
   `char_size`를 사용하던 모든 함수들을 변경하였습니다.

2. 이전 branch에서 만든 Equality 연산을 활용하여 `dp/dp_module.go` 에 `EqTable` 함수를
   작성하였습니다.

```go
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
```