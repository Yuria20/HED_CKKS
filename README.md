# banded_HED

### Go version

```
go version go1.23.3 linux/amd64
```

### Project folder

```
./lattigo/examples/singleparty/ckks_bootstrapping/slim
```

가끔 내부 module을 수정해야 할 때가 있어서, lattigo 폴더를 통째로 넣어놨습니다.
해당 폴더에서 아래 명령어를 실행하면 코드가 돌아갑니다.

```bash
go run .
```

---

## Main1 : Enc / Dec / Print 함수 만들기

### 1. 폴더 구조 변경

```
slim/
├── dp/
│   ├── dp_ideal/
│   └── dp_sim/
│       ├── DPtable.go
│       ├── equal.go
│       └── pack.go
├── he/
│   ├── btp.go
│   ├── context.go
│   └── poly.go
├── test/
│   ├── dp/pack_test.go
│   └── he/he_test.go
├── util/
│   ├── debug.go
│   ├── print.go
│   └── sampling.go
├── go.mod
├── go.sum
├── main_test.go
└── main.go
```

1. 이름이 맛이 없어서 `dp_module_ideal`, `dp_module_sim`, `he_module` 을 각각 `dp_ideal`, `dp_sim`, `he` 로 변경함.
2. `dp_ideal` 과 `dp_sim` 을 `dp` 폴더에 넣었음. 이후 브랜치에서 설명할 예정.
3. `he/util_poly.go` 를 `he/util.go` 로 이름 변경.
4. `he` 에 있던 `debug.go` 를 `util` 폴더에 집어넣었고, `he_test.go` 를 `test/he` 에 집어넣었음.
5. `pack.go`, `print.go`, `sampling.go` 는 이번 브랜치에서 추가되는 것들.

### 2. util/sampling.go

알고리즘을 짜기 전에 값을 암호화해야 하니까 샘플링하는 함수를 작성했음.
기존에 작성되었던 것들이 있어서, 주석이랑 이름 등만 변경하였음.

```go
func SampleStringIntMany(N int, m int, sz int) [][]int {
	valuesArr := make([][]int, N)
	for i := 0; i < N; i++ {
		valuesArr[i] = SampleStringInts(m, sz)
	}
	return valuesArr
}
```

실제로 사용되는 함수는 위에 있는 거로, 길이 `m` 인 string을 `N` 개만큼 sampling 하며 각 문자는 `[0, sz)` 사이에서 sampling 됨.

### 3. dp/dp_sim/pack.go

이제 sampling 된 string 들을 암호화하고 복호화할 함수를 `pack.go` 에 만들 예정.

```go
type StrCipher struct {
	Length int
	Cipher []*rlwe.Ciphertext
}
```

일단 Go는 클래스가 없으므로 깔쌈하게 구조체 선언. 이름은 아무렇게나 지었음.
`StrCipher` 안의 `Cipher` 에서 암호문들을 저장할 거임.

자 이제 `StrEncrypt` 와 `StrDecrypt` 라는 함수를 `StrCipher` 안의 멤버 함수로 넣어주고 끝내자.
얘네는 `[][]int` 배열과 `[]*rlwe.Ciphertext` 로 암호화와 복호화를 수행할 거임.

이 프로젝트의 경우, 길이 `m` 의 string을 `N` 개 처리할 건데,
`string[N][m]` 이 이렇게 저장되어 있고,
`N` = 암호문의 slot 갯수가 되게끔 setting 해서
`m` 개의 암호문에 `N` 개의 string이 담길 예정.

> 즉, `i` 번째 암호문의 `j` 번째 slot에는 `string[j][i]` 의 정수값이 저장됨.

### 4. util/print.go

이제 string들을 print하는 함수를 하나 만들어준다.
string 갯수가 워낙 많기 때문에, `0`, `1` 그리고 `N-1` 번지만 출력해준다.
`tag` 를 달아서 `"uint8"` 이면 정수 그대로 출력하고, `"str"` 이면 string 값으로 출력된다.

```go
func PrintString(strpt [][]int, tag string) {

	// strpt : N x m
	m := len(strpt[0])
	N := len(strpt)
	for i := 0; i < 2; i++ {
		fmt.Printf("str[%d] : ", i)
		for j := 0; j < m; j++ {
			if tag == "uint8" {
				fmt.Print(strpt[i][j])
				fmt.Print(" ")
			} else if tag == "str" {
				fmt.Print(string(strpt[i][j]))
			}
		}
		fmt.Println()
	}
	fmt.Println("......")
	for i := N - 1; i < N; i++ {
		fmt.Printf("str[%d] : ", i)
		for j := 0; j < m; j++ {
			if tag == "uint8" {
				fmt.Print(strpt[i][j])
				fmt.Print(" ")
			} else if tag == "str" {
				fmt.Print(string(strpt[i][j]))
			}
		}
		fmt.Println()
	}
}
```

실제 디버깅할 때는 복호화를 해서 값을 봐야 하는 경우가 많기 때문에,
복호화 + 출력까지 되는 함수를 추가로 만들어준다.

```go
func StrDecryptAndPrint(cc he.Context, strct dp.StrCipher, tag string) {
	PrintString(strct.StrDecrypt(cc), tag)
}
```

### 5. test 하기

```go
LogN := 13

var cc he.Context
fmt.Print("Generate Context... ")
cc.GenContext(LogN, 10, 4, 1)
fmt.Println("Done!")

str_len := 8
char_size := 127

var pack_ct dp.StrCipher
pack_ct.Length = str_len

strArr_want := util.SampleStringIntMany(cc.Params.MaxSlots(), str_len, char_size)
util.PrintString(strArr_want, "uint8")

pack_ct.StrEncrypt(cc, strArr_want)

util.StrDecryptAndPrint(cc, pack_ct, "uint8")

strArr_test := pack_ct.StrDecrypt(cc)

equality_test(strArr_want, strArr_test)
```

일단 `logN = 16` 이 표준 파라미터지만, 기다리다 지치므로 test는 `logN = 13` 으로 실행해준다.
그리고 context 생성해주고, 문자열 길이 (`str_len = m`), 문자 공간의 크기 (`char_size`) 설정해 둔 다음에
만들어진 평문과 암복호화된 평문이 서로 일치하는지 test해주면 끝.

전체 코드는 `test/dp/pack_test.go` 에 있다.

### 6. Appendix

원래 이 정도까지 폴더 나누고, 함수 쪼개고 test 함수도 따로 다 빼서 하진 않는데,
이 기록은 FHE 딸수들이 보는 것이기 때문에 FM으로 진행한다.
여러분들도 익숙해지면 가라를 섞어서 코드 짜길 바란다...
 