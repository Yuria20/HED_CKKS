# banded_HED

### Go version 
`go version go1.23.3 linux/amd64`

### Project folder 
`./lattigo/examples/singleparty/ckks_bootstrapping/slim`
가끔 내부 module 수정해야할 때가 있어서, lattigo 폴더 통째로 넣어놨습니다.
해당 폴더에서 다음 명령어 실행하면, 코드가 실행되도록 되있습니다. 
`go run . ` 


### Setup 
동형암호 코드를 작성하기 전에, 작업을 어떤 구조로 진행할지 먼저 정리한다.
 
이번 작업은 **edit distance의 동형 알고리즘**이다.
이때 **Dynamic Programming (DP) Table**을 사용해야 하므로, 이와 관련된 모듈들은 `dp_module`에 넣는다.
 
그리고 동형암호, 즉 **HE와 관련된 모듈들은 `he_module`**에 넣는다.
 
마지막으로 `ideal`과 `simulator`는 각각 **평문 알고리즘**과 **동형 알고리즘**을 실행하기 위한 폴더이다.
 
---
 
## 폴더 구조
 
| 폴더 | 역할 |
|------|------|
| `dp_module` | DP Table 관련 모듈 |
| `he_module` | 동형암호(HE) 관련 모듈 |
| `ideal` | 평문 알고리즘 실행 |
| `simulator` | 동형 알고리즘 실행 |
 
---
 
## 디버깅 방식
 
우리는 디버깅의 편의를 위해 **`ideal`과 `simulator`를 같이 돌릴 예정**이다.
 
같은 시점에 대해
- 암호문을 복호화한 평문값과
- 평문 알고리즘의 중간 평문값을
서로 비교한다.
 
필자의 경우 원래 귀찮아서 이런 거 안 하는데...
정신건강상 이게 **58000% 이득**이긴 하다.
