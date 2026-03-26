## `dp_full.py` 결과 해석

`dp_full.py`를 실행하면 전체 anti-diagonal DP 결과가 `.npz` 파일로 저장된다.  
이 파일을 열어서 출력하면 다음과 같은 항목들이 보인다.

### 예시 출력

```text
pid = r0
L = 512
mode = full
i0 = [0 0 0 ... 510 511 512]
Hsym = [...]
Vsym = [...]
dist_true = 20
H = [...]
V = [...]
dist = 20
dist_sym = 20
```

### 각 항목의 의미

- `pid`  
  입력 파일에서 읽은 pair ID이다.

- `L`  
  문자열 길이이다.  
  현재 구현에서는 두 문자열 길이가 같다고 가정한다.

- `mode`  
  계산 모드이다.  
  `full`은 band를 자르지 않고 전체 anti-diagonal DP를 계산했다는 뜻이다.

- `i0`  
  각 anti-diagonal의 시작 row index이다.  
  anti-diagonal `k = i + j`에 대해, `i0[k]`는 그 대각선에서 저장된 첫 번째 row 위치를 의미한다.

### 심볼 도메인 결과

- `Hsym`, `Vsym`  
  심볼 매핑된 DP 값들을 anti-diagonal 단위로 저장한 결과이다.

  이 구현에서는 다음과 같은 심볼 매핑을 사용한다.

  - `g(-1) = 0.5 + 0j`
  - `g(0) = 0 + 0j`
  - `g(1) = 0 + 0.5j`

  따라서 다음과 같이 해석할 수 있다.

  - `0.5 + 0j` → 정수값 `-1`
  - `0 + 0j` → 정수값 `0`
  - `0 + 0.5j` → 정수값 `1`

  즉 `Hsym`, `Vsym`는 아직 정수로 복원되기 전의 복소수 심볼 상태라고 보면 된다.

### 정수 도메인 결과

- `H`, `V`  
  `Hsym`, `Vsym`에 inverse mapping `g_inv`를 적용해서 얻은 정수값 결과이다.

  각 anti-diagonal의 값은 `{-1, 0, 1}` 중 하나로 양자화된다.

### 거리 관련 항목

- `dist_true`  
  일반적인 full DP로 계산한 정확한 edit distance이다.  
  즉 기준이 되는 정답이다.

- `dist`  
  심볼 기반 recurrence로 계산한 뒤, 마지막에 정수값으로 복원해서 얻은 edit distance이다.

- `dist_sym`  
  심볼 기반 계산 결과를 확인하기 위한 값이다.  
  현재 full mode에서는 `dist`와 같은 의미로 보면 된다.

### 올바른 결과인지 확인하는 방법

`dp_full.py`가 정상적으로 동작했다면 다음이 같아야 한다.

```text
dist_true == dist == dist_sym
```

예를 들어

```text
dist_true = 20
dist = 20
dist_sym = 20
```

라고 나오면 다음을 의미한다.

1. 일반 full DP의 정답이 20이다.
2. 심볼 매핑 기반 recurrence로 복원한 값도 20이다.
3. 따라서 full 구현이 내부적으로 일관되게 동작하고 있다.

### 참고

- `dp_full.py`는 정답 검증용 baseline 코드이다.
- `dp_patched_layout.py` 결과는 이 baseline과 비교해서 해석해야 한다.
- 특히 `dp_patched_layout.py`를 `w = L + 1`처럼 사실상 full width로 실행했을 때는, 복원된 거리값이 `dp_full.py`의 `dist_true`와 같아야 한다.
