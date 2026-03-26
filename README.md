# banded_HED

### Go version 
`go version go1.23.3 linux/amd64`
제가 알기론 버전 달라도 가능한걸로 아는데요. 만약 있으면 위의 버전으로 하시는게 좋을 것 같습니다.

### Project folder 
`./lattigo/examples/singleparty/ckks_bootstrapping/slim`
가끔 내부 module 수정해야할 때가 있어서, lattigo 폴더 통째로 넣어놨습니다.
해당 폴더에서 다음 명령어 실행하면, 코드가 실행되도록 되있습니다. 
`go run . ` 

### To Do Lists 
- Equality checking function  
- To generate an EQ table of size $m \times m$, rearrange the positions of the $1 \times m$ string.
- Interpolation polynomial for $H, V$
- Algorithm for update Horizon, Vertical Table
