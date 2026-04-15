from typing import List
import rlwe
import numpy as np

def ED_ideal(a: List[int], b: List[int]) -> int:
    n, m = len(a), len(b)

    dp = [[0] * (m + 1) for _ in range(n + 1)]

    # 초기 조건
    for i in range(n + 1):
        dp[i][0] = i
    for j in range(m + 1):
        dp[0][j] = j

    # DP 계산
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0 if a[i-1] == b[j-1] else 1

            dp[i][j] = min(
                dp[i-1][j] + 1,      # 삭제
                dp[i][j-1] + 1,      # 삽입
                dp[i-1][j-1] + cost  # 치환
            )

    return dp[n][m]

from typing import List

def myers_with_tables(a: List[int], b: List[int]):
    m = len(a)
    mask = (1 << m) - 1

    # Eq 생성
    Eq = {}
    for val in set(a):
        Eq[val] = 0
    for i, val in enumerate(a):
        Eq[val] |= 1 << i

    Pv = mask
    Mv = 0
    score = m

    # 테이블 저장
    H_table = []  # (Ph, Mh)
    V_table = []  # (Pv, Mv)

    for val in b:
        Eq_v = Eq.get(val, 0)

        Xv = Eq_v | Mv
        Xh = (((Eq_v & Pv) + Pv) ^ Pv) | Eq_v

        Ph = Mv | ~(Xh | Pv)
        Mh = Pv & Xh

        # score 업데이트
        if Ph & (1 << (m - 1)):
            score += 1
        elif Mh & (1 << (m - 1)):
            score -= 1

        # shift
        Ph = (Ph << 1) & mask
        Mh = (Mh << 1) & mask

        # update
        Pv = (Mh | ~(Xv | Ph)) & mask
        Mv = (Ph & Xv) & mask

        # 저장 (bit → list 형태)
        H_table.append((
            [(Ph >> i) & 1 for i in range(m)],
            [(Mh >> i) & 1 for i in range(m)]
        ))

        V_table.append((
            [(Pv >> i) & 1 for i in range(m)],
            [(Mv >> i) & 1 for i in range(m)]
        ))

    return H_table, V_table, score

def ED_sim(cc : rlwe, a : List[int],b : List[int]) :
    debug = False
    m = cc.num_slots

    A = cc.SetCiphertext(a)
    B = cc.SetCiphertext(b)

    # B : automorphism invariant
    arrA = re_ordering(cc,A)
    
    if debug : 
        print("re-ordering A")
        print(arrA)

    # make Equality Table    
    E = np.zeros((m,m),dtype=np.float64)
    for i in range(m) :
        E[i] = cc.eq(B,arrA[i])
        E[i] = -1*E[i]
    
    if debug :
        print("Equaltiy Table")
        print(E)

    # Compute DP table
    H = np.zeros(m,dtype=np.float64)
    V = np.zeros(m,dtype=np.float64)
    H[0] = 1.0
    V[0] = 1.0
    acc = np.zeros(m,dtype=np.float64)

    for i in range(1,m) :
        mask = make_mask_foward(cc,i)
        E_m = cc.mul(E[i-1],mask)
        
        min = cc.min(E_m,H,V)
        T = 1-V+min
        V = 1-H+min
        T = cc.mul(T,mask)
        V = cc.mul(V,mask)
        
        H = cc.rot(T,1)
        
        one1,one2 = make_ones(cc,i+1)
        H = cc.add(H,one1)
        V = cc.add(V,one2)
        
        #print(H,V)
    
    #acc = cc.add(acc,V)
    
    for i in range(m,2*m) :
        mask = make_mask_backward(cc,i)
        E_m = E[(i-1)%m]
        E_m = cc.mul(E[(i-1)%m],mask)
        
        min = cc.min(E_m,H,V)
        T = 1-V+min
        V = 1-H+min
        H = cc.rot(T,1)
        
        H = cc.mul(H,mask)
        V = cc.mul(V,mask)
    
        acc = cc.add(acc,  cc.rot(V,(m-i)))

    return m + acc[0]

def ED_sim_fold(cc : rlwe, a : List[int],b : List[int]) :
    m = cc.num_slots
    A = cc.SetCiphertext(a)
    B = cc.SetCiphertext(b)

    # B : automorphism invariant
    arrA = re_ordering(cc,A)

    # make Equality Table    
    E = np.zeros((m,m),dtype=np.float64)
    for i in range(m) :
        E[i] = cc.eq(B,arrA[i])
        E[i] = -1*E[i]
    

    # Compute DP table
    H = np.zeros(m,dtype=np.float64)
    V = np.zeros(m,dtype=np.float64)
    H[0] = 1.0
    V[0] = 1.0
    H2 = np.zeros(m,dtype=np.float64)
    V2 = np.zeros(m,dtype=np.float64)
    H2[m-1] = 1.0
    V2[m-1] = 1.0
    
    acc = np.zeros(m,dtype=np.float64)
    acc2 = np.zeros(m,dtype=np.float64)

    for i in range(1,m) :
        mask = make_mask_foward(cc,i)
        E_m = cc.mul(E[i-1],mask)
        
        minEHV = cc.min(E_m,H,V)
        T = 1-V+minEHV
        V = 1-H+minEHV
        T = cc.mul(T,mask)
        V = cc.mul(V,mask)
        
        H = cc.rot(T,1)
        
        one1,one2 = make_ones(cc,i+1)
        H = cc.add(H,one1)
        V = cc.add(V,one2)
        
        
        mask2 = mask[::-1]
        E2_m = cc.mul(E[m-i-1],mask2)
        min2 = cc.min(E2_m,H2,V2)
        T2 = 1-V2+min2
        V2 = 1-H2+min2
        T2 = cc.mul(T2,mask2)
        V2 = cc.mul(V2,mask2)
        H2 = cc.rot(T2,-1)
        
        one1,one2 = make_ones(cc,i+1)
        one1 = one1[::-1]
        one2 = one2[::-1]

        if i!=m-1:
            H2 = cc.add(H2,one1)
            V2 = cc.add(V2,one2)
        acc = cc.add(acc, cc.rot(V,m-i))
        acc2 = cc.add(acc2, cc.rot(V2,-(m-i-1)))
        
    
    temp = acc
    for i in range(m,m+1) :
        mask = make_mask_backward(cc,i)
        E_m = E[(i-1)%m]
        E_m = cc.mul(E[(i-1)%m],mask)
        
        minEHV = cc.min(E_m,H,V)
        T = 1-V+minEHV
        V = 1-H+minEHV
        H = cc.rot(T,1)
        
        H = cc.mul(H,mask)
        V = cc.mul(V,mask)
        temp = cc.add(temp,V)
    
    #print("acc")
    #print(acc)
    #print(acc2)
   
    val = np.arange(1, m+1)[::-1]

    temp = cc.add(temp,val)

    val[0]=0
    val2 = np.arange(0,m)
    #print(val)
    #print(val2)
    acc = cc.add(acc,val)
    acc2 = cc.add(acc2,val2)
    #print(acc,acc2)

    #acc2[0]=0.0
    #print(acc)
    #print(acc2)
   
  
    mid1 = cc.add(acc,acc2)
    mid1[0]=1000.0
    mid2 = cc.add(temp,acc2)
    #print(temp,acc2)
    #print(mid2)
    #mid1[0] = 1000.0
    #mid2 = cc.add(acc,cc.rot(acc2,-1))
    #mid2[m-1] = 1000.0
    
   
    #print(mid1)
    #print(mid2)
    
    #mid3 = cc.min2(mid1,mid2)
    rtn = 10000.0
    import math
    for i in range(1,mid1.shape[0]) :
        rtn = min(mid1[i],rtn)
    
    #for i in range(1,mid1.shape[0]) :
    #    rtn = min(mid2[i],rtn)
    
    return rtn

def ED_sim_pipe(cc : rlwe, a : List[int], b_list : List[List[int]]) :
    debug = False
    m = cc.num_slots
    len = len(b_list)
    A = cc.SetCiphertext(a)
    B = cc.SetCiphertext(b_list[0])


    # B : automorphism invariant
    arrA = re_ordering(cc,A)
    
    if debug : 
        print("re-ordering A")
        print(arrA)

    # make Equality Table    
    E = np.zeros((m,m),dtype=np.float64)
    for i in range(m) :
        E[i] = cc.eq(B,arrA[i])
        E[i] = -1*E[i]
    
    if debug :
        print("Equaltiy Table")
        print(E)

    # Compute DP table
    H = np.zeros(m,dtype=np.float64)
    V = np.zeros(m,dtype=np.float64)
    H[0] = 1.0
    V[0] = 1.0
    acc = np.zeros(m,dtype=np.float64)

    for i in range(1,m) :
        mask = make_mask_foward(cc,i)
        E_m = cc.mul(E[i-1],mask)
        
        min = cc.min(E_m,H,V)
        T = 1-V+min
        V = 1-H+min
        T = cc.mul(T,mask)
        V = cc.mul(V,mask)
        
        H = cc.rot(T,1)
        
        one1,one2 = make_ones(cc,i+1)
        H = cc.add(H,one1)
        V = cc.add(V,one2)
        
        #print(H,V)
    
    for l in range(1,len) :
        for i in range(m,2*m) :
            mask = make_mask_backward(cc,i)
            E_m = E[(i-1)%m]
            #E_m = cc.mul(E[(i-1)%m],mask)
            
            min = cc.min(E_m,H,V)
            T = 1-V+min
            V = 1-H+min
            H = cc.rot(T,1)
            
            H = cc.mul(H,mask)
            V = cc.mul(V,mask)
        
            acc = cc.add(acc,  cc.rot(V,(m-i)))

    
    for i in range(m,2*m) :
        mask = make_mask_backward(cc,i)
        E_m = E[(i-1)%m]
        #E_m = cc.mul(E[(i-1)%m],mask)
        
        min = cc.min(E_m,H,V)
        T = 1-V+min
        V = 1-H+min
        H = cc.rot(T,1)
        
        H = cc.mul(H,mask)
        V = cc.mul(V,mask)
    
        acc = cc.add(acc,  cc.rot(V,(m-i)))

    return m + acc[0]


def re_ordering(cc : rlwe, A : np.ndarray) :
    m = cc.num_slots
    arrA = np.zeros((m,m),dtype=np.float64)
    for i in range(m) :
        if i == 0 :
            arrA[m-1][:i] = A[:i][::-1]
            arrA[m-1][i:] = A[i:][::-1]
        else :
            arrA[i-1][:i] = A[:i][::-1]
            arrA[i-1][i:] = A[i:][::-1]
      
    return arrA

def make_ones(cc : rlwe, idx :int) :
    one1 = cc.SetCiphertext([0]*cc.num_slots)
    one1[0] = 1.0
    one2 = cc.SetCiphertext([0]*cc.num_slots)
    one2[idx-1] = 1.0

    return one1,one2

def make_mask_foward(cc : rlwe, idx : int) : 
    mask = np.zeros(cc.num_slots, dtype=np.float64)
    for i in range(idx) :
        mask[i] = 1.0
    return mask

def make_mask_backward(cc : rlwe, idx : int) : 
    back_idx = idx - cc.num_slots 
    mask = np.ones(cc.num_slots, dtype=np.float64)
    for i in range(back_idx) :
        mask[i] = 0.0
    return mask


