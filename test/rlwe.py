import numpy as np
import copy as cp

class rlwe:
    def __init__(self, num_slots) :
        self.num_slots = num_slots 

    def exception_haddling(self, a : np.ndarray) :
        if a.shape[0] != self.num_slots : 
            raise ValueError(f"num_slot mismatch : {self.num_slots} vs {a.shape[0]}")
        return 

    def SetCiphertext(self, v : list) :
        if len(v) != self.num_slots : 
             raise ValueError(f"num_slot mismatch : {self.num_slots} vs {len(v)}")
        
        return np.array(v, dtype=np.float64)

    def add(self, a: np.ndarray, b: np.ndarray) :
        self.exception_haddling(a)
        self.exception_haddling(b)
        return a+b
    
    def mul(self,a: np.ndarray, b: np.ndarray) :
        self.exception_haddling(a)
        self.exception_haddling(b)
        return a*b
    
    def rot(self,a:np.ndarray, r_idx : int) :
        self.exception_haddling(a)
        a_cp = cp.copy(a)
        return np.roll(a_cp,r_idx)
    
    def conj(self,a:np.ndarray) : 
        self.exception_haddling(a)
        return a.conj()
    
    def eq(self,a: np.ndarray, b: np.ndarray) :
        self.exception_haddling(a)
        self.exception_haddling(b)
        return (a == b).astype(np.int32)
    
    def min(self, a:np.ndarray, b:np.ndarray, c:np.ndarray) :
        self.exception_haddling(a)
        self.exception_haddling(b)
        self.exception_haddling(c)
        
        r = np.minimum(a, b)
        r = np.minimum(r, c)

        return r
    
    def min2(self, a:np.ndarray, b:np.ndarray) :
        self.exception_haddling(a)
        self.exception_haddling(b)
    
        
        r = np.minimum(a, b)

        return r


