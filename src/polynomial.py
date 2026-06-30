'''
    This file contains code for the polynomial class.
'''
from __future__ import annotations
import numpy as np

class Polynomial:
    '''
        This class represents a polynomial in the ring Z_q[x]/(x^n + 1)
        Note that the coefficients here are in increasing order of degrees of x, [constant, x, x**2, ..., x**N]
    '''

    def __init__(self, N, coeffs: np.ndarray):
        self.N = N
        self.coeffs = coeffs

    def __str__(self):
        '''
            Print the ciphertext as an equation.
        '''
        prtstr = ""
        nonzeroct = 0
        for i, c in enumerate(self.coeffs):
            if c:
                if nonzeroct:
                    prtstr += " + "
                prtstr += f"{c}x^{i}"
                nonzeroct += 1

        return prtstr

    def polynomial_const_multiply(self, c:np.int32):
        '''
            Multiply the current polynomial with the given constant integer and return it as a new polynomial.
        '''
        new_p = Polynomial(self.N, np.multiply(c, self.coeffs, dtype=np.int32))
        return new_p
    
    def polynomial_multiply(self, p2: Polynomial):
        '''
            Multiply the current polynomial with the given poylnomial and return the new polynomial.
        '''
        # get the N from p2
        p2N = p2.N
        # multiply and pad the results to have 2N-1 (thats the longest polynomial you can get)
        # the np.polymul function requires the coefficients in decreasing order of powers or x, se we reverse the current coeff array and then reverse the final output's back again.
        poly_prod = np.polymul(self.coeffs[::-1], p2.coeffs[::-1])[::-1]
        # creating the final answer with padded coeffs
        poly_padded = np.zeros(2*self.N - 1, dtype=np.int32)
        poly_padded[:poly_prod.shape[0]] = poly_prod

        # now, taking modulus of the result (poly_padded) wrt x^N + 1
        result = poly_padded[:self.N]
        result[:-1] -= poly_padded[self.N:]
        return Polynomial(self.N, result)
    
    def polynomial_add(self, p2: Polynomial):
        '''
            This function adds a polyonmial to the given polynomial and returns the output as a new polynomial.
        '''
        add_coeffs = np.add(self.coeffs, p2.coeffs, dtype=np.int32)
        return Polynomial(self.N, add_coeffs)
    
    def polynomial_subtract(self, p2: Polynomial):
        '''
            This function subtracts the given polynomial from the current one and returns the new polynomial.
            Note: Its self - p2 and not p2 - self.
        '''
        sub_coeffs = np.subtract(self.coeffs, p2.coeffs, dtype=np.int32)
        return Polynomial(self.N, sub_coeffs)
    
def make_zero_polynomial(N:int):
    '''
        Declare a zero polynomial and return it. Basically a polynomial where all the coefficients are zeros.
    '''
    return Polynomial(N = N, coeffs = np.zeros(N, dtype=np.int32))

def build_monomial(N, i, c):
    '''
        This function builds a monomial c * x^i in the ring R[x]/(x^N + 1)
        NOTE: Rememeber that the coeffs are in increasing order of powers of x. [constant, x, x**2, ..., x**N]
    '''
    coeffs = np.zeros(N, dtype=np.int32)

    # find k such that 0 <= i + k*N < N
    i_mod_n = i % N
    k = (i_mod_n - i)//N

    # if k is odd then the monomial gets a negative sign because
    # x^i = (-1)^k + x^(i + k*N) = (-1)^k + x^(i % N)
    sign = 1 if k % 2 == 0 else -1
    coeffs[i_mod_n] = sign * c

    return Polynomial(N, coeffs)