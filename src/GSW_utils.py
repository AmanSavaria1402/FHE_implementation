'''
    This file contains utility code for the GSW encryption.
'''
import numpy as np
from typing import Sequence
from src.polynomial import *

def get_k(q, p):
    '''
        Given q (ciphertext modulus) and p (not plaintext modulus, but a different parameter specifically for base p representation used for ciphertext-ciphertext multiplication of an RLWE and GWE cipherterxts), get the value of k, which is the sequence length of base-p representation of any input.
    '''
    return np.int32(np.log2(q)/np.log2(p))

def convert_array_to_base_p(a: np.ndarray, p, q):
    '''
        This function converts each element into a its base p-representation, given the values of q (ciphertext modulus) and p (not plaintext modulus). Note that the representation for each element in a will be a sequence of a fixed length.
    '''
    log_p = np.log2(p).astype(np.int32)
    k = get_k(q = q, p = p)
    p_half = p//2
    offset = p_half * np.sum([p**i for i in range(k)])
    mask = p - 1 # this in binary is exactly log_p ones.

    a_offset = np.uint32(a + offset) # need to set this to uint 32 as the the p-bit representation will be calculated from unsigned integers and then converted to signed
    # return a_offset
 
    output = [] # this will store the output base-p representation for each element
    # NOTE: output is a list of k numpy arrays. Each numpy array stores the ith p-bit representation of the corresponding original input.
    for i in range(k):
        output.append(
            (np.right_shift(a_offset, i * log_p) & mask).astype(np.int32) - np.int32(p_half)
            )
        
    return output

def convert_base_p_to_int_array(a: Sequence[np.array], p: np.int32):
    '''
        This function takes as input the p-bit representations of an array, and converts it back into integer. a should be a list or a iterator of ndarrays where each ndarray contains the ith p-bit representation of the integer.
    '''
    # calculate x from the representation: x = x_0 + x_1 p + ... + x_k-1 p^k-1
    decimal_repr = np.zeros(len(a[0]), dtype=np.int32)
    for i in range(len(a)):
        decimal_repr += (a[i] * (p**i)).astype(np.int32)
    return decimal_repr

# defining base p representation functions for polynomials
def get_base_p_from_polynomial(poly: Polynomial, p: np.int32, q: np.int32) -> Sequence[Polynomial]:
    '''
        This function converts each coefficient of the coefficients of the given polynomial and returns them as a set of polynomials (sequence of the polynomials is very important).
    '''
    coeffs = poly.coeffs
    # converting coeffs into base p representation
    coeffs_base_p = convert_array_to_base_p(coeffs, p, q)
    # create a polynomial for each array
    base_p_poly = [Polynomial(N=poly.N, coeffs = ic) for ic in coeffs_base_p]

    return base_p_poly

def get_poly_from_base_p(base_p_polys: Sequence[Polynomial], p) -> Polynomial:
    '''
        This function calculates the polynomial from the given base p representation polynomials.
    '''
    # get all the coefficients from the sequence
    base_p_coeffs = [p.coeffs for p in base_p_polys]
    # get the integer coeffs from the list
    int_coeffs = convert_base_p_to_int_array(base_p_coeffs, p)
    return Polynomial(N=base_p_polys[0].N, coeffs = int_coeffs)