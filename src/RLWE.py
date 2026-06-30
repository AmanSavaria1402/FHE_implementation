'''
    This file contains all the classes and utility functions needed in for RLWE encryption.
    NOTE TO SELF: This file also has the same 'problem' as the LWE file.
'''
from __future__ import annotations
import numpy as np
from src.polynomial import Polynomial, make_zero_polynomial, build_monomial
from src.LWE import LWEEncryptionKey


class RLWEConfig:
    '''
        This class defines the RLWE config.
    '''
    def __init__(self, N, p, q, sigma):
        self.N = N
        self.p = p
        self.q = q
        self.sigma = sigma

class RLWEPlaintext:
    '''
        This class implements the RLWE Plaintext. This just stores the encoded message polynomial. Functionally, it isnt really needed but Ive just added it to make my code consistent with how TFHE libraries work.
    '''
    def __init__(self, p: Polynomial, config: RLWEConfig):
        self.p = p
        self.config = config


# ---- utility functions ----
# utility functions for RLWE
def encode_rlwe(p: Polynomial, config: RLWEConfig):
    '''
        This function encodes the cleartext polynomial and returns it as a Plaintext.
    '''
    # get the delta
    delta = config.q//config.p
    # scale each coefficient of the polynomial by multiplying it with delta
    scaled_coeffs = np.multiply(p.coeffs, delta, dtype=np.int32)
    return RLWEPlaintext(Polynomial(p.N, scaled_coeffs), config=config)

def decode_coeff(c, delta, p):
    '''
        This function decodes a specific coefficient within the polynomial.
    '''
    # divide by delta and round
    decoded_c = int(np.rint(c/delta))
    # calculate decoded_c mod p and you get the final decoded value
    out = ((decoded_c + (p//2)) % p) - (p//2)
    return out

def decode_rlwe(p: Polynomial, config: RLWEConfig):
    '''
        This function decodes the given polynomial to remove the noise and get the final decrypted cleartext.
    '''
    # calcualte the delta
    delta = config.q//config.p
    raw_coeffs = p.coeffs
    # decoded each coeff and get the output polynomial
    decoded_coeffs = np.array([decode_coeff(c, delta,config.p) for c in raw_coeffs])

    return Polynomial(config.N, decoded_coeffs)

def rlwe_generate_uniform_sample(config: RLWEConfig) -> Polynomial:
    '''
        This function generates a polynomial of the degree specified in the config where each coefficient is an integer sampled from a uniform random distribution.
    '''
    # get the limits
    minval = -config.q//2
    maxval = config.q//2

    # generate N values between the min (inclusive) and max (exclusive)
    coeffs = np.random.randint(minval, maxval, size=config.N, dtype=np.int32)

    return Polynomial(config.N, coeffs=coeffs)

def rlwe_generate_normal_sample(config: RLWEConfig) -> Polynomial:
    '''
        This function generates a polynomial of the degree as per the config where each coefficient is drawn from a normal distribution with mean 0 and standard distribution given in the config.
    '''
    coeffs = np.int32(np.multiply(np.random.normal(0, config.sigma, size=config.N), config.q//2))
    return Polynomial(config.N, coeffs)

# ---- utility functions over -----
class RLWECiphertext:
    '''
        This class implements the RLWE ciphertext.
    '''
    def __init__(self, a: Polynomial, b: Polynomial, config: RLWEConfig):
        self.config = config
        self.a = a
        self.b = b

    def add_ciphertext(self, other: RLWECiphertext):
        '''
            Add another ciphertext <other> to the current ciphertext and return the result as a new ciphertext. Two add two ciphertexts, you just add their respective a and b polynomials.
        '''
        # adding a and b, and creating a new ciphertext
        added_a = self.a.polynomial_add(other.a)
        added_b = self.b.polynomial_add(other.b)

        return RLWECiphertext(added_a, added_b, self.config)
    
    def sub_ciphertext(self, other: RLWECiphertext):
        '''
            Subtract another ciphertext <other> from the current one and return the result as a new ciphertext. This operation subtracts <other> from the <self> and not the other way round.
        '''
        # subtracting a and b and returning a new ciphertext
        subbed_a = self.a.polynomial_subtract(other.a)
        subbed_b = self.b.polynomial_subtract(other.b)

        return RLWECiphertext(subbed_a, subbed_b, self.config)
    
    def mul_plaintext(self, cx: RLWEPlaintext):
        '''
            This function multiplys the given plaintext to self and returns the product as a new ciphertext.
        '''
        # for this, you multiply both a and b with the plaintext
        mul_a = self.a.polynomial_multiply(cx.p)
        mul_b = self.b.polynomial_multiply(cx.p)

        return RLWECiphertext(mul_a, mul_b, self.config)


class RLWEEncryptionKey:
    '''
        This class defines the RLWE encryption/secret key and performs encryption and decryption of ciphertexts.
    '''
    def __init__(self, config:RLWEConfig):
        self.config = config
        self.key = Polynomial # this time, the key is also a polynomial.

    def generate_key(self):
        '''
            Generate the key. The key is a random polynomial of degree N with binary coefficients.
        '''
        self.key = Polynomial(
            N=self.config.N,
            coeffs = np.random.randint(
                low = 0,
                high = 2,
                size = self.config.N,
                dtype = np.int32
            )
        )
    
    def encrypt(self, mx: RLWEPlaintext):
        '''
            Encrypt the input Polynomial p into an RLWE Ciphertext.
        '''
        # generate a and noise e
        a = rlwe_generate_uniform_sample(self.config)
        e = rlwe_generate_normal_sample(self.config)

        # encrypting
        a_times_s = a.polynomial_multiply(self.key) # <a, s>
        as_plus_m = a_times_s.polynomial_add(mx.p) # <a, s> + m (m is already encoded, which is why I wrote m and not delta * m)
        b = as_plus_m.polynomial_add(e)

        return RLWECiphertext(a, b, self.config)
    
    def decrypt(self, cx: RLWECiphertext) -> Polynomial:
        '''
            Decrypt the given ciphertext, decode it and return the cleartext polynomial.
        '''
        # separate a and b
        ca = cx.a
        cb = cx.b
        # subtracting <a, s> from b
        subbed = cb.polynomial_subtract(ca.polynomial_multiply(self.key))
        # decode: divide by delta to remove the noise and take mod wrt p
        decoded = decode_rlwe(subbed, self.config)

        return decoded
    
def get_trivial_rlwe_ciphertext(mx: Polynomial, config: RLWEConfig):
    '''
        This function returns the trivial RLWE ciphertext for the given config. Trivial RLWE ciphertext is basically just an RLWE ciphertext where a is just a zero polynomial and b is polynomial itself.
        NOTE: Im not very sure if we should encode the input polynomial, meaning, if we take a cleartext input or plaintext input.
    '''
    # sanity check
    if len(mx.coeffs) != config.N:
        raise ValueError("The polynomial has a degree thats different from the degree provided in the config file. Please ensure that theyre the same.")
    # creating a
    a = make_zero_polynomial(config.N)
    return RLWECiphertext(a = a, b = mx, config=config)

def lwe_to_rlwe_key(lwe_key: LWEEncryptionKey, config: RLWEConfig):
    '''
        This function creates an RLWE key for/from the given RLWE key.
    '''
    # creating a polynomial with the coeffs as the lwe_key
    if lwe_key.config.n != config.N:
        raise ValueError("The size of LWEKey and degree of RLWEConfig do not match.")
    key_poly = Polynomial(N=config.N, coeffs=lwe_key.key)
    rlwe_key = RLWEEncryptionKey(config=config)
    rlwe_key.key = key_poly

    return rlwe_key