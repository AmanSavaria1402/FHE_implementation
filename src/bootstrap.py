'''
    This file contains code for bootstrapping
'''
import numpy as np
from src.RLWE import *
from src.GSW import *
from src.LWE import *
from src.cmux import *



class BootstrapKey:
    '''
        This class implements the bootstrap keys and the methods associated with it, like BlindRotate.
    '''
    def __init__(self, config: GSWConfig):
        self.config = config
        self.key = Sequence[GSWCiphertext] # each component (s_i) from the LWE key will be a GSW ciphertext

    def generate_bootstrap_key(self, lwe_key: LWEEncryptionKey, gsw_key: GSWEncryptionKey):
        '''
            generate the bootstrapping keys. The bootstrapping keys are just GSW encryption of each component of the LWE_key using the gsw_key provided.
            NOTE TO SELF: Im debating whether to make the lwe_key and gsw_key members of this class or not. Right now I wont, but if there are other operations that require these keys in addition to the bootstrap keys, I will.
        '''
        bs_keys = [] # this will store the keys
        # looping across each lwe key component and encrypting it
        for si in lwe_key.key:
            # first, create a monomial with the si being constant
            si_monomial = build_monomial(self.config.rwlecfg.N, 0, si)
            si_gsw = gsw_key.encrypt(si_monomial)
            bs_keys.append(si_gsw)
        self.key = bs_keys

    def blindrotate(self, FX: RLWECiphertext, I: LWECiphertext) -> RLWECiphertext:
        '''
            This function implements the Rotate(I, FX) function homomorphically. This function, given i, and a polynomial f(x) (I = LWE(i), FX = RLWE(f(x)), calculates x^i . f(x) as an RLWECiphertext.
        '''
        # scale the lwe ciphertext 
        scaled_I_a = np.int32(np.rint(np.multiply(I.a, ((2 * self.config.rwlecfg.N)/ self.config.rwlecfg.q)))) 
        scaled_I_b = np.int32(np.rint(np.multiply(I.b, ((2 * self.config.rwlecfg.N)/ self.config.rwlecfg.q)))) # NOTE: This is going to be a constant since this is an LWE ciphertext

        # rotate operation
        # initialize g0
        xb_poly = build_monomial(self.config.rwlecfg.N, scaled_I_b, 1)
        xb_plain = RLWEPlaintext(xb_poly, FX.config)        
        rotated_poly = FX.mul_plaintext(xb_plain) # this will updated every timestep
        # loop N times to finally get the rotated output
        for i, ai in enumerate(scaled_I_a):
            # get g_i-1 (x) * x^a_i
            poly = RLWEPlaintext(build_monomial(self.config.rwlecfg.N, -ai, 1), self.config.rwlecfg) # just converting it to a plaintext because the plaintext multiply in RLWE wants a plaintext
            right = rotated_poly.mul_plaintext(poly)
            # update rotated_poly using the cmux operation
            rotated_poly = cmux(self.key[i], rotated_poly, right)
        
        return rotated_poly
    

def extract_sample(i: int, FX: RLWECiphertext) -> LWECiphertext:
    '''
        This function extracts the ith coefficient of the polynomial encrypted in the input FX and returns the ith coefficient as an LWE Polynomial. It calculates the ai (a for the LWE ciphertext) using Conv(a(x), i) and the ith coefficient of b(x). The resulting LWE ciphertext is encrypted using the LWE key obtained from the RLWE key (you can convert them into one another).
    '''
    # first get the conv(a(x), i)
    # conv(a(x), i) is just going to be coeffs of a(x), but 'rotated' around x.
    # that means, if a(x) was [a0, a1, a2, ...., a(N-1)]
    # then conv(a(x), i) is [ai, a_i-1, ..., a1, a0, -aN-1, -aN-1, ..., -ai+1]
    ax = FX.a.coeffs
    conv_ax = np.hstack([
        ax[:i+1][::-1],
        -1 * ax[i+1:][::-1]])
    bi = FX.b.coeffs[i] # just need the ith coefficient's b

    # creating an LWECiphertext for this
    return LWECiphertext(
        a = conv_ax,
        b = bi
    )


def bootstrap(i: LWECiphertext, bsk: BootstrapKey, R:RLWECiphertext, scale_lwe: LWECiphertext):
    '''
        This function implements the bootstrap function to calculate the step operation defined above for the input i. The bsk is the bootstrapping key.
        This function returns 0 if i is in (-q/4, q/4], else, it returns an lwe encryption of scale.
        NOTE [IMPORTANT]: Note that if scale is the cleartext np.int32 value, scale_lwe encrypts Encode(scale)/scale.
    '''
    # blind rotate R
    rotated_rlwe = bsk.blindrotate(R, i)
    # extract sample to get the coefficient of constant
    coeff_lwe = extract_sample(0, rotated_rlwe)
    # adding scale_lwe to the output to get the offset-ed result of bootstrapping
    result_lwe = scale_lwe.add_ciphertext(coeff_lwe)

    return result_lwe

def generate_tx(N) -> Polynomial:
    '''
        This function generates the polynomial i(t).
    '''
    # building the polynomial
    tx_coeffs = np.ones(N, dtype=np.int32)
    # updating the first half to be -1
    tx_coeffs[:N//2] = -1
    tx_poly = Polynomial(N, tx_coeffs)
    return tx_poly

def homomorphic_nand(b0: LWECiphertext, b1: LWECiphertext, constCipher: LWECiphertext, scale_lwe: LWECiphertext, txR:RLWECiphertext, bsk:BootstrapKey) -> LWECiphertext:
    '''
        This function implements the nand operation between two inputs b0 and b1 homomorphically.
        The way it does it is by first evaluating the value of F(b0, b1) = LWE(-3) - b0 - b1 and then running the bootstrap that implements the step function which evaluates the final nand function.
        Remember that 0 is equivalent to F while 2 is equivalent to T here. The output will be an LWE ciphertext encrypting 0 or 2 as well.
        constCipher is just the LWE Encryption of -3.
    '''
    # calculating constCipher - b0 - b1
    cc_minus_b0 = constCipher.sub_ciphertext(b0)
    F_lwe = cc_minus_b0.sub_ciphertext(b1)
    # call bootstrapping on this
    nand_lwe = bootstrap(F_lwe, bsk, txR, scale_lwe)

    return nand_lwe

def create_trivial_rlwe_ciphertext(mx: Polynomial, config: RLWEConfig):
    '''
        This function creates a trivial RLWE ciphertext. A trivial RLWE ciphertext is one where the polynomial is one where b is the polynomial itself while a is a zero polynomial.
    '''
    a = make_zero_polynomial(config.N)
    return RLWECiphertext(a, mx, config)

def create_trivial_lwe_ciphertext(m: LWEPlaintext, config: LWEConfig):
    '''
        This function creates a trivial LWE ciphertext, which is an LWECiphertext with a being a vector of zeros and b being the plaintext itself.
    '''
    a = np.zeros(config.n, dtype=np.int32)
    return LWECiphertext(a, m.m)