'''
    This file contains code for running GSW encryption.
'''
import numpy as np
from src.polynomial import *
from src.RLWE import *
from src.GSW_utils import *

class GSWConfig:
    '''
        This class implements the GSW Config. It takes the RLWE config, as most of the configurations are identical. Additionally, it also takes gsw_p, which will be used for the base p representation.
    '''
    def __init__(self, rlwecfg: RLWEConfig, gsw_p: int):
        self.rwlecfg = rlwecfg
        self.gsw_p = gsw_p

class GSWPlaintext:
    '''
        This class represents the gsw plaintext.
    '''
    def __init__(self, config: GSWConfig, message: Polynomial):
        self.config = config
        self.p = message


class GSWCiphertext:
    '''
        This class implements the gsw ciphertext.
    '''
    def __init__(self, config: GSWConfig, polys: Sequence[RLWECiphertext]):
        self.config = config
        self.polys = polys


class GSWEncryptionKey:
    '''
        This class implements the GSW encryption key and methods to encrypt and (possibly) decrypt a gsw ciphertext.
    '''
    def __init__(self, config: GSWConfig):
        self.config = config
        self.key = Polynomial # this is just going to be the RLWE encryption key s(x)
        self.rlwe_key = RLWEEncryptionKey

    def generate_key_from_rlwe_key(self, rlwe_key: RLWEEncryptionKey):
        '''
            This function creates the GSWEncryptionKey from the RLWEEncryption key. To be honest, its just the RLWEEncryption key as it is, but for the sake of structure of code and to keep it consistent with other classes. Its a little silly to save both the rlwe_key object and the key polynomail as well, but eh...
        '''
        self.key = rlwe_key.key
        self.rlwe_key = rlwe_key

    def get_rlwe_key_from_gsw(self):
        '''
            This function returns the rlwe key from the current gsw_key. Again, its the same, but consistency, formality, something, something...
        '''
        rlwe_key = RLWEEncryptionKey(config=self.config.rwlecfg)
        rlwe_key.key = self.key
        return rlwe_key
    
    
    def encrypt(self, fx: Polynomial):
        '''
            Encrypt a given polynomial as a GSWCiphertext.
        '''
        # first, we need to get the value of k
        k = get_k(self.config.rwlecfg.q, self.config.gsw_p)
        gsw_p = np.int32(self.config.gsw_p) # i dont want to write the whole thing
        # make Z: the zero encryption array of shape (2k, 2), note that each row is a RLWE ciphertext, the 2 just means (a(x), b(x))
        zero_poly = make_zero_polynomial(self.config.rwlecfg.N)
        zero_poly_plain = encode_rlwe(zero_poly, self.config.rwlecfg) # while for 0 config, the coefficients for both the plaintext and ciphertext are going to be the same, the encryption for RLWE only accepts an RLWEPlaintext and not a Polynomial
        Z = [
            self.rlwe_key.encrypt(zero_poly_plain) for _ in range(2*k)
        ] # this has shpae (2k, 2)
        # create the GSW ciphertext
        # its created using the equation: fx * B_p + Z
        # i will implicitly divide the Z matrix in two halves 
        # here, ill iterate for k times, and for each iteration, the corresponding index's RLWE of the first half's a gets fx * (p**i) added and the corresponding index of the second half's b gets the same added
        for i in range(k):
            # multiply fx with the correct power of p
            scaled_fx = fx.polynomial_const_multiply(gsw_p**i)
            # add it to first half's a
            Z[i].a = Z[i].a.polynomial_add(scaled_fx)
            # add it to second half's b
            Z[i + k].b = Z[i + k].b.polynomial_add(scaled_fx)

        return GSWCiphertext(self.config, Z)
    
def ciphertext_multiplication(G: GSWCiphertext, R: RLWECiphertext) -> RLWECiphertext:
    '''
        This function performs the ciphertext-ciphertext multiplication between G and R. G is of shape (2k, 2), basically 2k RLWE Ciphertexts and each Rlwe ciphertext is (a, b).
        So a little explanation because this function looks a little confusing (but in reality it really isnt once you look at it hard enough keeping in mind how the multiplication is actually done). Im going to define in the docstring how does this operation actually works. 
        NOTE TO SELF: If its been long since this code was studied last, I practically beg you to please read the thing after so in its entirety. Otherwise, there is a real chance that something in your brain explodes (or it starts to hurt) and I dont want that because I care for you :)
        Also, I know that im calling using the term vector of polynomials, which might sound wierd becasue a polynomial itself is a vector of coefficients, so a vector of polynomials is just a higher dimension matrix, but in FHE, it seems that polynomials are referred to as a single entity, adressed like constant numbers, so Im just keeping it consistent.
        So:
        First, G is a series of 2k RLWECiphertexts. Essentially, theyre all zero ciphertexts to which fx (the polynomial with which we actually want to achieve multiplication with R) is added. For the first half of G [0:k-1], fx is added to a(x) of each RLWE 'row', and is added to b(x) for the second half [k: 2k-1].
        Second, a(x) and b(x) of R are converted to corresponding base-p representation. Now R was basically a vector/sequence/whatever of 2 polynomials, but now it becomes a vector of 2k polynomials (k polynomials for base p representation of a(x) and b(x)). 
        To perform c-c multiplication, you basically do something akin to matrix multiplication, (Basep(a(x)), Basep(b(x))) @ G, which gives you a single rlwe polynomial, and its shape will be (1, 2), a single (a(x), b(x)) ciphertext. Here, each multiplication is polynomial multiplication followed by polynomial addition for matmul.
    '''
    gsw_config = G.config
    rlwe_config = R.config
    # first get the base p representation for R, for both a(x) and b(x) and create it into a sequence of 2k polynomials
    R_base_p_repr = get_base_p_from_polynomial(R.a, p=gsw_config.gsw_p, q=rlwe_config.q) + get_base_p_from_polynomial(R.b, p=gsw_config.gsw_p, q=rlwe_config.q) # the output of the function is a list so + just appends
    
    # so matmul here is implemented using a for loop, first, an RLWECiphertext that is basically all zeros is initialized and then the corresponding products a.G and b.G for each R_base_p_repr is added to a and be of this new RLWECiphertext respectively
    res_rlwe_cipher = RLWECiphertext(
        config = rlwe_config,
        a = make_zero_polynomial(rlwe_config.N),
        b = make_zero_polynomial(rlwe_config.N)
    )

    for i, rp in enumerate(R_base_p_repr):
        # multiply the current r to both a and b of the current row of g.a
        a_mul = rp.polynomial_multiply(G.polys[i].a)
        # adding
        res_rlwe_cipher.a = res_rlwe_cipher.a.polynomial_add(a_mul)

        # doing the same for b
        b_mul = rp.polynomial_multiply(G.polys[i].b)
        # adding
        res_rlwe_cipher.b = res_rlwe_cipher.b.polynomial_add(b_mul)

    return res_rlwe_cipher