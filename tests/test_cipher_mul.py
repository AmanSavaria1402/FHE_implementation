'''
    Test the ciphertext ciphertext multiplication.
'''
from src.polynomial import *
from src.GSW import *
from src.GSW_utils import *
from src.RLWE import *


def test_mul():
    '''
        Test the ciphertext-ciphertext multiplication
    '''
    # key generation
    config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
    gsw_config = GSWConfig(config, 2**8)
    rlwe_key = RLWEEncryptionKey(config)
    rlwe_key.generate_key()
    gsw_key = GSWEncryptionKey(gsw_config)
    gsw_key.generate_key_from_rlwe_key(rlwe_key)

    # creating gsw polynomial
    fx = build_monomial(gsw_config.rwlecfg.N, 3, 2)
    fx_gsw_cipher = gsw_key.encrypt(fx)
    # creating rlwe polynomial
    f2x = build_monomial(config.N, i=2, c=1)
    plain_poly = encode_rlwe(f2x, config)
    cipher_poly = rlwe_key.encrypt(plain_poly)
    mult_cipher = ciphertext_multiplication(fx_gsw_cipher, cipher_poly)
    decrypt_poly = rlwe_key.decrypt(mult_cipher)
    assert decrypt_poly.coeffs[5] == 2