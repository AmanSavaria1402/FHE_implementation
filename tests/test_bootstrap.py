'''
    This function tests bootstrapping code
'''

import numpy as np
from src.LWE import *
from src.RLWE import *
from src.GSW_utils import *
from src.GSW import *
from src.bootstrap import *

def test_blindrotate():
    # configs
    lweconfig = LWEConfig(p=1<<4, q=1<<32, noise_level=2**-20, dimension=1024)
    config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
    gsw_config = GSWConfig(config, 2**8)
    lwe_key = LWEEncryptionKey(config=lweconfig)
    lwe_key.generate_key()
    rlwe_key = lwe_to_rlwe_key(lwe_key, config)
    # get the gsw key
    gsw_key = GSWEncryptionKey(gsw_config)
    gsw_key.generate_key_from_rlwe_key(rlwe_key)
    # get the bootstrap key
    boostrap_key = BootstrapKey(gsw_config)
    boostrap_key.generate_bootstrap_key(lwe_key=lwe_key, gsw_key=gsw_key)
    
    # generating polynomial to rotate
    N = config.N
    fx = Polynomial(N=N, coeffs=np.ones(N, dtype=np.int32))
    fx.coeffs[: N // 2] = -1
    # encrypt fx
    fx_plain = encode_rlwe(fx, config)
    fx_rlwe = rlwe_key.encrypt(fx_plain)
    i_plain = encode_plaintext(6, lweconfig)
    i_lwe = lwe_key.encrypt(i_plain)    
    # rotating
    rotated_rlwe = boostrap_key.blindrotate(fx_rlwe, i_lwe)
    rotated_clear = rlwe_key.decrypt(rotated_rlwe)
    assert rotated_clear.coeffs[0] == 1
    assert rotated_clear.coeffs[N // 2] == -1
    assert rotated_clear.coeffs[-1] == -1

def test_extract_sample():
    config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
    lweconfig = LWEConfig(p=1<<4, q=1<<32, noise_level=2**-20, dimension=1024)
    lwe_key = LWEEncryptionKey(config=lweconfig)
    lwe_key.generate_key()
    # get rlwe key from this
    rlwe_key = lwe_to_rlwe_key(lwe_key, config)
    fx = build_monomial(config.N, 1, -2) # 2x
    fx_plain = encode_rlwe(fx, config)
    fx_rlwe = rlwe_key.encrypt(fx_plain)
    # extracting the index at i=1
    extracted_lwe = extract_sample(i = 1, FX = fx_rlwe)
    # decrypting
    plain_lwe = lwe_key.decrypt(extracted_lwe)
    assert plain_lwe.m == -2

def test_nand():
    # declaring the configs
    lweconfig = LWEConfig(p=1<<3, q=1<<32, noise_level=2**-20, dimension=1024)
    rlweconfig = RLWEConfig(N=1024, sigma=2**(-24), p = 1<<3, q=2**32)
    gswconfig = GSWConfig(rlweconfig, 2**8)
    # NOTE: While its a good idea to have sigma/noise_level the same for lwe and rlwe, but here, im going to use the lwe key to generate rlwe key, so the noise_level of lwe will be used everywhere so its okay
    # generating the keys
    lwe_key = LWEEncryptionKey(lweconfig)
    lwe_key.generate_key()
    # get the rlwe key from this
    rlwe_key = lwe_to_rlwe_key(lwe_key, rlweconfig)
    # get the gsw key
    gsw_key = GSWEncryptionKey(gswconfig)
    gsw_key.generate_key_from_rlwe_key(rlwe_key)
    # get the bootstrap key
    bsk =BootstrapKey(gswconfig)
    bsk.generate_bootstrap_key(lwe_key=lwe_key, gsw_key=gsw_key)
    # bootstrapping
    # first, we need to create the test polynomial and encrypt it
    tx = generate_tx(rlweconfig.N)
    # create Encode(2)
    scale_plain = encode_plaintext(2, lweconfig)
    # we are actually working with scale_plain//2
    scale_plain.m = scale_plain.m//2
    # use it to encrypt tx as a trivial rlwe ciphertext
    tx_rlwe = create_trivial_rlwe_ciphertext(tx.polynomial_const_multiply(scale_plain.m), rlweconfig)
    # create trivial lwe ciphertext for scale_plain
    scale_lwe = create_trivial_lwe_ciphertext(scale_plain, lweconfig)

    # create encrypted inputs: 0==false, 2==true
    l0_plain = encode_plaintext(2, lweconfig) # true or false
    l1_plain = encode_plaintext(2, lweconfig) # false or false
    # ciphertexts
    l0_cipher = lwe_key.encrypt(l0_plain)
    l1_cipher = lwe_key.encrypt(l1_plain)
    # constant cipher
    const_plain = encode_plaintext(-3, lweconfig)
    const_cipher = lwe_key.encrypt(const_plain)

    nand_cipher = homomorphic_nand(l0_cipher, l1_cipher, const_cipher, scale_lwe, tx_rlwe, bsk)
    # decrypting the nand_cipher
    nand_clear = lwe_key.decrypt(nand_cipher)
    assert nand_clear.m == 0
