'''
    Code to test cmux
'''
from src.polynomial import *
from src.RLWE import *
from src.GSW import *
from src.cmux import *
from src.bootstrap import *

def test_cmux():
    # config and keygen
    config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
    gsw_config = GSWConfig(config, 2**8)
    rlwe_key = RLWEEncryptionKey(config)
    rlwe_key.generate_key()
    gsw_key = GSWEncryptionKey(gsw_config)
    gsw_key.generate_key_from_rlwe_key(rlwe_key)
    # The selector bit is b=1
    b = build_monomial(c=1, i=0, N=config.N)

    # The lines are: l_0(x) = x, l_1(x) = 2x
    l0 = build_monomial(c=1, i=1, N=config.N)
    l1 = build_monomial(c=2, i=3, N=config.N)
    l0_plaintext = encode_rlwe(l0, config)
    l1_plaintext = encode_rlwe(l1, config)

    # encrypt them
    b_gsw_cipher = gsw_key.encrypt(b)
    l0_cipher = rlwe_key.encrypt(l0_plaintext)
    l1_cipher = rlwe_key.encrypt(l1_plaintext)

    cmux_res = cmux(b_gsw_cipher, l0_cipher, l1_cipher)
    cmux_decrypt = rlwe_key.decrypt(cmux_res)
    assert cmux_decrypt.coeffs[3] == 2