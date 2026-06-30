import numpy as np
from src.RLWE import *
from src.polynomial import *
from src.LWE import *
from src.bootstrap import *

# config
config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
lweconfig = LWEConfig(p=1<<4, q=1<<32, noise_level=2**-20, dimension=1024)
lwe_key = LWEEncryptionKey(config=lweconfig)
lwe_key.generate_key()
# get rlwe key from this
rlwe_key = lwe_to_rlwe_key(lwe_key, config)
fx = build_monomial(config.N, 1, -2) # 2x
fx_plain = encode_rlwe(fx, config)
fx_rlwe = rlwe_key.encrypt(fx_plain)

extracted_lwe = extract_sample(i = 1, FX = fx_rlwe)
plain_lwe = lwe_key.decrypt(extracted_lwe)
print(plain_lwe)