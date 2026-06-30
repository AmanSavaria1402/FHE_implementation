'''
    This file contains code for cmux. And any (future) cmux related things.
'''

import numpy as np
from src.polynomial import Polynomial
from src.RLWE import *
from src.GSW import *
from src.GSW_utils import *

def cmux(b: GSWCiphertext, l0: RLWECiphertext, l1:RLWECiphertext) -> RLWECiphertext:
    '''
        This function implements the multiplexer for ciphertexts.
        The homomorphic multiplexer operation homomorphically computes: b(l1 - l0) + l0.
        So its basically chaining the following operations: CAdd(CMul(b, CSub(l1, l0)), l0)
        the function returns l0 is b == 0, else l1.
    '''
    # calculating l1 - l0
    subbed = l1.sub_ciphertext(l0)
    # calculating product of subbed and b
    prod = ciphertext_multiplication(b, subbed)
    # adding l0 to get the cmux result
    mux_cipher = prod.add_ciphertext(l0)

    return mux_cipher