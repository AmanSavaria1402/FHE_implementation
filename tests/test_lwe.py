'''
    Testing the lwe function.
'''
import numpy as np
from src.LWE import *

def test_encrypt_decrypt():
    lweconfig = LWEConfig(p=1<<4, q=1<<32, noise_level=2**-20, dimension=1024)
    # encoding
    plain_m = encode_plaintext(-1, lweconfig)
    # secret key
    secret_key = LWEEncryptionKey(lweconfig)
    secret_key.generate_key()
    # encrypting
    m_cipher = secret_key.encrypt(plain_m)
    # decrypting
    m_decrypt = secret_key.decrypt(m_cipher)
    assert m_decrypt.m == -1

def test_plain_add():
    lweconfig = LWEConfig(p=1<<4, q=1<<32, noise_level=2**-20, dimension=1024)
    # encoding
    plain_m = encode_plaintext(-1, lweconfig)
    plain_m2 = encode_plaintext(4, lweconfig)
    # secret key
    secret_key = LWEEncryptionKey(lweconfig)
    secret_key.generate_key()
    # encrypting
    m_cipher = secret_key.encrypt(plain_m)
    # adding plaintext
    added_cipher = m_cipher.add_plaintext(plain_m2)
    # decrypting
    added_decrypt = secret_key.decrypt(added_cipher)
    assert added_decrypt.m == 3

def test_plain_mul():
    lweconfig = LWEConfig(p=1<<4, q=1<<32, noise_level=2**-20, dimension=1024)
    # encoding
    plain_m = encode_plaintext(-1, lweconfig)
    # secret key
    secret_key = LWEEncryptionKey(lweconfig)
    secret_key.generate_key()
    # encrypting
    m_cipher = secret_key.encrypt(plain_m)
    # adding plaintext
    added_cipher = m_cipher.mul_plaintext(4)
    # decrypting
    added_decrypt = secret_key.decrypt(added_cipher)
    assert added_decrypt.m == -4

def test_cipher_add():
    lweconfig = LWEConfig(p=1<<4, q=1<<32, noise_level=2**-20, dimension=1024)
    # encoding
    plain_m = encode_plaintext(-1, lweconfig)
    plain_m2 = encode_plaintext(4, lweconfig)
    # secret key
    secret_key = LWEEncryptionKey(lweconfig)
    secret_key.generate_key()
    # encrypting
    m_cipher = secret_key.encrypt(plain_m)
    m2_cipher = secret_key.encrypt(plain_m2)
    # adding plaintext
    added_cipher = m_cipher.add_ciphertext(m2_cipher)
    # decrypting
    added_decrypt = secret_key.decrypt(added_cipher)
    assert added_decrypt.m == 3

def test_cipher_sub():
    lweconfig = LWEConfig(p=1<<4, q=1<<32, noise_level=2**-20, dimension=1024)
    # encoding
    plain_m = encode_plaintext(-1, lweconfig)
    plain_m2 = encode_plaintext(4, lweconfig)
    # secret key
    secret_key = LWEEncryptionKey(lweconfig)
    secret_key.generate_key()
    # encrypting
    m_cipher = secret_key.encrypt(plain_m)
    m2_cipher = secret_key.encrypt(plain_m2)
    # adding plaintext
    added_cipher = m_cipher.sub_ciphertext(m2_cipher)
    # decrypting
    added_decrypt = secret_key.decrypt(added_cipher)
    assert added_decrypt.m == -5