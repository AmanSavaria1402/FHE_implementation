'''
    Test the RLWE ciphertext code.
'''
from src.RLWE import *

def test_encrypt_decrypt():
    # config
    config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
    # generating the encryption key
    secret_key = RLWEEncryptionKey(config=config)
    secret_key.generate_key()
    # encrypting a polynomial, say 4x^2
    clear_poly = build_monomial(config.N, i=2, c=4)
    # encoding the polynomial
    plain_poly = encode_rlwe(clear_poly, config)
    # encrypt the ciphertext
    cipher_poly = secret_key.encrypt(plain_poly)
    # decrypt the ciphertext
    decrypt_poly = secret_key.decrypt(cipher_poly)
    assert decrypt_poly.coeffs[2] == 4

def test_cipher_add():
    # config
    config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
    # generating the encryption key
    secret_key = RLWEEncryptionKey(config=config)
    secret_key.generate_key()
    # encrypting a polynomial, say 4x^2
    clear_poly = build_monomial(config.N, i=2, c=4)
    clear_poly2 = build_monomial(config.N, i=3, c=3)
    # encoding the polynomial
    plain_poly = encode_rlwe(clear_poly, config)
    plain_poly2 = encode_rlwe(clear_poly2, config)
    # encrypt the ciphertext
    cipher_poly = secret_key.encrypt(plain_poly)
    cipher_poly2 = secret_key.encrypt(plain_poly2)
    # adding
    added_poly = cipher_poly.add_ciphertext(cipher_poly2)
    # decrypt the ciphertext
    decrypt_poly = secret_key.decrypt(added_poly)
    assert decrypt_poly.coeffs[2]==4 and decrypt_poly.coeffs[3]==3

def test_cipher_sub():
    # config
    config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
    # generating the encryption key
    secret_key = RLWEEncryptionKey(config=config)
    secret_key.generate_key()
    # encrypting a polynomial, say 4x^2
    clear_poly = build_monomial(config.N, i=2, c=4)
    clear_poly2 = build_monomial(config.N, i=3, c=3)
    # encoding the polynomial
    plain_poly = encode_rlwe(clear_poly, config)
    plain_poly2 = encode_rlwe(clear_poly2, config)
    # encrypt the ciphertext
    cipher_poly = secret_key.encrypt(plain_poly)
    cipher_poly2 = secret_key.encrypt(plain_poly2)
    # adding
    added_poly = cipher_poly.sub_ciphertext(cipher_poly2)
    # decrypt the ciphertext
    decrypt_poly = secret_key.decrypt(added_poly)
    assert decrypt_poly.coeffs[2]==4 and decrypt_poly.coeffs[3]==-3

def test_plain_mul():
    # config
    config = RLWEConfig(N=1024, sigma=2**(-24), p = 16, q=2**32)
    # generating the encryption key
    secret_key = RLWEEncryptionKey(config=config)
    secret_key.generate_key()
    # encrypting a polynomial, say 4x^2
    clear_poly = build_monomial(config.N, i=2, c=2)
    clear_poly2 = build_monomial(config.N, i=3, c=3)
    # encoding the polynomial
    plain_poly = encode_rlwe(clear_poly, config)
    plain_poly2 = RLWEPlaintext(p=clear_poly2, config=config)
    # plain_poly2 = encode_rlwe(clear_poly2, config)
    # encrypt the ciphertext
    cipher_poly = secret_key.encrypt(plain_poly)
    # cipher_poly2 = secret_key.encrypt(plain_poly2)
    # adding
    added_poly = cipher_poly.mul_plaintext(plain_poly2)
    # decrypt the ciphertext
    decrypt_poly = secret_key.decrypt(added_poly)
    assert decrypt_poly.coeffs[5]==6
