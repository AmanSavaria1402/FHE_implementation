'''
    This file contains the classes for LWE ciphertext.
    NOTE TO SELF: I REALLY want to separate the utiliy functions and classes in separate files, because it makes sense and looks better. But theyre coupled too tightly, and Ill have to do extra work to get it like that.
    TODO: Do what I said the the above NOTE.
'''
from __future__ import annotations
import numpy as np

# --- THE UTILITY FUNCTIONS ---
def generate_uniform_sample(config: LWEConfig):
    '''
        Generate an array of {size} elements using a uniform random distribution in the required range
    '''
    # get the limits
    low = (config.q//2) * -1
    high = (config.q//2) # not subtracting -1 because high is exclusive in random.randint
    # get the random array
    a = np.random.randint(low=low, high=high, size=config.n, dtype=np.int32)
    return a

def generate_normal_sample(config: LWEConfig):
    '''
        This function generates the normal sample required for encryption within the given range. Since this is LWE, we only need 1 element here as m = 1.
    '''
    ei = np.int32(np.random.normal(loc=0.0, scale=config.sigma) * (config.q//2))
    return ei

def encode_plaintext(m, config):
    '''
        This function encodes the message before its enrypted.
    '''
    # get the delta
    delta = config.q//config.p
    encoded = np.int32(delta * m)
    return LWEPlaintext(encoded)

def decode_plaintext(m, config):
    '''
        This function decodes the decrypted plaintext to remove the noise and get the final output
    '''
    # get the delta
    delta = config.q/config.p
    # round to remove the noise
    decoded = int(np.rint(m/delta))
    # need to take a mod wrt p centered around 0 to get the final output
    final = ((decoded + (config.p/2)) % config.p) - (config.p/2)
    return final


# creating the config class
class LWEConfig:
    '''
        This class defines the configurations for encryption
    '''
    def __init__(self, p, q, noise_level, dimension):
        '''
            p: plaintext modulus
            q: ciphertext modulus
            noise_level: The noise level (used as standard deviation for generating the noise)
            dimension: m (number of elements in the sk)
        '''
        self.p = p
        self.q = q
        self.sigma = noise_level
        self.n = dimension

# creating the class to generate plaintext
class LWEPlaintext:
    '''
        This class contains the plaintext to be encrypted.
    '''
    def __init__(self, m):
        '''
            This class (for now), just stores the message wrapping it in a class. It might contain code for encoding it too later but i dont know
        '''
        self.m = m
        
    def __str__(self):
        return str(self.m)
    
    def __repr__(self):
        return str(self.m)
    

class LWECiphertext:
    '''
        This class defines the ciphertext.
    '''
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def add_ciphertext(self, other: LWECiphertext):
        '''
            This function adds a given ciphertext to the current ciphertext and returns a new ciphertext.
        '''
        # adding their as and bs
        added_a = np.add(self.a, other.a, dtype=np.int32)
        added_b = np.add(self.b, other.b, dtype=np.int32)

        # # taking mods just to see if this works, but we have to take mod wrt q
        # q = config.q
        # added_a = ((added_a + (q//2)) % q) - (q//2)
        # added_b = ((added_b + (q//2)) % q) - (q//2) # im getting overflow errors when doing this

        return LWECiphertext(added_a, added_b)
    
    def sub_ciphertext(self, other:LWECiphertext):
        '''
            This function subtracts the other ciphertext from the current one.
        '''
        sub_a = np.subtract(self.a, other.a, dtype=np.int32)
        sub_b = np.subtract(self.b, other.b, dtype=np.int32)
        return LWECiphertext(sub_a, sub_b)
    
    def mul_plaintext(self, other:np.int32):
        '''
            Multiply the other integer to current ciphertext
        '''
        # we just need to multiply the give plaintext to both a and b
        mul_a = np.multiply(other, self.a, dtype=np.int32)
        mul_b = np.multiply(other, self.b, dtype=np.int32)
        return LWECiphertext(mul_a, mul_b)
    
    def add_plaintext(self, other:LWEPlaintext):
        '''
            Add the other plaintext to the current ciphertext. This function assumes that the input plaintext message is already encoded to R_q.
        '''
        # we just need to scale other and then add it to b
        addb = np.add(self.b, other.m, dtype=np.int32)
        return LWECiphertext(self.a, addb)
    
class LWEEncryptionKey:
    '''
        This class implements the LWEEncryption key. The encryption key implements creating the key, which is just a random binary vector or n elements and encrypting a ciphertext form a plaintext.
    '''

    def __init__(self, config:LWEConfig):
        self.config = config
        self.key = np.ndarray # this just initializes key to be an array, which will be then populated
    
    def generate_key(self):
        '''
            This function generates the secret key.
        '''
        self.key = np.random.randint(low=0, high=2, size=self.config.n)

    def encrypt(self, plaintext:LWEPlaintext):
        '''
            Encrypt the plaintext message
        '''
        # generate a and noise
        a = generate_uniform_sample(self.config)
        e = generate_normal_sample(self.config)

        # get b
        b = np.dot(a, self.key) + plaintext.m + e 

        return LWECiphertext(a, b)
    
    def decrypt(self, cipher:LWECiphertext):
        '''
            This function decrypts the ciphertext generated using the same secret key.
        '''
        b = cipher.b
        a = cipher.a

        inter = b - np.dot(a, self.key) # gives m + e, and not just m
        # return LWEPlaintext(inter)

        # round and get the final decrypted output
        final = decode_plaintext(inter, self.config)

        return LWEPlaintext(int(final))