"""
Benchmark key generation, encryption and decryption.

"""

import random
import resource
import time
import phe.paillier as paillier
import os.path as osp
import os
try:
    import gmpy2
    HAVE_GMP = True
except ImportError:
    HAVE_GMP = False


def getprimeover(N):
    """Return a random N-bit prime number using the System's best
    Cryptographic random source.

    Use GMP if available, otherwise fallback to PyCrypto
    """
    randfunc = random.SystemRandom()

    if HAVE_GMP:        
        r = gmpy2.mpz(randfunc.getrandbits(N))
        r = gmpy2.bit_set(r, N - 1)
        return int(gmpy2.next_prime(r))    
    else:        
        n = randfunc.randrange(2**(N-1), 2**N) | 1
        while not is_prime(n):
            n += 2
        return n


def generate_pq(n_length):
    """Return p and q
    """
    
    p = q = n = None
    n_len = 0
    while n_len != n_length:
        p = getprimeover(n_length // 2)
        q = p
        while q == p:
            q = getprimeover(n_length // 2)
        n = p * q
        n_len = n.bit_length()
    
    return p,q


def int_to_bytes(num):
    return bytes([num])


def save_bytes_to_file(data,path):
    f = open(path, 'wb+')
    f.write(data)
    f.close()


key_sizes = [64, 128, 256, 512, 1024, 2048]
snap_dir = './generated_keys/'

if not osp.exists(snap_dir):
    os.makedirs(snap_dir)

for key_size in key_sizes:

    p,q = generate_pq(key_size)

    print('p',p)
    print(hex(p))
    print('q',q)
    print(hex(q))

    data = p.to_bytes(key_size//16,'little') + q.to_bytes(key_size//16,'little')
    pq_file = osp.join(snap_dir,'key_{}'.format(key_size))
    save_bytes_to_file(data,pq_file)

    f = open(pq_file, 'rb')
    data = f.read()
    p = int.from_bytes(data[:key_size//16],'little')
    q = int.from_bytes(data[key_size//16:],'little')

    print('p',p)
    print(hex(p))
    print('q',q)
    print(hex(q))

    print('')

    