"""
Benchmark key generation, encryption and decryption.

"""

import random
import resource
import time
import phe.paillier as paillier
import gmpy2
import sys


def generate_paillier_keypair(p,q):
    n = p*q

    public_key = paillier.PaillierPublicKey(n)
    private_key = paillier.PaillierPrivateKey(public_key, p, q)

    return public_key, private_key


def bench_encrypt(pk, nums, r_values, res_enc):
    for idx in range(len(nums)):
        res_enc.append(pk.raw_encrypt(nums[idx], r_value=r_values[idx]))


def bench_decrypt(vk, nums,res_dec):
    for idx in range(len(nums)):
        res_dec.append(vk.raw_decrypt(nums[idx]))


def bench_encrypt_add(pk,inputs_enc1,inputs_enc2,res_enc):
    for input_enc1,input_enc2 in zip(inputs_enc1,inputs_enc2):
        res_enc.append( (input_enc1*input_enc2) % pk.nsquare)


def bench_add(nums1, nums2,res):
    for num1, num2 in zip(nums1, nums2):
        res.append(num1 + num2)


def bench_mul(nums1, nums2):
    for num1, num2 in zip(nums1, nums2):
        num1 * num2


def time_method(method, *args):
    start = time.time()
    method(*args)
    return time.time() - start


def bench_time(test_size,key_bytes,pk,vk, data_size=4):

    print('Paillier Benchmarks with key size of {} bits'.format(key_bytes*8))

    input1 = []
    input2 = []
    input3 = []
    input4 = []
    pks = [] 
    vks  = []
    r_values1 = []
    r_values2 = []
    inputs1 = [] 
    inputs2 = [] 

    half_key_bytes = key_bytes//2


    
    
    for i in range(test_size*data_size):
        input1.append((i*17+11)%256)   

    for i in range(test_size):
        inputs1.append(int.from_bytes(bytes(input1[i*data_size:(i+1)*data_size]),'little'))
    
    for i in range(test_size*half_key_bytes):
        input2.append((i*13+11)%256)
    
    for i in range(test_size):
        r_values1.append(int.from_bytes(bytes(input2[i*half_key_bytes:(i+1)*half_key_bytes]),'little'))


    
    
    for i in range(test_size*data_size):
        input3.append((i*23+11)%256)   

    for i in range(test_size):
        inputs2.append(int.from_bytes(bytes(input3[i*data_size:(i+1)*data_size]),'little'))

    for i in range(test_size*half_key_bytes):
        input4.append((i*19+11)%256)
    
    for i in range(test_size):
        r_values2.append(int.from_bytes(bytes(input4[i*half_key_bytes:(i+1)*half_key_bytes]),'little'))
    
    
    '''print('input1: ', hex(inputs1[0]))
    print('input2: ', hex(inputs2[0]))        
    print('rand: ', hex(r_values[0]))    
    print('raw encrypted1: ', hex(pk.raw_encrypt(inputs[0],r_value = r_values[0])))
    print('raw encrypted2: ', hex(pk.raw_encrypt(inputs2[0],r_value = r_values[0])))
    print('encrypted sum: ', hex(pk.raw_encrypt(inputs[0],r_value = r_values[0]) + pk.raw_encrypt(inputs2[0],r_value = r_values[0])))
    print('decrypted: ', hex(vk.raw_decrypt(pk.raw_encrypt(inputs[0],r_value = r_values[0]) + pk.raw_encrypt(inputs2[0],r_value = r_values[0]))))
    print('for check: ', hex(inputs[0] + inputs2[0]))'''
    
    res_enc1 = []
    res_enc2 = []
    res_enc = []
    res_dec = []
    for_check = []
    
    ones = [1.0 for _ in range(test_size)]

    times = [
        time_method(bench_encrypt, pk, inputs1, r_values1, res_enc),
        time_method(bench_decrypt, vk, res_enc, res_dec)
        #time_method(bench_encrypt, pk, inputs1, r_values1, res_enc1),
        #time_method(bench_encrypt, pk, inputs2, r_values2, res_enc2),
        #time_method(bench_encrypt_add, pk, res_enc1, res_enc2, res_enc),
        #time_method(bench_decrypt, vk, res_enc, res_dec),
        #time_method(bench_add, inputs1,inputs2,for_check)'''
    ]
    times = [t / test_size for t in times]
    ops = [int(1.0 / t) for t in times]

    print(
        'operation: time in seconds (# operations per second)\n'
        'encrypt: {:.6f} s ({} ops/s)\n'
        'decrypt: {:.6f} s ({} ops/s)\n'.format(
            times[0], ops[0], times[1], ops[1]
        )
    )

    '''print(res_enc1)
    print(res_enc2)
    print(res_enc)
    print(res_dec)
    print(for_check)'''


    return times

if __name__ == '__main__':    
    test_size = int(sys.argv[1])

    #key_sizes = [64, 128,256,512,1024,2048]
    key_sizes = [1024]
    

    pq_files = {
        64:'./generated_keys/key_64',
        128:'./generated_keys/key_128',
        256:'./generated_keys/key_256',
        512:'./generated_keys/key_512',
        1024:'./generated_keys/key_1024',
        2048:'./generated_keys/key_2048'
    }
    
    #key_sizes = [16, 32, 64, 128, 256]

    
    for key_size in key_sizes:
        key_bytes = key_size//8
        f = open(pq_files[key_size], 'rb')
        data = f.read()
        p = int.from_bytes(data[:key_size//16],'little')
        q = int.from_bytes(data[key_size//16:],'little')

        print('p: ',hex(p))
        print('q: ',hex(q))
        pk, vk = generate_paillier_keypair(p,q)
        bench_time(test_size,key_bytes,pk,vk,data_size=4)

