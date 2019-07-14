/***

Copyright (c) 2018-2019, NVIDIA CORPORATION.  All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

***/


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <cuda.h>
#include <gmp.h>
#include <sys/time.h>
#include "time.h"
#include "cgbn/cgbn.h"
//#include "../utility/cpu_support.h"
//#include "../utility/cpu_simple_bn_math.h"
#include "../utility/gpu_support.h"
#include <curand.h>
#include <curand_kernel.h>

/************************************************************************************************
 *  This example performs component-wise addition of two arrays of 1024-bit bignums.
 *
 *  The example uses a number of utility functions and macros:
 *
 *    random_words(uint32_t *words, uint32_t count)
 *       fills words[0 .. count-1] with random data
 *
 *    add_words(uint32_t *r, uint32_t *a, uint32_t *b, uint32_t count) 
 *       sets bignums r = a+b, where r, a, and b are count words in length
 *
 *    compare_words(uint32_t *a, uint32_t *b, uint32_t count)
 *       compare bignums a and b, where a and b are count words in length.
 *       return 1 if a>b, 0 if a==b, and -1 if b>a
 *    
 *    CUDA_CHECK(call) is a macro that checks a CUDA result for an error,
 *    if an error is present, it prints out the error, call, file and line.
 *
 *    CGBN_CHECK(report) is a macro that checks if a CGBN error has occurred.
 *    if so, it prints out the error, and instance information
 *
 ************************************************************************************************/
 
// IMPORTANT:  DO NOT DEFINE TPI OR BITS BEFORE INCLUDING CGBN
#define TPI 32
#define BITS 2048
//#define BITS 2048  
#define INSTANCES 100000
#define MAX_RAND_SEED 4294967295U
#define WINDOW_BITS 5
#include <time.h> 
/*
typedef struct {
  cgbn_mem_t<BITS> cgbn;
  mpz_t mpz;
} mixed_mpz_t;
*/

// helpful typedefs for the kernel
typedef cgbn_context_t<TPI>         context_t;
typedef cgbn_env_t<context_t, BITS> env_t;
typedef cgbn_mem_t<BITS> gpu_mpz; 

void store2dev(cgbn_mem_t<BITS> *address,  mpz_t z) {
  size_t words;
  if(mpz_sizeinbase(z, 2)>BITS) {
    printf("error mpz_sizeinbase:%d\n", mpz_sizeinbase(z, 2));
    exit(1);
  }

  mpz_export((uint32_t *)address, &words, -1, sizeof(uint32_t), 0, 0, z);
  while(words<(BITS+31)/32)
    ((uint32_t *)address)[words++]=0;
}

void store2gmp(mpz_t z, cgbn_mem_t<BITS> *address ) {
  mpz_import(z, (BITS+31)/32, -1, sizeof(uint32_t), 0, 0, (uint32_t *)address);
}


void getprimeover(mpz_t rop, int bits, int &seed_start){
  gmp_randstate_t state;
  gmp_randinit_default(state);
  gmp_randseed_ui(state, seed_start);
  seed_start++;
  mpz_t rand_num;
  mpz_init(rand_num);
  mpz_urandomb(rand_num, state, bits);
  gmp_printf("rand_num:%Zd\n", rand_num);
  mpz_setbit(rand_num, bits-1);
  mpz_nextprime(rop, rand_num); 
  mpz_clear(rand_num);
}

void invert(mpz_t rop, mpz_t a, mpz_t b) {
  mpz_invert(rop, a, b);
}

class PaillierPublicKey {
 public:
  cgbn_mem_t<BITS> g;
  cgbn_mem_t<BITS> n;
  cgbn_mem_t<BITS> nsquare;
  cgbn_mem_t<BITS> max_int;
  void init(mpz_t &n, mpz_t g) {
    mpz_t nsquare, max_int; 
    mpz_init(nsquare);
    mpz_init(max_int);
    mpz_add_ui(g, n,1); 
    mpz_mul(nsquare, n, n);
    mpz_div_ui(max_int, n, 3);
    mpz_sub_ui(max_int, max_int, 1);
    store2dev(&this->g, g); 
    store2dev(&this->n, n); 
    store2dev(&this->nsquare, nsquare); 
    store2dev(&this->max_int, max_int); 
    mpz_clear(nsquare);
    mpz_clear(max_int);
  }
};

class PaillierPrivateKey {
 public:
  cgbn_mem_t<BITS> p;
  cgbn_mem_t<BITS> q;
  cgbn_mem_t<BITS> psquare;
  cgbn_mem_t<BITS> qsquare;
  cgbn_mem_t<BITS> q_inverse;
  cgbn_mem_t<BITS> hp;
  cgbn_mem_t<BITS> hq;

  void h_func_gmp(mpz_t rop, mpz_t g, mpz_t x, mpz_t xsquare) {
    mpz_t tmp;
    mpz_init(tmp);
    mpz_sub_ui(tmp, x, 1);
    mpz_powm(rop, g, tmp, xsquare); 
    mpz_sub_ui(rop, rop, 1);
    mpz_div(rop, rop, x);
    invert(rop, rop, x);
    mpz_clear(tmp); 
  }
  void init(PaillierPublicKey pub_key, mpz_t g, mpz_t raw_p, mpz_t raw_q) {
    // TODO: valid publick key    
    mpz_t p, q, psquare, qsquare, q_inverse, hp, hq;
    mpz_init(p);
    mpz_init(q);
    mpz_init(psquare);
    mpz_init(qsquare);
    mpz_init(q_inverse);
    mpz_init(hp);
    mpz_init(hq);
    if(mpz_cmp(raw_q, raw_p) < 0) {
      mpz_set(p, raw_q);
      mpz_set(q, raw_p);
    } else {
      mpz_set(p, raw_p);
      mpz_set(q, raw_q);
    }
    mpz_mul(psquare, p, p);
    mpz_mul(qsquare, q, q);
    invert(q_inverse, q, p);
    h_func_gmp(hp, g, p, psquare); 
    h_func_gmp(hq, g, q, qsquare); 

    gmp_printf("hp:%Zd\n", hp);
    gmp_printf("hq:%Zd\n", hq);
    store2dev(&this->p, p);
    store2dev(&this->q, q);
    store2dev(&this->psquare, psquare);
    store2dev(&this->qsquare, qsquare);
    store2dev(&this->q_inverse, q_inverse);
    store2dev(&this->hp, hp);
    store2dev(&this->hq, hq);

    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(psquare);
    mpz_clear(qsquare);
    mpz_clear(q_inverse);
    mpz_clear(hp);
    mpz_clear(hq);
  }
};

__constant__  PaillierPublicKey gpu_pub_key[1];
__constant__  PaillierPrivateKey gpu_priv_key[1]; 

//code for CGBN sample
__device__ __forceinline__ void fixed_window_powm_odd(env_t &bn_env, env_t::cgbn_t &result, const env_t::cgbn_t &x, const env_t::cgbn_t &power, const env_t::cgbn_t &modulus) {
  env_t::cgbn_t       t;
  env_t::cgbn_local_t window[1<<WINDOW_BITS];
  int32_t    index, position, offset;
  uint32_t   np0;

  // conmpute x^power mod modulus, using the fixed window algorithm
  // requires:  x<modulus,  modulus is odd

  // compute x^0 (in Montgomery space, this is just 2^BITS - modulus)
  cgbn_negate(bn_env, t, modulus);
  cgbn_store(bn_env, window+0, t);
  
  // convert x into Montgomery space, store into window table
  np0=cgbn_bn2mont(bn_env, result, x, modulus);
  cgbn_store(bn_env, window+1, result);
  cgbn_set(bn_env, t, result);
  
  // compute x^2, x^3, ... x^(2^WINDOW_BITS-1), store into window table
  #pragma nounroll
  for(index=2;index<(1<<WINDOW_BITS);index++) {
    cgbn_mont_mul(bn_env, result, result, t, modulus, np0);
    cgbn_store(bn_env, window+index, result);
  }

  // find leading high bit
  position= BITS - cgbn_clz(bn_env, power);

  // break the exponent into chunks, each WINDOW_BITS in length
  // load the most significant non-zero exponent chunk
  offset=position % WINDOW_BITS;
  if(offset==0)
    position=position-WINDOW_BITS;
  else
    position=position-offset;
  index=cgbn_extract_bits_ui32(bn_env, power, position, WINDOW_BITS);
  cgbn_load(bn_env, result, window+index);
  
  // process the remaining exponent chunks
  while(position>0) {
    // square the result WINDOW_BITS times
    #pragma nounroll
    for(int sqr_count=0;sqr_count<WINDOW_BITS;sqr_count++)
      cgbn_mont_sqr(bn_env, result, result, modulus, np0);
    
    // multiply by next exponent chunk
    position=position-WINDOW_BITS;
    index=cgbn_extract_bits_ui32(bn_env, power, position, WINDOW_BITS);
    cgbn_load(bn_env, t, window+index);
    cgbn_mont_mul(bn_env, result, result, t, modulus, np0);
  }
  
  // we've processed the exponent now, convert back to normal space
  cgbn_mont2bn(bn_env, result, result, modulus, np0);
}

__device__ __forceinline__ void sliding_window_powm_odd(env_t &bn_env, env_t::cgbn_t &result, const env_t::cgbn_t &x, const env_t::cgbn_t &power, const env_t::cgbn_t &modulus) {
  env_t::cgbn_t         t, starts;
  int32_t      index, position, leading;
  uint32_t     mont_inv;
  env_t::cgbn_local_t   odd_powers[1<<WINDOW_BITS-1];

  // conmpute x^power mod modulus, using Constant Length Non-Zero windows (CLNZ).
  // requires:  x<modulus,  modulus is odd
      
  // find the leading one in the power
  leading=BITS-1-cgbn_clz(bn_env, power);
  if(leading>=0) {
    // convert x into Montgomery space, store in the odd powers table
    mont_inv=cgbn_bn2mont(bn_env, result, x, modulus);
    
    // compute t=x^2 mod modulus
    cgbn_mont_sqr(bn_env, t, result, modulus, mont_inv);
    
    // compute odd powers window table: x^1, x^3, x^5, ...
    cgbn_store(bn_env, odd_powers, result);
    #pragma nounroll
    for(index=1;index<(1<<WINDOW_BITS-1);index++) {
      cgbn_mont_mul(bn_env, result, result, t, modulus, mont_inv);
      cgbn_store(bn_env, odd_powers+index, result);
    }

    // starts contains an array of bits indicating the start of a window
    cgbn_set_ui32(bn_env, starts, 0);

    // organize p as a sequence of odd window indexes
    position=0;
    while(true) {
      if(cgbn_extract_bits_ui32(bn_env, power, position, 1)==0)
        position++;
      else {
        cgbn_insert_bits_ui32(bn_env, starts, starts, position, 1, 1);
        if(position+WINDOW_BITS>leading)
          break;
        position=position+WINDOW_BITS;
      }
    }

    // load first window.  Note, since the window index must be odd, we have to
    // divide it by two before indexing the window table.  Instead, we just don't
    // load the index LSB from power
    index=cgbn_extract_bits_ui32(bn_env, power, position+1, WINDOW_BITS-1);
    cgbn_load(bn_env, result, odd_powers+index);
    position--;
    
    // Process remaining windows 
    while(position>=0) {
      cgbn_mont_sqr(bn_env, result, result, modulus, mont_inv);
      if(cgbn_extract_bits_ui32(bn_env, starts, position, 1)==1) {
        // found a window, load the index
        index=cgbn_extract_bits_ui32(bn_env, power, position+1, WINDOW_BITS-1);
        cgbn_load(bn_env, t, odd_powers+index);
        cgbn_mont_mul(bn_env, result, result, t, modulus, mont_inv);
      }
      position--;
    }
    
    // convert result from Montgomery space
    cgbn_mont2bn(bn_env, result, result, modulus, mont_inv);
  }
  else {
    // p=0, thus x^p mod modulus=1
    cgbn_set_ui32(bn_env, result, 1);
  }
}



__device__  __forceinline__ void l_func(env_t &bn_env, env_t::cgbn_t &out, env_t::cgbn_t &cipher_t, env_t::cgbn_t &x_t, env_t::cgbn_t &xsquare_t, env_t::cgbn_t &hx_t) {
  env_t::cgbn_t  tmp, tmp2, cipher_lt;
  env_t::cgbn_wide_t  tmp_wide;
  cgbn_sub_ui32(bn_env, tmp2, x_t, 1);  
  if(cgbn_compare(bn_env, cipher_t, xsquare_t) >= 0) {
    cgbn_rem(bn_env, cipher_lt, cipher_t, xsquare_t);
    //cgbn_modular_power(bn_env, tmp, cipher_lt, tmp2, xsquare_t);
    sliding_window_powm_odd(bn_env,tmp,cipher_lt,tmp2,xsquare_t);
  } else {
    //cgbn_modular_power(bn_env, tmp, cipher_t, tmp2, xsquare_t);
    sliding_window_powm_odd(bn_env,tmp,cipher_t,tmp2,xsquare_t);
  }
  cgbn_sub_ui32(bn_env, tmp, tmp, 1);
  cgbn_div(bn_env, tmp, tmp, x_t);
  cgbn_mul_wide(bn_env, tmp_wide, tmp, hx_t);
  cgbn_rem_wide(bn_env, tmp, tmp_wide, x_t);
  cgbn_set(bn_env, out, tmp);
}

__device__ __forceinline__ void powmod(env_t &bn_env, env_t::cgbn_t &r, env_t::cgbn_t &a, env_t::cgbn_t &b, env_t::cgbn_t &c) {
  if(cgbn_compare(bn_env, a, b) >= 0) {
    cgbn_rem(bn_env, r, a, c);
  } 
  cgbn_modular_power(bn_env, r, r, b, c);
}


__global__ void setup_kernel(curandState *state){
  int idx = threadIdx.x+blockDim.x*blockIdx.x;
  curand_init(1234, idx, 0, &state[idx]);
}

   // cuda rand ?
  // reuse obfuscated random value ?
__global__ __noinline__ void apply_obfuscator(cgbn_error_report_t *report, gpu_mpz *ciphers, gpu_mpz *obfuscators, int count, curandState *state) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  int tid= idx/TPI;
  if(tid>=count)
    return;

  context_t      bn_context(cgbn_report_monitor, report, tid);  
  env_t          bn_env(bn_context.env<env_t>());                   
  env_t::cgbn_t  n, nsquare,cipher, r, tmp; 
  env_t::cgbn_wide_t tmp_wide;     

  curandState localState = state[idx];
  unsigned int rand_seed = curand_uniform(&localState) * MAX_RAND_SEED;                  
  state[idx] = localState;
  cgbn_set_ui32(bn_env, r, rand_seed); // TODO: new rand or reuse
  cgbn_load(bn_env, n, &gpu_pub_key[0].n);      
  cgbn_load(bn_env, nsquare, &gpu_pub_key[0].nsquare);
  cgbn_load(bn_env, cipher, &ciphers[tid]);
  cgbn_modular_power(bn_env,tmp, r, n, nsquare); 
  cgbn_mul_wide(bn_env, tmp_wide, cipher, tmp); 
  cgbn_rem_wide(bn_env, r, tmp_wide, nsquare); 
  cgbn_store(bn_env, obfuscators + tid, r);   // store r into sum
   
}
__global__ __noinline__ void raw_encrypt(cgbn_error_report_t *report, gpu_mpz *plains, gpu_mpz *ciphers,int count, int rand_seed ) {
  int tid=(blockIdx.x*blockDim.x + threadIdx.x)/TPI;
  if(tid>=count)
    return;
  context_t      bn_context(cgbn_report_monitor, report, tid);  
  env_t          bn_env(bn_context.env<env_t>());                   
  env_t::cgbn_t  n, nsquare, plain,  tmp, max_int, neg_plain, neg_cipher, cipher;               
  cgbn_load(bn_env, n, &gpu_pub_key[0].n);      
  cgbn_load(bn_env, plain, plains + tid);
  cgbn_load(bn_env, nsquare, &gpu_pub_key[0].nsquare);
  cgbn_load(bn_env, max_int, &gpu_pub_key[0].max_int);
  cgbn_load(bn_env, plain, plains + tid);
  cgbn_sub(bn_env, tmp, n, max_int); 
  if(cgbn_compare(bn_env, plain, tmp) >= 0 &&  cgbn_compare(bn_env, plain, n) < 0) {
    // Very large plaintext, take a sneaky shortcut using inverses
    cgbn_sub(bn_env, neg_plain, n, plain);
    cgbn_mul(bn_env, neg_cipher, n, neg_plain);
    cgbn_add_ui32(bn_env, neg_cipher, neg_cipher, 1);
    cgbn_rem(bn_env, neg_cipher, neg_cipher, nsquare);
    cgbn_modular_inverse(bn_env, cipher, neg_cipher, nsquare);
  } else {
    cgbn_mul(bn_env, cipher, n, plain);
    cgbn_add_ui32(bn_env, cipher, cipher, 1);
    cgbn_rem(bn_env, cipher, cipher, nsquare);
  }

  cgbn_store(bn_env, ciphers + tid, cipher);   // store r into sum

}

__global__ __noinline__ void raw_encrypt_with_obfs(cgbn_error_report_t *report, gpu_mpz *plains, gpu_mpz *ciphers,int count, curandState *state) {
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  int tid= idx/TPI;
  if(tid>=count)
    return;
  context_t      bn_context(cgbn_report_monitor, report, tid);  
  env_t          bn_env(bn_context.env<env_t>());                   
  env_t::cgbn_t  n, nsquare, plain,  tmp, max_int, neg_plain, neg_cipher, cipher;               
  cgbn_load(bn_env, n, &gpu_pub_key[0].n);      
  cgbn_load(bn_env, plain, plains + tid);
  cgbn_load(bn_env, nsquare, &gpu_pub_key[0].nsquare);
  cgbn_load(bn_env, max_int, &gpu_pub_key[0].max_int);
  cgbn_load(bn_env, plain, plains + tid);
  cgbn_sub(bn_env, tmp, n, max_int); 
  if(cgbn_compare(bn_env, plain, tmp) >= 0 &&  cgbn_compare(bn_env, plain, n) < 0) {
    // Very large plaintext, take a sneaky shortcut using inverses
    cgbn_sub(bn_env, neg_plain, n, plain);
    cgbn_mul(bn_env, neg_cipher, n, neg_plain);
    cgbn_add_ui32(bn_env, neg_cipher, neg_cipher, 1);
    cgbn_rem(bn_env, neg_cipher, neg_cipher, nsquare);
    cgbn_modular_inverse(bn_env, cipher, neg_cipher, nsquare);
  } else {
    cgbn_mul(bn_env, cipher, n, plain);
    cgbn_add_ui32(bn_env, cipher, cipher, 1);
    cgbn_rem(bn_env, cipher, cipher, nsquare);
  }

  env_t::cgbn_t r; 
  env_t::cgbn_wide_t tmp_wide;     

  curandState localState = state[idx];
  unsigned int rand_seed = curand_uniform(&localState) * MAX_RAND_SEED;                  
  state[idx] = localState;

  cgbn_set_ui32(bn_env, r, rand_seed); // TODO: new rand or reuse
  
  cgbn_modular_power(bn_env,tmp, r, n, nsquare); 
  cgbn_mul_wide(bn_env, tmp_wide, cipher, tmp); 
  cgbn_rem_wide(bn_env, r, tmp_wide, nsquare); 
  cgbn_store(bn_env, ciphers + tid, r);   // store r into sum
}
 
__global__ __noinline__ void raw_add(cgbn_error_report_t *report, gpu_mpz *ciphers_r, gpu_mpz *ciphers_a, gpu_mpz *ciphers_b,int count) {
  int tid=(blockIdx.x*blockDim.x + threadIdx.x)/TPI;
  if(tid>=count)
    return;
  context_t      bn_context(cgbn_report_monitor, report, tid);  
  env_t          bn_env(bn_context.env<env_t>());                   
  env_t::cgbn_t  nsquare, r, a, b;   
  env_t::cgbn_wide_t  r_wide;               
  cgbn_load(bn_env, nsquare, &gpu_pub_key[0].nsquare);      
  cgbn_load(bn_env, a, ciphers_a + tid);      
  cgbn_load(bn_env, b, ciphers_b + tid);
  cgbn_mul_wide(bn_env, r_wide, a, b);
  cgbn_rem_wide(bn_env, r, r_wide, nsquare);

/*    
 uint32_t np0;

// convert a and b to Montgomery space
np0=cgbn_bn2mont(bn_env, a, a, nsquare);
cgbn_bn2mont(bn_env, b, b, nsquare);

cgbn_mont_mul(bn_env, r, a, b, nsquare, np0);

// convert r back to normal space
cgbn_mont2bn(bn_env, r, r, nsquare, np0);
*/
  cgbn_store(bn_env, ciphers_r + tid, r);
}

__global__ void raw_mul(cgbn_error_report_t *report, gpu_mpz *ciphers_r, gpu_mpz *ciphers_a, gpu_mpz *plains_b,int count) {
  int tid=(blockIdx.x*blockDim.x + threadIdx.x)/TPI;
  if(tid>=count)
    return;
  context_t      bn_context(cgbn_report_monitor, report, tid);  
  env_t          bn_env(bn_context.env<env_t>());                   
  env_t::cgbn_t  n,max_int, nsquare, r, cipher, plain, neg_c, neg_scalar,tmp;               

  cgbn_load(bn_env, n, &gpu_pub_key[0].n);      
  cgbn_load(bn_env, max_int, &gpu_pub_key[0].max_int);      
  cgbn_load(bn_env, nsquare, &gpu_pub_key[0].nsquare);      
  cgbn_load(bn_env, cipher, ciphers_a + tid);      
  cgbn_load(bn_env, plain, plains_b + tid);

  cgbn_sub(bn_env, tmp, n, max_int); 
 if(cgbn_compare(bn_env, plain, tmp) >= 0 ) {
    // Very large plaintext, take a sneaky shortcut using inverses
    cgbn_modular_inverse(bn_env,neg_c, cipher, nsquare);
    cgbn_sub(bn_env, neg_scalar, n, plain);
    sliding_window_powm_odd(bn_env, r, neg_c, neg_scalar, nsquare);
  } else {
    sliding_window_powm_odd(bn_env, r, cipher, plain, nsquare); 
  }

  cgbn_store(bn_env, ciphers_r + tid, r);
}

  
__global__ void raw_decrypt(cgbn_error_report_t *report, gpu_mpz *plains, gpu_mpz *ciphers, int count) {
  int tid=(blockIdx.x*blockDim.x + threadIdx.x)/TPI;
  if(tid>=count)
    return;
  context_t      bn_context(cgbn_report_monitor, report, tid);
  env_t          bn_env(bn_context.env<env_t>());
  env_t::cgbn_t  mp, mq, tmp, q_inverse, n, p, q, hp, hq, psquare, qsquare, cipher;  
  cgbn_load(bn_env, cipher, ciphers + tid);
  cgbn_load(bn_env, q_inverse, &gpu_priv_key[0].q_inverse);
  cgbn_load(bn_env, n, &gpu_pub_key[0].n);
  cgbn_load(bn_env, p, &gpu_priv_key[0].p);
  cgbn_load(bn_env, q, &gpu_priv_key[0].q);
  cgbn_load(bn_env, hp, &gpu_priv_key[0].hp);
  cgbn_load(bn_env, hq, &gpu_priv_key[0].hq);
  cgbn_load(bn_env, psquare, &gpu_priv_key[0].psquare);
  cgbn_load(bn_env, qsquare, &gpu_priv_key[0].qsquare);
  l_func(bn_env, mp, cipher, p, psquare, hp); 
  l_func(bn_env, mq, cipher, q, qsquare, hq); 
  cgbn_sub(bn_env, tmp, mp, mq);
  cgbn_mul(bn_env, tmp, tmp, q_inverse); 
  cgbn_rem(bn_env, tmp, tmp, p);
  cgbn_mul(bn_env, tmp, tmp, q);
  cgbn_add(bn_env, tmp, mq, tmp);
  cgbn_rem(bn_env, tmp, tmp, n);
  cgbn_store(bn_env, plains + tid, tmp);
} 

void generate_keypair(PaillierPublicKey &pub_key, PaillierPrivateKey &priv_key) {
  mpz_t p;
  mpz_t q;    
  mpz_t n;    
  mpz_init(p);
  mpz_init(q);
  mpz_init(n);
  int n_len = 0;
  srand((unsigned)time(NULL));
  //int seed_start = rand();
  int seed_start = 2;
  int key_len = BITS/2;
  while(n_len != key_len) {
    getprimeover(p, key_len / 2, seed_start);
    mpz_set(q, p);
    while(mpz_cmp(p, q) == 0){
      getprimeover(q, key_len / 2, seed_start);
      mpz_mul(n, p, q);
      n_len = mpz_sizeinbase(n, 2);
    }
  }
  
  mpz_t g;
  mpz_init(g);
  printf("rand bits2:%d\n",mpz_sizeinbase(n, 2));
  pub_key.init(n, g);
  priv_key.init(pub_key,g, p, q);
  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(n);
  mpz_clear(g);
}

long mil_time(){
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*1000000 + tv.tv_usec;
}

void log_time(const char* label, long &start){
  printf("%s --time--: %ld\n",label, mil_time() - start);
  start = mil_time();
}


void log_time(const char* label,int count, long &start){
  long op_time = mil_time() - start;
  double ops = (double(count) / double(op_time)) * 1000.0;
  printf("%s k ops/s: %.4f, time: %ld\n",label, ops, op_time);
  start = mil_time();
}



int main(int argc,char **argv) {
  int32_t              TPB=128;
  int32_t              IPB=TPB/TPI;
  long start;
  start = mil_time();
  cgbn_error_report_t *report;
  PaillierPrivateKey priv_key;
  PaillierPublicKey pub_key;
 
  generate_keypair(pub_key, priv_key);
  int count = 1000*100;
  if(argc > 1){
    count = atoi(argv[1]);
  }
  //int count = 10;
  int mem_size = sizeof(gpu_mpz) * count;
  gpu_mpz *plains = (gpu_mpz*)malloc(mem_size); 
  gpu_mpz *plains2 = (gpu_mpz*)malloc(mem_size); 
  gpu_mpz *ciphers = (gpu_mpz*)malloc(mem_size); 
  gpu_mpz *obfs = (gpu_mpz*)malloc(mem_size); 
  gpu_mpz *gpu_plains;
  gpu_mpz *gpu_plains2;
  gpu_mpz *gpu_ciphers;
  curandState *d_state;
  
  int block_size = (count+IPB-1)/IPB;
  int thread_size = TPB;

  cudaSetDevice(0);

  cudaMalloc(&d_state, sizeof(curandState) * block_size * thread_size);
  setup_kernel<<<block_size, thread_size>>>(d_state);

  cudaMalloc((void **)&gpu_plains, mem_size); 
  cudaMalloc((void **)&gpu_plains2, mem_size); 
  cudaMalloc((void **)&gpu_ciphers, mem_size);

  cudaMemcpyToSymbol(gpu_priv_key, &priv_key, sizeof(priv_key)); 
  cudaMemcpyToSymbol(gpu_pub_key, &pub_key, sizeof(pub_key)); 
  
  CUDA_CHECK(cgbn_error_report_alloc(&report));
  for(int i = 0; i < count; i++){
    mpz_t n;
    mpz_init(n);
    mpz_set_ui(n, i);
    store2dev(plains + i, n);  
    //gmp_printf("input:%Zd\n", n);  
    mpz_clear(n); 
  }
  log_time("prepare", start);
  cudaMemcpy(gpu_plains, plains, mem_size, cudaMemcpyHostToDevice); 

  log_time("time1", start);
  raw_encrypt_with_obfs<<<block_size, thread_size>>>(report,  gpu_plains, gpu_ciphers, count,d_state); 
  CUDA_CHECK(cudaDeviceSynchronize());
  log_time("enc time",count, start);
  
  cudaMemcpy(ciphers, gpu_ciphers, mem_size, cudaMemcpyDeviceToHost);   
  for(int i = 0; i <  count; i++){
    mpz_t n;
    mpz_init(n);
    store2gmp(n, ciphers + i);
    mpz_clear(n); 
  }

  log_time("time2", start);
  raw_decrypt<<<block_size, thread_size>>>(report, gpu_plains2, gpu_ciphers, count);
  CUDA_CHECK(cudaDeviceSynchronize());
  log_time("dec time",count, start);


  cudaMemcpy(plains2, gpu_plains2, mem_size, cudaMemcpyDeviceToHost); 
  for(int i = 0; i < count; i++){
    mpz_t n;
    mpz_init(n);
    store2gmp(n, plains2 + i);
    //gmp_printf("output:%Zd\n", n);
    mpz_clear(n); 
  }


  log_time("time3", start);
  raw_add<<<block_size, thread_size>>>(report, gpu_plains, gpu_ciphers, gpu_ciphers, count);
  CUDA_CHECK(cudaDeviceSynchronize());
  log_time("add time",count, start);


  raw_decrypt<<<block_size, thread_size>>>(report, gpu_plains2, gpu_plains, count);
  CUDA_CHECK(cudaDeviceSynchronize());
  cudaMemcpy(plains2, gpu_plains2, mem_size, cudaMemcpyDeviceToHost); 
  for(int i = 0; i < count; i++){
    mpz_t n;
    mpz_init(n);
    store2gmp(n, plains2 + i);
    //gmp_printf("add output:%Zd\n", n);
    mpz_clear(n); 
  }


  log_time("time4", start);
  raw_mul<<<block_size, thread_size>>>(report, gpu_plains, gpu_ciphers, gpu_plains2, count);
  CUDA_CHECK(cudaDeviceSynchronize());
  log_time("mul time",count, start);

  raw_decrypt<<<block_size, thread_size>>>(report, gpu_plains2, gpu_plains, count);
  CUDA_CHECK(cudaDeviceSynchronize());
  cudaMemcpy(plains2, gpu_plains2, mem_size, cudaMemcpyDeviceToHost); 
  for(int i = 0; i < count; i++){
    mpz_t n;
    mpz_init(n);
    store2gmp(n, plains2 + i);
    //gmp_printf("mul output:%Zd\n", n);
    mpz_clear(n); 
  }

  CGBN_CHECK(report);
  CUDA_CHECK(cgbn_error_report_free(report));

  free(plains);
  free(plains2);
  free(ciphers);
  cudaFree(gpu_plains);
  cudaFree(gpu_plains2);
  cudaFree(gpu_ciphers);

}

