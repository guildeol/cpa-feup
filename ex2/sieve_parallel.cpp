#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <omp.h>

#include "stdint.h"

#include "papi_helper.h"

using namespace std;

#ifndef _OPENMP
    #error "Cannot build file: no OpenMP support"
#endif

// n >> 3   -> Division by 8 (gets the byte)
// n & 0x07 -> Remainder of division by 8 (gets the bit offset within the byte)
#define GET_BIT(_arr, _b)     (_arr[(_b) >> 3] & (1 << ((_b) & 0x07)))
#define CLEAR_BIT(_arr, _b)   (_arr[(_b) >> 3] &= ~(1 << ((_b) & 0x07)))

#define BLOCK_LOW(id, p, n)         ((id) * (n)/(p))
#define BLOCK_HIGH(id, p, n)        (BLOCK_LOW((id) + 1, (p), (n)) - 1)
#define BLOCK_SIZE(id, p, n)        (BLOCK_HIGH((id), (p), (n)) - BLOCK_LOW((id), (p), (n)))
#define BLOCK_OWNER(index, p, n)    (((p) * ((index) + 1) - 1) / (n))

uint8_t *primes;

// Third implementation - we use the second implementation and reorganize computation so that the cache misses are reduced
// by searching several seed numbers in the same data block
void sieve_parallel(uint64_t n, uint32_t threads) 
{
    PapiHelper ph;
    uint64_t k = 3;

    ph.startCounting();

    do
    {
        #pragma omp parallel num_threads(threads)
        {
            uint64_t first;
            uint64_t remainder;
            uint64_t id = omp_get_thread_num();
            uint64_t low_value = 2 + BLOCK_LOW(id, threads, n - 1);
            uint64_t high_value = 2 + BLOCK_HIGH(id, threads, n -1);
            uint64_t block_size = BLOCK_SIZE(id, threads, n - 1);

            // printf("ID: %lu\n", id);
            // printf("low_value: %lu\n", low_value);
            // printf("high_value: %lu\n", high_value);
            // printf("block_size: %lu\n", block_size);

            if (id == 0 && (block_size * block_size) < n)
            {
                cout << "Block size too small! Errors will happen!" << endl;
            }

            if ((k * k) > low_value)
            {
                first = (k * k) - low_value;
            }
            else
            {
                remainder = low_value % k;

                if (remainder == 0)
                    first = 0;
                else
                    first = k - remainder;
            }

            for (long long j = first; j < block_size; j += 2*k)
            {
                CLEAR_BIT(primes, (low_value + j) >> 1);
            }
        }
        
        do
        {
            k += 2;
        } while (k*k <= n && !GET_BIT(primes, k >> 1));
    } while (k*k <= n - 1);

    ph.stopCounting();

    ph.report();
}

int main (int argc, char *argv[])
{
    uint64_t n;
    uint32_t threads;
    
    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <pow_n> <threads>" << endl;
        return -1;
    }
    
    n = atoi(argv[1]);
    n = pow(10, n);

    threads = atoi(argv[2]);

    cout << "PARALLEL IMPLEMENTATION";
    cout << endl;

    uint32_t size = (n / 8) + 1;
    primes = new uint8_t[size];

    // Start by marking all values as primes
    memset(primes, 0xFF, size * sizeof(uint8_t));

    sieve_parallel(n, threads);

    // cout << "2 ";
    // for (int i = 3; i < n; i += 2)
    //     if (GET_BIT(primes, i >> 1))
    //         cout << i << " ";

    // cout << endl;

    delete[] primes;
}