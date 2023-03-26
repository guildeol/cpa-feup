#include <iostream>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstring>

#include "stdint.h"

#include "papi_helper.h"

using namespace std;

#ifndef _OPENMP
    #error "Cannot build file: no OpenMP support"
#endif

// n >> 3   -> Division by 8 (gets the byte)
// n & 0x08 -> Remainder of division by 8 (gets the bit offset within the byte)
#define GET_BIT(_arr, _b)     (_arr[(_b) >> 3] & (1 << ((_b) & 0x07)))
#define CLEAR_BIT(_arr, _b)   (_arr[(_b) >> 3] &= ~(1 << ((_b) & 0x07)))

// Third implementation - we use the second implementation and reorganize computation so that the cache misses are reduced
// by searching several seed numbers in the same data block
void sieve_cache_friendly(uint64_t n, uint64_t block_size) 
{
    PapiHelper ph;
    uint64_t k = 3;
    uint32_t size = (n / 8) + 1;
    uint8_t *primes = new uint8_t[size];

    // Start by marking all values as primes
    memset(primes, 0xFF, size * sizeof(uint8_t));

    ph.startCounting();

    #pragma omp parallel
    {
        do
        {
            for (long long j = k*k ; j < n ; j += 2*k)
            {
                CLEAR_BIT(primes, j >> 1);
            }
            
            do
            {
                k += 2;
            } while (k*k <= n && !GET_BIT(primes, k >> 1));
            
        } while (k*k <= n);
    }

    ph.stopCounting();

    // cout << "2 ";
    // for (int i = 3; i < n; i += 2)
    //     if (GET_BIT(primes, i >> 1))
    //         cout << i << " ";

    // cout << endl;

    ph.report();

    delete[] primes;
}

int main (int argc, char *argv[])
{
    PapiHelper ph;
    uint64_t n;
    
    cout << "Power of 10: ";
    cin >> n;
 
    n = pow(10,n);
    
    cout << "THIRD IMPLEMENTATION";
    cout << endl;
    sieve_cache_friendly(n, n);
    cout << endl;
}