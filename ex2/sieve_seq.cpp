#include <iostream>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstring>

#include "stdint.h"

#include "papi_helper.h"

using namespace std;

// n >> 3   -> Division by 8 (gets the byte)
// n & 0x08 -> Remainder of division by 8 (gets the bit offset within the byte)
#define GET_BIT(_arr, _b)     (_arr[(_b) >> 3] & (1 << ((_b) & 0x07)))
#define CLEAR_BIT(_arr, _b)   (_arr[(_b) >> 3] &= ~(1 << ((_b) & 0x07)))

// First implementation - we divide each element by k and by checking if the remainder of the division is zero
void sieve_remainder(uint64_t n)
{
    PapiHelper ph;
    uint64_t k = 3;
    bool *primes = new bool[n];

    // Start by marking all values as primes
    memset(primes, true, n * sizeof(bool));

    ph.startCounting();
    
    for (int k = 2; (k * k) < n; k++)
    {
        for (int j = (k * k); j < n; j++)
        {
            if (primes[j] == false)
                continue;

            if (j % k == 0)
            {
                primes[j] = false;

                for(int m = j; m < n; m += j)
                    primes[m] = false;
            }
        }
    }

    ph.stopCounting();

    // cout << "2 ";
    // for (int i = 3; i < n; i += 2)
    //     if (primes[i])
    //         cout << i << " ";
    
    // cout << endl;

    ph.report();

    delete[] primes;
}


// Second implementation - we use a fast marking where only the multiples are computed (2k, 3k, 4k, etc.)
void sieve_fast_marking(uint64_t n)
{
    PapiHelper ph;
    uint64_t k = 3;
    bool *primes = new bool[n];

    // Start by marking all values as primes
    memset(primes, true, n * sizeof(bool));

    ph.startCounting();

    do
    {
        for (long long j = k*k ; j<n ; j += 2*k)
        {
            primes[j >> 1] = false;
        }
        
        do
        {
            k += 2;
        }while (k*k <= n && !primes[k>>1]);
        
    } while (k*k <= n);
    
    ph.stopCounting();

    // cout << "2 ";
    // for (int i = 3; i < n; i += 2)
    //     if (primes[i >> 1])
    //         cout << i << " ";

    // cout << endl;

    ph.report();

    delete[] primes;
}

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
    
    // cout << "FIRST IMPLEMENTATION";
    // cout << endl;
    // sieve_remainder(n);
    // cout << endl;

    cout << "SECOND IMPLEMENTATION";
    cout << endl;
    sieve_fast_marking(n);
    cout << endl;

    cout << "THIRD IMPLEMENTATION";
    cout << endl;
    sieve_cache_friendly(n, n);
    cout << endl;
}