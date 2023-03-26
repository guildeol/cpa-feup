#include <iostream>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstring>

#include "stdint.h"

#include "papi_helper.h"

using namespace std;

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

    cout << "2 ";
    for (int i = 3; i < n; i += 2)
        if (primes[i])
            cout << i << " ";
    
    cout << endl;

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

    cout << "2 ";
    for (int i = 3; i < n; i += 2)
        if (primes[i >> 1])
            cout << i << " ";

    cout << endl;

    ph.report();

    delete[] primes;
}

// Third implementation - we use the second implementation and reorganize computation so that the cache misses are reduced
// by searching several seed numbers in the same data block
void sieve_blocks(uint64_t n, uint64_t block_size) {
    PapiHelper ph;
    uint64_t k = 3;
    bool *primes = new bool[n];

    // Start by marking all values as primes
    memset(primes, true, n * sizeof(bool));

    // Mark 0 and 1 as not primes
    primes[0] = primes[1] = false;

    ph.startCounting();

    // Loop for blocks of numbers
    for (uint64_t block_start = 2; block_start <= n; block_start += block_size) {
        // Calculate end of the block
        uint64_t block_end = min(block_start + block_size - 1, n);

        // Calculate maximum value of k for this block, which is equal to the square root of the end of the block
        uint64_t max_k = sqrt(block_end);

        // Loop that handles marking multiples of primes
        for (uint64_t i = 2; i <= max_k; i++){
            if (primes[i]) {
                uint64_t start = (block_start + i - 1) / i;
                uint64_t j = max(i, start) * i - block_start;

                for (; j <= block_end - block_start; j += i){
                    primes[j + block_start] = false;
                }
            }
        }

        // Print primes in current block
        for (uint64_t i = block_start; i <= block_end; i++){
            if (primes[i]){
                cout << i << " ";
            }
        }
    }

    ph.stopCounting();

    cout << endl;

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
    
    cout << "FIRST IMPLEMENTATION";
    cout << endl;
    sieve_remainder(n);
    cout << endl;

    cout << "SECOND IMPLEMENTATION";
    cout << endl;
    sieve_fast_marking(n);
    cout << endl;

    cout << "THIRD IMPLEMENTATION";
    cout << endl;
    sieve_blocks(n, 10);
    cout << endl;
}