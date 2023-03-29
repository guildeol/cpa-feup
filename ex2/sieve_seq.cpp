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

#define BLOCK_LOW(id, p, n)         ((id) * (n)/(p))
#define BLOCK_HIGH(id, p, n)        (BLOCK_LOW((id) + 1, (p), (n)) - 1)
#define BLOCK_SIZE(id, p, n)        (BLOCK_HIGH((id), (p), (n)) - BLOCK_LOW((id), (p), (n)))
#define BLOCK_OWNER(index, p, n)    (((p) * ((index) + 1) - 1) / (n))

#define IS_ODD(_num_) ((_num_) & 0x01)

bool print = false;
bool skip_first_impl = false;

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
void sieve_cache_friendly(uint64_t n, uint32_t num_blocks) 
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
        for (uint64_t id = 0; id < num_blocks; id++)
        {
            uint64_t first;
            uint64_t remainder;

            uint64_t low_value = BLOCK_LOW(id, num_blocks, n - 1);
            uint64_t high_value = 2 + BLOCK_HIGH(id, num_blocks, n -1);
            uint64_t block_size = BLOCK_SIZE(id, num_blocks, n - 1);

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

            for (long long j = (low_value + first); j < high_value; j += 2*k)
            {   
                if (!IS_ODD(j))
                    j += k;

                CLEAR_BIT(primes, j >> 1);
            }
        }
        
        do
        {
            k += 2;
        } while (k*k <= n && !GET_BIT(primes, k >> 1));
        
    } while (k*k <= n);

    ph.stopCounting();

    if (print)
    {
        cout << "2 ";
        for (int i = 3; i < n; i += 2)
            if (GET_BIT(primes, i >> 1))
                cout << i << " ";

        cout << endl;
    }

    ph.report();

    delete[] primes;
}

int main (int argc, char *argv[])
{
    PapiHelper ph;
    uint64_t n;
    uint32_t num_blocks;

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <pow_2> <num_block> [print (0/1)] [skip impl. 1 (0/1)]" << endl;
        return -1;
    }

    n = atoi(argv[1]);
    n = pow(2, n);

    num_blocks = atoi(argv[2]);

    if (argc > 3)
        print = (bool)atoi(argv[3]);

    if (argc > 4)
        skip_first_impl = (bool)atoi(argv[4]);
    
    if(!skip_first_impl)
    {
        cout << "FIRST IMPLEMENTATION";
        cout << endl;
        sieve_remainder(n);
        cout << endl;
    }

    cout << "SECOND IMPLEMENTATION";
    cout << endl;
    sieve_fast_marking(n);
    cout << endl;

    cout << "THIRD IMPLEMENTATION";
    cout << endl;
    sieve_cache_friendly(n, num_blocks);
    cout << endl;
}