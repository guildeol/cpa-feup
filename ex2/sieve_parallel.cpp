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

#define IS_ODD(_num_) ((_num_) & 0x01)

uint8_t *primes;

void sieve_parallel(uint64_t n, uint32_t threads) 
{
    PapiHelper ph;
    uint64_t k = 3;

    ph.startCounting();

    do
    {
        #pragma omp parallel num_threads(threads) firstprivate(k)
        {
            uint64_t first;
            uint64_t remainder;
            uint64_t id = omp_get_thread_num();

            uint64_t low_value = BLOCK_LOW(id, threads, n - 1);
            uint64_t high_value = 2 + BLOCK_HIGH(id, threads, n -1);
            uint64_t block_size = BLOCK_SIZE(id, threads, n - 1);

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
                // if (id == 1)
                // {
                //     if (j == (low_value + first))
                //     {
                //         cout << "Marking multiples of " << k << endl;
                //         cout << "This first multiple of k in this block is " << (low_value + first) << endl;
                //     }
                
                //     cout << "Thread " << id << " clearing bit " << j <<" [" << k << "]" << endl;
                // }

                // Sometimes we might end up in an even multiple of k, so if that is the case, skip it.
                if (!IS_ODD(j))
                    j += k;

                CLEAR_BIT(primes, j >> 1);
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
    uint32_t size;
    bool print = false;

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <pow_10> <threads> [print (0/1)]" << endl;
        return -1;
    }
    
    n = atoi(argv[1]);
    n = pow(10, n);

    threads = atoi(argv[2]);

    if (argc > 3)
        print = (bool)atoi(argv[3]);

    cout << "PARALLEL IMPLEMENTATION";
    cout << endl;

    size = (n / 8) + 1;
    primes = new uint8_t[size];

    // Start by marking all values as primes
    memset(primes, 0xFF, size * sizeof(uint8_t));

    sieve_parallel(n, threads);

    if (print)
    {
        cout << "2 ";
        for (int i = 3; i < n; i += 2)
            if (GET_BIT(primes, i >> 1))
                cout << i << " ";

        cout << endl;
    }

    delete[] primes;
}