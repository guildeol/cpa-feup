#include <iostream>
#include <vector>

#include <cstdio>
#include <cmath>
#include <cstring>

#include "stdint.h"

#include "papi_helper.h"

using namespace std;

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

int main (int argc, char *argv[])
{
    PapiHelper ph;
    uint64_t n;
    
    cout << "Power of 10: ";
    cin >> n;
 
    n = pow(10,n);
    
    sieve_remainder(n);
    sieve_fast_marking(n);
}