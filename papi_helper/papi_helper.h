#ifndef PAPI_HELPER_H
#define PAPI_HELPER_H

#include <chrono>

constexpr int CACHE_MISS_COUNTERS = 2;
constexpr int L1_MISSES_IDX = 0;
constexpr int L2_MISSES_IDX = 1;

using namespace std;

class PapiHelper
{
    private:
        int EventSet;
        uint64_t counters[CACHE_MISS_COUNTERS];
        chrono::time_point<chrono::high_resolution_clock> countStarted, countStopped;

    public:
        PapiHelper();
        ~PapiHelper();
        void startCounting();
        void stopCounting();
        void destroyCounters();
        void report();
};

#endif