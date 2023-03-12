#ifndef PAPI_HELPER_H
#define PAPI_HELPER_H

constexpr int CACHE_MISS_COUNTERS = 2;
constexpr int L1_MISSES_IDX = 2;
constexpr int L2_MISSES_IDX = 1;

class PapiHelper
{
    private:
        int EventSet;
        long long counters[CACHE_MISS_COUNTERS];
        clock_t countStarted, countStopped;

    public:
        PapiHelper();
        ~PapiHelper();
        void startCounting();
        void stopCounting();
        void destroyCounters();
        void report();
};

#endif