#include <iostream>
#include <chrono>
#include <cstdio>
#include <cstring>

#include <stdint.h>
#include <papi.h>

#include "papi_helper.h"

using namespace std;

PapiHelper::PapiHelper()
{
  int ret;

#ifndef DISABLE_PAPI_COUNTERS
  this->EventSet = PAPI_NULL;

  ret = PAPI_library_init(PAPI_VER_CURRENT);
  if (ret != PAPI_VER_CURRENT)
    std::cout << "FAIL " << ret << endl;

  ret = PAPI_create_eventset(&this->EventSet);
  if (ret != PAPI_OK)
    cout << "ERROR: create this->EventSet: " << ret << endl;

  ret = PAPI_add_event(this->EventSet, PAPI_L1_DCM);
  if (ret != PAPI_OK)
    cout << "ERROR: PAPI_L1_DCM " << ret << endl;

  ret = PAPI_add_event(this->EventSet, PAPI_L2_DCM);
  if (ret != PAPI_OK)
    cout << "ERROR: PAPI_L2_DCM " << ret << endl;
#endif
}

PapiHelper::~PapiHelper()
{
  int ret;
  
#ifndef DISABLE_PAPI_COUNTERS
  ret = PAPI_remove_event(this->EventSet, PAPI_L1_DCM);
  if (ret != PAPI_OK)
    std::cout << "FAIL remove event " << ret << endl;


  ret = PAPI_remove_event(this->EventSet, PAPI_L2_DCM);
  if (ret != PAPI_OK)
    std::cout << "FAIL remove event " << ret << endl;


  ret = PAPI_destroy_eventset(&this->EventSet);
  if (ret != PAPI_OK)
    std::cout << "FAIL destroy " << ret << endl;
#endif
}

void PapiHelper::startCounting()
{
  int ret;

  this->countStarted = chrono::high_resolution_clock::now();

#ifndef DISABLE_PAPI_COUNTERS
  ret = PAPI_reset(this->EventSet);
  if (ret != PAPI_OK)
    std::cout << "FAIL reset " << ret << endl;

  // Start counting
  ret = PAPI_start(this->EventSet);
  if (ret != PAPI_OK)
    cout << "ERROR: Start PAPI " << ret << endl;
#endif
}

void PapiHelper::stopCounting()
{
  int ret;

#ifndef DISABLE_PAPI_COUNTERS
  ret = PAPI_stop(this->EventSet, (long long*)this->counters);
  if (ret != PAPI_OK)
    cout << "ERROR: Stop PAPI " << ret << endl;
#endif

  this->countStopped = chrono::high_resolution_clock::now();
}

void PapiHelper::report()
{
  std::chrono::duration<double, std::milli> ellapsed_ms = (this->countStopped - this->countStarted);

  cout << "Performance Report:" << endl;
  cout << "\tTime: " <<  ellapsed_ms.count() << " ms"  << endl;
#ifndef DISABLE_PAPI_COUNTERS
  cout << "\tL1 DCM: " << this->counters[L1_MISSES_IDX] << endl;
  cout << "\tL2 DCM: " << this->counters[L2_MISSES_IDX] << endl;
#endif
}
