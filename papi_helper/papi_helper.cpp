#include <iostream>
#include <cstdio>
#include <ctime>
#include <cstring>

#include <stdint.h>
#include <papi.h>

#include "papi_helper.h"

using namespace std;

PapiHelper::PapiHelper()
{
  int ret;

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

}

PapiHelper::~PapiHelper()
{
  int ret;

  ret = PAPI_remove_event(this->EventSet, PAPI_L1_DCM);
  if (ret != PAPI_OK)
    std::cout << "FAIL remove event " << ret << endl;


  ret = PAPI_remove_event(this->EventSet, PAPI_L2_DCM);
  if (ret != PAPI_OK)
    std::cout << "FAIL remove event " << ret << endl;


  ret = PAPI_destroy_eventset(&this->EventSet);
  if (ret != PAPI_OK)
    std::cout << "FAIL destroy " << ret << endl;

}

void PapiHelper::startCounting()
{
  int ret;

  ret = PAPI_reset(this->EventSet);
  if (ret != PAPI_OK)
    std::cout << "FAIL reset " << ret << endl;


  this->countStarted = clock();

  // Start counting
  ret = PAPI_start(this->EventSet);
  if (ret != PAPI_OK)
    cout << "ERROR: Start PAPI " << ret << endl;

}

void PapiHelper::stopCounting()
{
  int ret;

  ret = PAPI_stop(this->EventSet, (long long*)this->counters);
  if (ret != PAPI_OK)
    cout << "ERROR: Stop PAPI " << ret << endl;


  this->countStopped = clock();
}

void PapiHelper::report()
{
  double time_diff = (this->countStopped - this->countStarted)/CLOCKS_PER_SEC;

  cout << "Performance Report:" << endl;
  cout << "\tTime: " <<  time_diff << " s"  << endl;
  cout << "\tL1 DCM: " << this->counters[L1_MISSES_IDX] << endl;
  cout << "\tL2 DCM: " << this->counters[L2_MISSES_IDX] << endl;
}
