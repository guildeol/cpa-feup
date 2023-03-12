// #include <omp.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <papi.h>

using namespace std;

#define SYSTEMTIME clock_t

int EventSet = PAPI_NULL;

void createCounters()
{
  int ret;

  ret = PAPI_library_init(PAPI_VER_CURRENT);
  if (ret != PAPI_VER_CURRENT)
    std::cout << "FAIL" << endl;

  ret = PAPI_create_eventset(&EventSet);
  if (ret != PAPI_OK)
    cout << "ERRO: create eventset" << endl;

  ret = PAPI_add_event(EventSet, PAPI_L1_DCM);
  if (ret != PAPI_OK)
    cout << "ERRO: PAPI_L1_DCM" << endl;

  ret = PAPI_add_event(EventSet, PAPI_L2_DCM);
  if (ret != PAPI_OK)
    cout << "ERRO: PAPI_L2_DCM" << endl;
}

void startCounting()
{
  int ret;

  ret = PAPI_reset(EventSet);
  if (ret != PAPI_OK)
    std::cout << "FAIL reset" << endl;

  // Start counting
  ret = PAPI_start(EventSet);
  if (ret != PAPI_OK)
    cout << "ERRO: Start PAPI" << endl;
}

long long *stopCounting()
{
  int ret;
  static long long values[2];

  ret = PAPI_stop(EventSet, values);
  if (ret != PAPI_OK)
    cout << "ERRO: Stop PAPI" << endl;

  return values;
}

void report(int dimension, double time, long long *values)
{
  printf("Report Matrix %d x %d:\n", dimension, dimension);
  printf("\tTime: %3.3lf s\n", time);
  printf("\tL1 DCM: %lld \n", values[0]);
  printf("\tL2 DCM: %lld \n", values[1]);
}

void destroyCounters()
{
  int ret;

  ret = PAPI_remove_event(EventSet, PAPI_L1_DCM);
  if (ret != PAPI_OK)
    std::cout << "FAIL remove event" << endl;

  ret = PAPI_remove_event(EventSet, PAPI_L2_DCM);
  if (ret != PAPI_OK)
    std::cout << "FAIL remove event" << endl;

  ret = PAPI_destroy_eventset(&EventSet);
  if (ret != PAPI_OK)
    std::cout << "FAIL destroy" << endl;
}

void OnMult(int m_ar, int m_br)
{
  SYSTEMTIME Time1, Time2;

  char st[100];
  double temp;
  int i, j, k;
  long long *values;

  double *pha, *phb, *phc;

  pha = (double *)malloc((m_ar * m_ar) * sizeof(double));
  phb = (double *)malloc((m_ar * m_ar) * sizeof(double));
  phc = (double *)malloc((m_ar * m_ar) * sizeof(double));

  for (i = 0; i < m_ar; i++)
    for (j = 0; j < m_ar; j++)
      pha[i * m_ar + j] = (double)1.0;

  for (i = 0; i < m_br; i++)
    for (j = 0; j < m_br; j++)
      phb[i * m_br + j] = (double)(i + 1);

  Time1 = clock();
  startCounting();

  for (i = 0; i < m_ar; i++)
  {
    for (j = 0; j < m_br; j++)
    {
      temp = 0;
      for (k = 0; k < m_ar; k++)
      {
        temp += pha[i * m_ar + k] * phb[k * m_br + j];
      }

      phc[i * m_ar + j] = temp;
    }
  }

  values = stopCounting();
  Time2 = clock();

  report(m_ar, (double)(Time2 - Time1)/CLOCKS_PER_SEC, values);

  cout << "\tphc[0] = "<< phc[0] << endl;

  free(pha);
  free(phb);
  free(phc);
}

void OnMultLine(int m_ar, int m_br)
{
  SYSTEMTIME Time1, Time2;

  char st[100];
  double temp;
  int i, j, k;
  long long *values;

  double *pha, *phb, *phc;

  pha = (double *)calloc((m_ar * m_ar), sizeof(double));
  phb = (double *)calloc((m_ar * m_ar), sizeof(double));
  phc = (double *)calloc((m_ar * m_ar), sizeof(double));

  for (i = 0; i < m_ar; i++)
    for (j = 0; j < m_ar; j++)
      pha[i * m_ar + j] = (double)1.0;

  for (i = 0; i < m_br; i++)
    for (j = 0; j < m_br; j++)
      phb[i * m_br + j] = (double)(i + 1);

  Time1 = clock();
  startCounting();

  for (i = 0; i < m_ar; i++)
  {
    for (k = 0; k < m_ar; k++)
    {
      for (j = 0; j < m_br; j++)
        phc[i * m_ar + j] += pha[i * m_ar + k] * phb[k * m_br + j];
    }
  }

  values = stopCounting();
  Time2 = clock();

  report(m_ar, (double)(Time2 - Time1)/CLOCKS_PER_SEC, values);

  cout << "\tphc[0] = "<< phc[0] << endl;

  free(pha);
  free(phb);
  free(phc);
}


int main(int argc, char *argv[])
{
  int min_dimension, max_dimension, step;
  int op = 1;
  void (*multAlgorithm)(int, int);

  createCounters();

  do
  {
    cout << endl
         << "1. Multiplication" << endl;
    cout << "2. Line Multiplication" << endl;
    cout << "Selection?: ";

    cin >> op;
    if (op == 0)
      break;

    printf("Please enter: Initial Matrix Dimension, Final Matrix Dimension, and Step Size: ");
    cin >> min_dimension >> max_dimension >> step;

    switch (op)
    {
      case 1:
        multAlgorithm = OnMult;
      break;

      case 2:
        multAlgorithm = OnMultLine;
      break;
    }

    for (auto dimension = min_dimension; dimension <= max_dimension; dimension += step)
    {
      multAlgorithm(dimension, dimension);
    }

  } while (op != 0);

  destroyCounters();

  return 0;
}