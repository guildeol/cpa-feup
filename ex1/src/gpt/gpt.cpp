#include <iostream>
#include <cstdlib>
#include <ctime>
#include <papi.h>

void matrix_multiply(int* A, int* B, int* C, int n) {
    int event_set = PAPI_NULL;
    long long values[2];

    // Initialize PAPI and add the L1 and L2 cache miss events to the event set
    PAPI_library_init(PAPI_VER_CURRENT);
    PAPI_create_eventset(&event_set);
    PAPI_add_event(event_set, PAPI_L1_DCM);
    PAPI_add_event(event_set, PAPI_L2_DCM);

    // Start counting cache misses
    PAPI_start(event_set);

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                C[i*n + j] += A[i*n + k] * B[k*n + j];
            }
        }
    }

    // Stop counting cache misses and read the results
    PAPI_stop(event_set, values);

    // Print the results
    std::cout << "L1 cache misses: " << values[0] << std::endl;
    std::cout << "L2 cache misses: " << values[1] << std::endl;
}

int main() {
    int n;

    // Prompt the user to input the value of n
    std::cout << "Enter the value of n: ";
    std::cin >> n;

    // Allocate memory for matrices A, B, and C using a 1D array
    int* A = new int[n * n];
    int* B = new int[n * n];
    int* C = new int[n * n]();

    // Initialize matrices A and B
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i*n + j] = 1;
            B[i*n + j] = i+1;
        }
    }

    // Measure the execution time of the matrix_multiply function
    clock_t start_time = clock();
    matrix_multiply(A, B, C, n);
    clock_t end_time = clock();

    // Calculate the execution time in seconds
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    // Print the first element of the result matrix C
    std::cout << "C[0][0] = " << C[0] << std::endl;

    std::cout << "Execution time: " << total_time << " seconds" << std::endl;

    // Deallocate memory
    delete[] A;
    delete[] B;
    delete[] C;

    return 0;
}
