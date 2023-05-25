#include <iostream>
#include <mpi.h>
#include <cmath>

void printMatrix(int* matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << matrix[i * cols + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int* allocateMatrix(int rows, int cols) {
    return new int[rows * cols]();
}

void deallocateMatrix(int* matrix) {
    delete[] matrix;
}

int main(int argc, char** argv) {
    double global_start, global_end;
    double communication_start, communication_end;
    double communication_accumulated;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <matrix_size>" << std::endl;
        return 1;
    }

    // Assume A and B are square matrices of the same size
    // Size of each matrix
    int matrix_size = std::stoi(argv[1]);

    MPI_Init(&argc, &argv);

    global_start = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int size_squared = round(sqrt(size));
    int submatrix_size = matrix_size / size_squared;
    int* matrix_A = nullptr;
    int* matrix_B = nullptr;
    int* submatrix_A = allocateMatrix(submatrix_size, submatrix_size);
    int* submatrix_B = allocateMatrix(submatrix_size, submatrix_size);
    int* submatrix_C = allocateMatrix(submatrix_size, submatrix_size);
    int submatrix_size_squared = submatrix_size * submatrix_size;

    // Initialize the original matrices A and B
    matrix_A = allocateMatrix(matrix_size, matrix_size);
    matrix_B = allocateMatrix(matrix_size, matrix_size);

    if (rank == 0) {

        for (int i = 0; i < matrix_size; ++i) {
            for (int j = 0; j < matrix_size; ++j) {
                matrix_A[i * matrix_size + j] = 1;
                matrix_B[i * matrix_size + j] = (i + 1);
            }
        }

        // Print Matrix A and Matrix B at the beginning
        // std::cout << "Matrix A:" << std::endl;
        // printMatrix(matrix_A, matrix_size, matrix_size);
        // std::cout << std::endl;
        // std::cout << "Matrix B:" << std::endl;
        // printMatrix(matrix_B, matrix_size, matrix_size);
        // std::cout << std::endl;
    }

    //Count communication time
    communication_start = MPI_Wtime();
    // Send the entire matrix to all processes
    MPI_Bcast(matrix_A, matrix_size * matrix_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(matrix_B, matrix_size * matrix_size, MPI_INT, 0, MPI_COMM_WORLD);
    communication_end = MPI_Wtime();
    communication_accumulated += communication_end - communication_start;

    // Extract the n square of matrix A and matrix B on each process
    for (int i = 0; i < submatrix_size; i++) {
        for (int j = 0; j < submatrix_size; j++) {
            submatrix_A[i * submatrix_size + j] = matrix_A[(i + submatrix_size * (rank / size_squared)) * matrix_size + submatrix_size * (rank % size_squared) + j];
            submatrix_B[i * submatrix_size + j] = matrix_B[(i + submatrix_size * (rank / size_squared)) * matrix_size + submatrix_size * (rank % size_squared) + j];
        }
    }

    // Perform SUMMA algorithm for matrix multiplication

    // Temporary matrix 
    int* temp_matrix = allocateMatrix(submatrix_size, submatrix_size);

    for (int k = 0; k < size; ++k) {
        // Count communication time
        communication_start = MPI_Wtime();
        // Broadcast submatrix_A to all processes
        MPI_Bcast(submatrix_A, submatrix_size_squared, MPI_INT, k, MPI_COMM_WORLD);
        communication_end = MPI_Wtime();
        communication_accumulated += communication_end - communication_start;

        // Perform local matrix multiplication
        for (int i = 0; i < submatrix_size; ++i) {
            for (int j = 0; j < submatrix_size; ++j) {
                temp_matrix[i * submatrix_size + j] = 0;
                for (int p = 0; p < submatrix_size; ++p) {
                    temp_matrix[i * submatrix_size + j] += submatrix_A[i * submatrix_size + p] * submatrix_B[p * submatrix_size + j];
                }
            }
        }

        // Accumulate the results of local matrix multiplication into submatrix_C
        for (int i = 0; i < submatrix_size; ++i) {
            for (int j = 0; j < submatrix_size; ++j) {
                submatrix_C[i * submatrix_size + j] += temp_matrix[i * submatrix_size + j]/sqrt(size); 
            }
        }

        // Rotate submatrix_B in each row
        int source = (rank + 1) % size;
        int destination = (rank + size - 1) % size;
        communication_start = MPI_Wtime();
        MPI_Sendrecv_replace(submatrix_B, submatrix_size_squared, MPI_INT, destination, 0, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        communication_end = MPI_Wtime();
        communication_accumulated += communication_end - communication_start;
    }

    int* recv_counts = new int[size];
    int* displacements = new int[size];
    for (int i = 0; i < size; ++i) {
        recv_counts[i] = submatrix_size_squared;
        displacements[i] = submatrix_size_squared * i;
    }

    int* matrix_C = nullptr;

    if (rank == 0) {
        matrix_C = allocateMatrix(matrix_size, matrix_size);
    }   

    // Count communication time
    communication_start = MPI_Wtime();
    // Gather all submatrices C to the master process
    MPI_Gatherv(submatrix_C, submatrix_size_squared, MPI_INT, matrix_C, recv_counts, displacements, MPI_INT, 0, MPI_COMM_WORLD);
    communication_end = MPI_Wtime();
    communication_accumulated += communication_end - communication_start;

    // Print the final matrix C on the master process
    if (rank == 0) {
        // std::cout << "Matrix C:" << std::endl;
        // printMatrix(matrix_C, matrix_size, matrix_size);

        std::cout << "Matrix C first cell: ";
        std::cout << matrix_C[0] << std::endl;

        deallocateMatrix(matrix_C);
    }

    // Deallocate memory
    deallocateMatrix(submatrix_A);
    deallocateMatrix(submatrix_B);
    deallocateMatrix(submatrix_C);
    deallocateMatrix(temp_matrix);
    if (rank == 0) {
        global_end = MPI_Wtime();
        printf("Total Ellapsed time: %.10lf seconds\n", global_end - global_start);
        printf("Communication time: %lf\n", communication_accumulated);

        deallocateMatrix(matrix_A);
        deallocateMatrix(matrix_B);
    }

    delete[] recv_counts;
    delete[] displacements;

    MPI_Finalize();

    return 0;
}
