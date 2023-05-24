#include <iostream>
#include <mpi.h>

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
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Assume A and B are square matrices of the same size
    int matrix_size = 8; // Size of each matrix

    int submatrix_size = (matrix_size *matrix_size)/ size;
    int* matrix_A = nullptr;
    int* matrix_B = nullptr;
    int* submatrix_A = allocateMatrix(submatrix_size, submatrix_size);
    int* submatrix_B = allocateMatrix(submatrix_size, submatrix_size);
    int* submatrix_C = allocateMatrix(submatrix_size, submatrix_size);

    if (rank == 0) {
        // Initialize the original matrices A and B
        matrix_A = allocateMatrix(matrix_size, matrix_size);
        matrix_B = allocateMatrix(matrix_size, matrix_size);

        for (int i = 0; i < matrix_size; ++i) {
            for (int j = 0; j < matrix_size; ++j) {
                matrix_A[i * matrix_size + j] = 1;
                matrix_B[i * matrix_size + j] = (i + 1);
            }
        }

        std::cout << "Matrix A:" << std::endl;
        printMatrix(matrix_A, matrix_size, matrix_size);
        std::cout << std::endl;
        std::cout << "Matrix B:" << std::endl;
        printMatrix(matrix_B, matrix_size, matrix_size);
        std::cout << std::endl;
    }

    // Scatter the submatrices A and B to all processes
    MPI_Scatter(matrix_A, submatrix_size * matrix_size, MPI_INT, submatrix_A, submatrix_size * matrix_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(matrix_B, matrix_size * submatrix_size, MPI_INT, submatrix_B, matrix_size * submatrix_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Perform SUMMA algorithm for matrix multiplication
    int* temp_matrix = allocateMatrix(submatrix_size, submatrix_size);

    for (int k = 0; k < size; ++k) {
        // Broadcast submatrix_A to all processes in the same row
        MPI_Bcast(submatrix_A, submatrix_size * submatrix_size, MPI_INT, k, MPI_COMM_WORLD);

        // Perform local matrix multiplication
        for (int i = 0; i < submatrix_size; ++i) {
            for (int j = 0; j < submatrix_size; ++j) {
                temp_matrix[i * submatrix_size + j] = 0;
                for (int p = 0; p < matrix_size; ++p) {
                    temp_matrix[i * submatrix_size + j] += submatrix_A[i * matrix_size + p] * submatrix_B[p * submatrix_size + j];
                }
            }
        }

        // Accumulate the results of local matrix multiplication into submatrix_C
        for (int i = 0; i < submatrix_size; ++i) {
            for (int j = 0; j < submatrix_size; ++j) {
                submatrix_C[i * submatrix_size + j] += temp_matrix[i * submatrix_size + j];
            }
        }

        // Rotate submatrix_B in each row
        int source = (rank + 1) % size;
        int destination = (rank + size - 1) % size;
        MPI_Sendrecv_replace(submatrix_B, submatrix_size * submatrix_size, MPI_INT, destination, 0, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Gather all submatrices C to the master process
    int* recv_counts = new int[size];
    int* displacements = new int[size];
    for (int i = 0; i < size; ++i) {
        recv_counts[i] = submatrix_size * submatrix_size;
        displacements[i] = i * submatrix_size * submatrix_size;
    }

    int* matrix_C = nullptr;

    if (rank == 0) {
        matrix_C = allocateMatrix(matrix_size, matrix_size);
    }

    MPI_Gatherv(submatrix_C, submatrix_size * submatrix_size, MPI_INT, matrix_C, recv_counts, displacements, MPI_INT, 0, MPI_COMM_WORLD);

    // Print the final matrix C on the master process
    if (rank == 0) {
        std::cout << "Matrix C:" << std::endl;
        printMatrix(matrix_C, matrix_size, matrix_size);
        deallocateMatrix(matrix_C);
    }

    // Deallocate memory
    deallocateMatrix(submatrix_A);
    deallocateMatrix(submatrix_B);
    deallocateMatrix(submatrix_C);
    deallocateMatrix(temp_matrix);
    if (rank == 0) {
        deallocateMatrix(matrix_A);
        deallocateMatrix(matrix_B);
    }

    delete[] recv_counts;
    delete[] displacements;

    MPI_Finalize();

    return 0;
}
