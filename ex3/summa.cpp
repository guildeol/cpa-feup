#include <iostream>
#include <mpi.h>

void printMatrix(int** matrix, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int** allocateMatrix(int rows, int cols) {
    int** matrix = new int*[rows];
    for (int i = 0; i < rows; ++i) {
        matrix[i] = new int[cols];
    }
    return matrix;
}

void deallocateMatrix(int** matrix, int rows) {
    for (int i = 0; i < rows; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Assume A and B are square matrices of the same size
    int matrix_size = 4; // Size of each matrix

    int submatrix_size = matrix_size / size;
    int** submatrix_A = allocateMatrix(submatrix_size, submatrix_size);
    int** submatrix_B = allocateMatrix(submatrix_size, submatrix_size);
    int** submatrix_C = allocateMatrix(submatrix_size, submatrix_size);

    // Generate data for submatrices A and B
    for (int i = 0; i < submatrix_size; ++i) {
        for (int j = 0; j < submatrix_size; ++j) {
            submatrix_A[i][j] = 1.0;
            submatrix_B[i][j] = (double)(i + 1);
        }
    }

    // Perform SUMMA algorithm for matrix multiplication
    int** temp_matrix = allocateMatrix(submatrix_size, submatrix_size);
    
    for (int k = 0; k < size; ++k) {
        // Broadcast submatrix_A to all processes in the same column
        MPI_Bcast(&(submatrix_A[0][0]), submatrix_size * submatrix_size, MPI_INT, k, MPI_COMM_WORLD);

        // Perform local matrix multiplication
        for (int i = 0; i < submatrix_size; ++i) {
            for (int j = 0; j < submatrix_size; ++j) {
                temp_matrix[i][j] = 0;
                for (int p = 0; p < submatrix_size; ++p) {
                    temp_matrix[i][j] += submatrix_A[i][p] * submatrix_B[p][j];
                }
            }
        }

        // Accumulate the results of local matrix multiplication into submatrix_C
        for (int i = 0; i < submatrix_size; ++i) {
            for (int j = 0; j < submatrix_size; ++j) {
                submatrix_C[i][j] += temp_matrix[i][j];
            }
        }

        // Rotate submatrix_B in each row
        MPI_Sendrecv_replace(&(submatrix_B[0][0]), submatrix_size * submatrix_size, MPI_INT,
                             (rank + 1) % size, 0,
                             (rank + size - 1) % size, 0,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Gather all submatrices C to the master process
    int** matrix_C = nullptr;
    if (rank == 0) {
        matrix_C = allocateMatrix(matrix_size, matrix_size);
    }
    MPI_Gather(&(submatrix_C[0][0]), submatrix_size * submatrix_size, MPI_INT,
               &(matrix_C[0][0]), submatrix_size * submatrix_size, MPI_INT,
               0, MPI_COMM_WORLD);

    // Print the final matrix C on the master process
    if (rank == 0) {
        std::cout << "Matrix C:" << std::endl;
        printMatrix(matrix_C, matrix_size, matrix_size);
        deallocateMatrix(matrix_C, matrix_size);
    }

    // Deallocate memory
    deallocateMatrix(submatrix_A, submatrix_size);
    deallocateMatrix(submatrix_B, submatrix_size);
    deallocateMatrix(submatrix_C, submatrix_size);
    deallocateMatrix(temp_matrix, submatrix_size);

    MPI_Finalize();

    return 0;
}
