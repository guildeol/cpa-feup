#include <iostream>
#include <stdint.h>
#include <mpi.h>

using namespace std;

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

    printf("Size - ");
    printf("%d", size);
    std::cout << std::endl;
    printf("Rank - ");
    printf("%d", rank);
    std::cout << std::endl;


    // Assume A and B are square matrices of the same size
    int matrix_size = 8; // Size of each matrix

    int submatrix_size = matrix_size / size;
    int* submatrix_A = allocateMatrix(submatrix_size, matrix_size);
    int* submatrix_B = allocateMatrix(matrix_size, submatrix_size);
    int* submatrix_C = allocateMatrix(submatrix_size, submatrix_size);

    if(rank == 0){
        // Generate data for submatrices A and B
        // for (int i = 0; i < submatrix_size; ++i) {
        //     for (int j = 0; j < submatrix_size; ++j) {
        //         submatrix_A[i * submatrix_size + j] = 1.0;
        //         submatrix_B[i * submatrix_size + j] = (double)(i + 1);
        //     }
        // }

        // Generate data for submatrices A and B
        for (int i = 0; i < submatrix_size; ++i) {
            for (int j = 0; j < submatrix_size; ++j) {
                int aux = i * submatrix_size + j;
                printf("Position - %d", aux);
                std::cout << std::endl;
                submatrix_A[i * submatrix_size + j] = 1;
            }
        }

        for (int i = 0; i < matrix_size; ++i) {
            for (int j = 0; j < submatrix_size; ++j) {
                submatrix_B[i * submatrix_size + j] = (i + 1);
            }
        }
        std::cout << "Matrix A:" << std::endl;
        printMatrix(submatrix_A, matrix_size, matrix_size);
        std::cout << std::endl;
        std::cout << "Matrix B:" << std::endl;
        printMatrix(submatrix_B, matrix_size, matrix_size);
        std::cout << std::endl;
    }

    // Perform SUMMA algorithm for matrix multiplication
    int* temp_matrix = allocateMatrix(submatrix_size, submatrix_size);

    for (int k = 0; k < size; ++k) {
        // Broadcast submatrix_A to all processes in the same row
        // MPI_Bcast(submatrix_A, submatrix_size * submatrix_size, MPI_INT, k / submatrix_size, MPI_COMM_WORLD);
        MPI_Bcast(submatrix_A, submatrix_size * submatrix_size, MPI_INT, (k % size) / submatrix_size, MPI_COMM_WORLD);


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
                submatrix_C[i * submatrix_size + j] += temp_matrix[i * submatrix_size + j];
            }
        }

        // Rotate submatrix_B in each row
        int source = (rank + 1) % size;
        int destination = (rank + size - 1) % size;
        MPI_Sendrecv_replace(submatrix_B, submatrix_size * submatrix_size, MPI_INT, destination, 0, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }

    // Gather all submatrices C to the master process
    int* matrix_C = nullptr;
    if (rank == 0) {
        matrix_C = allocateMatrix(matrix_size, matrix_size);
    }
    MPI_Gather(submatrix_C, submatrix_size * submatrix_size, MPI_INT, matrix_C, submatrix_size * submatrix_size, MPI_INT, 0, MPI_COMM_WORLD);

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

    MPI_Finalize();

    return 0;
}























// #include <iostream>
// #include <mpi.h>

// void printMatrix(int** matrix, int rows, int cols) {
//     for (int i = 0; i < rows; ++i) {
//         for (int j = 0; j < cols; ++j) {
//             std::cout << matrix[i][j] << " ";
//         }
//         std::cout << std::endl;
//     }
//     std::cout << std::endl;
// }

// int** allocateMatrix(int rows, int cols) {
//     int** matrix = new int*[rows];
//     for (int i = 0; i < rows; ++i) {
//         matrix[i] = new int[cols];
//     }
//     return matrix;
// }

// void deallocateMatrix(int** matrix, int rows) {
//     for (int i = 0; i < rows; ++i) {
//         delete[] matrix[i];
//     }
//     delete[] matrix;
// }

// int main(int argc, char** argv) {
//     MPI_Init(&argc, &argv);

//     int rank, size;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);

//     // Assume A and B are square matrices of the same size
//     int matrix_size = 4; // Size of each matrix

//     int submatrix_size = matrix_size / size;
//     int** submatrix_A = allocateMatrix(submatrix_size, submatrix_size);
//     int** submatrix_B = allocateMatrix(submatrix_size, submatrix_size);
//     int** submatrix_C = allocateMatrix(submatrix_size, submatrix_size);

//     // Generate data for submatrices A and B
//     for (int i = 0; i < submatrix_size; ++i) {
//         for (int j = 0; j < submatrix_size; ++j) {
//             submatrix_A[i][j] = 1.0;
//             submatrix_B[i][j] = (double)(i + 1);
//         }
//     }

//     // Perform SUMMA algorithm for matrix multiplication
//     int** temp_matrix = allocateMatrix(submatrix_size, submatrix_size);
    
//     for (int k = 0; k < size; ++k) {
//         // Broadcast submatrix_A to all processes in the same column
//         MPI_Bcast(&(submatrix_A[0][0]), submatrix_size * submatrix_size, MPI_INT, k, MPI_COMM_WORLD);

//         // Perform local matrix multiplication
//         for (int i = 0; i < submatrix_size; ++i) {
//             for (int j = 0; j < submatrix_size; ++j) {
//                 temp_matrix[i][j] = 0;
//                 for (int p = 0; p < submatrix_size; ++p) {
//                     temp_matrix[i][j] += submatrix_A[i][p] * submatrix_B[p][j];
//                 }
//             }
//         }

//         // Accumulate the results of local matrix multiplication into submatrix_C
//         for (int i = 0; i < submatrix_size; ++i) {
//             for (int j = 0; j < submatrix_size; ++j) {
//                 submatrix_C[i][j] += temp_matrix[i][j];
//             }
//         }

//         // Rotate submatrix_B in each row
//         MPI_Sendrecv_replace(&(submatrix_B[0][0]), submatrix_size * submatrix_size, MPI_INT,
//                              (rank + 1) % size, 0,
//                              (rank + size - 1) % size, 0,
//                              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//     }

//     // Gather all submatrices C to the master process
//     int** matrix_C = nullptr;
//     if (rank == 0) {
//         matrix_C = allocateMatrix(matrix_size, matrix_size);
//     }
//     MPI_Gather(&(submatrix_C[0][0]), submatrix_size * submatrix_size, MPI_INT,
//                &(matrix_C[0][0]), submatrix_size * submatrix_size, MPI_INT,
//                0, MPI_COMM_WORLD);

//     // Print the final matrix C on the master process
//     if (rank == 0) {
//         std::cout << "Matrix C:" << std::endl;
//         printMatrix(matrix_C, matrix_size, matrix_size);
//         deallocateMatrix(matrix_C, matrix_size);
//     }

//     // Deallocate memory
//     deallocateMatrix(submatrix_A, submatrix_size);
//     deallocateMatrix(submatrix_B, submatrix_size);
//     deallocateMatrix(submatrix_C, submatrix_size);
//     deallocateMatrix(temp_matrix, submatrix_size);

//     MPI_Finalize();

//     return 0;
// }
