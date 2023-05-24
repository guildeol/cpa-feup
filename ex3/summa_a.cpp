#include <mpi.h>
#include <iostream>
#include <cstring>
#include <memory>
#include <functional>
#include <random>
#include <cmath>

static void MPI_Comm_Deleter(MPI_Comm *comm) {
    MPI_Comm_free(comm);
}

static void MPI_Datatype_Deleter(MPI_Datatype *p) {
    MPI_Type_free(p);
}

// template<typename T>
// using deleted_unique_ptr = std::unique_ptr<T, std::function<void(T *)>>;

template<typename T>
struct deleted_unique_ptr : std::unique_ptr<T, std::function<void(T *)>> {};


void summa(MPI_Comm comm_grid, const int row_A, const int col_A, const int row_B, const int col_B, const int nb,
           double A_local[], double B_local[], double C_local[]) {
    if (col_A % nb != 0 || row_B % nb != 0) {
        std::cerr << "k must be multiple of nb" << std::endl;
        MPI_Abort(comm_grid, 1);
        return;
    }
    // int n_block_per_process = k / nb;
    int nb_col_pp = col_A / nb;
    int nb_row_pp = row_B / nb;
    int rank, size, rank_row, rank_col, size_row, size_col;
    MPI_Comm_rank(comm_grid, &rank);
    MPI_Comm_size(comm_grid, &size);
    int dims[2];
    int periods[2];
    int my_coords[2];
    MPI_Cart_get(comm_grid, 2, dims, periods, my_coords);
    if (nb_col_pp * dims[1] != nb_row_pp * dims[0]) {
        std::cerr << "Not multiplicable matrices" << std::endl;
        MPI_Abort(comm_grid, 1);
        return;
    }
    int nb_global = nb_col_pp * dims[1];
    int remain_dims[2] = {0, 1};
    MPI_Comm row_comm;
    MPI_Comm col_comm;
    MPI_Cart_sub(comm_grid, remain_dims, &row_comm);
    deleted_unique_ptr<MPI_Comm> _row_comm(&row_comm, MPI_Comm_Deleter);
    MPI_Comm_rank(row_comm, &rank_row);
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(comm_grid, remain_dims, &col_comm);
    deleted_unique_ptr<MPI_Comm> _col_comm(&col_comm, MPI_Comm_Deleter);
    MPI_Comm_rank(col_comm, &rank_col);
    MPI_Comm_size(row_comm, &size_row);
    MPI_Comm_size(col_comm, &size_col);
    if (nb_col_pp * size_row != nb_row_pp * size_col) {
        std::cerr << "Not multiplicable matrices" << std::endl;
        MPI_Abort(comm_grid, 1);
        return;
    }
    memset(C_local, 0, row_A * col_B * sizeof(double));
    double *buffer_A = new double[nb * row_A];
    std::unique_ptr<double[]> _buffer_A(buffer_A);
    double *buffer_B = new double[nb * col_B];
    std::unique_ptr<double[]> _buffer_B(buffer_B);
    for (int k = 0; k < nb_global; ++k) {
        double *buff_ptr_A, *buff_ptr_B;
        int owner_A = k / nb_col_pp, owner_B = k / nb_row_pp;
        if (owner_A == rank_row) {
            // We still need to collect part of A to contiguous memory
            buff_ptr_A = buffer_A;
            int idx = 0;
            for (int i = 0; i < row_A; ++i) {
                for (int j = 0; j < nb; ++j) {
                    buff_ptr_A[idx++] = A_local[(k % nb_col_pp) * nb + i * col_A + j];
                }
            }
        } else {
            buff_ptr_A = buffer_A;
        }
        MPI_Bcast(buff_ptr_A, nb * row_A, MPI_DOUBLE, owner_A, row_comm);
        // Memory of B's part is already contiguous
        if (owner_B == rank_col) {
            buff_ptr_B = &B_local[(k % nb_row_pp) * nb * col_B];
        } else {
            buff_ptr_B = buffer_B;
        }
        MPI_Bcast(buff_ptr_B, nb * col_B, MPI_DOUBLE, owner_B, col_comm);
        for (int i = 0; i < row_A; ++i) {
            for (int l = 0; l < nb; ++l) {
                for (int j = 0; j < col_B; ++j) {
                    C_local[i * col_B + j] += buff_ptr_A[i * nb + l] * buff_ptr_B[l * col_B + j];
                }
            }
        }
    }
}

static void MPI_Datatype_Deleter(MPI_Datatype *p) {
    MPI_Type_free(p);
}

void matrix_gather(MPI_Comm comm, const double M[], const int n_row, const int n_col, double result[], int root) {
    int rank, size, pr, pc;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    if (rank == root) {
        int dims[2], periods[2], my_coords[2];
        MPI_Cart_get(comm, 2, dims, periods, my_coords);
        pr = dims[0];
        pc = dims[1];
        MPI_Datatype subarray, resized_subarray;
        int sizes[] = {pr * n_row, pc * n_col};
        int sub_sizes[] = {n_row, n_col};
        int starts[] = {0, 0};
        MPI_Type_create_subarray(2, sizes, sub_sizes, starts, MPI_ORDER_C, MPI_DOUBLE, &subarray);
        MPI_Type_commit(&subarray);
        deleted_unique_ptr<MPI_Datatype> subarray_owner(&subarray, MPI_Datatype_Deleter);
        MPI_Type_create_resized(subarray, 0, sizeof(double) * n_col, &resized_subarray);
        MPI_Type_commit(&resized_subarray);
        deleted_unique_ptr<MPI_Datatype> resized_subarray_owner(&resized_subarray, MPI_Datatype_Deleter);
        int displs[size];
        int recvcounts[size];
        for (int i = 0; i < size; ++i) {
            int coords[2];
            MPI_Cart_coords(comm, i, 2, coords);
            displs[i] = coords[0] * pc * n_row + coords[1];
            recvcounts[i] = 1;
        }
        MPI_Gatherv(M, n_col * n_row, MPI_DOUBLE, result, recvcounts, displs, resized_subarray, root, comm);
        // MPI_Type_free(&resized_subarray);
        // MPI_Type_free(&subarray);
    } else {
        MPI_Gatherv(M, n_col * n_row, MPI_DOUBLE, nullptr, nullptr, nullptr, MPI_DOUBLE, root, comm);
    }
}

void matrix_scatter(MPI_Comm comm, const double M[], const int n_row, const int n_col, double result[], int root) {
    int rank, size, pr, pc;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);
    if (rank == root) {
        int dims[2], periods[2], my_coords[2];
        MPI_Cart_get(comm, 2, dims, periods, my_coords);
        pr = dims[0];
        pc = dims[1];
        MPI_Datatype subarray, resized_subarray;
        int sizes[] = {pr * n_row, pc * n_col};
        int sub_sizes[] = {n_row, n_col};
        int starts[] = {0, 0};
        MPI_Type_create_subarray(2, sizes, sub_sizes, starts, MPI_ORDER_C, MPI_DOUBLE, &subarray);
        MPI_Type_commit(&subarray);
        deleted_unique_ptr<MPI_Datatype> subarray_owner(&subarray, MPI_Datatype_Deleter);
        MPI_Type_create_resized(subarray, 0, sizeof(double) * n_col, &resized_subarray);
        MPI_Type_commit(&resized_subarray);
        deleted_unique_ptr<MPI_Datatype> resized_subarray_owner(&resized_subarray, MPI_Datatype_Deleter);
        int displs[size];
        int sendcounts[size];
        for (int i = 0; i < size; ++i) {
            int coords[2];
            MPI_Cart_coords(comm, i, 2, coords);
            displs[i] = coords[0] * pc * n_row + coords[1];
            sendcounts[i] = 1;
        }
        MPI_Scatterv(M, sendcounts, displs, resized_subarray, result, n_row * n_col, MPI_DOUBLE, root, comm);
        // MPI_Type_free(&resized_subarray);
        // MPI_Type_free(&subarray);
    } else {
        MPI_Scatterv(nullptr, nullptr, nullptr, MPI_UNDEFINED, result, n_col * n_row, MPI_DOUBLE, root, comm);
    }
}

void print_matrix(const double M[], const int n_row, const int n_col) {
    for (int i = 0; i < n_row; ++i) {
        for (int j = 0; j < n_col; ++j) {
            std::cout << M[i * n_col + j] << '\t';
        }
        std::cout << std::endl;
    }
}

void init_matrix(double M[], const int n_row, const int n_col) {
#ifdef DEBUG
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    std::default_random_engine re((std::random_device()) ());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    double n = 1;
#endif
    for (int i = 0; i < n_row; ++i) {
        for (int j = 0; j < n_col; ++j) {
#ifdef DEBUG
            M[i * n_col + j] = rank;
#else
            M[i * n_col + j] = dist(re);
#endif
        }
    }
}

double validate_matrix(const double M[], const double N[], const int n_row, const int n_col) {
    double max_err = 0.0;
    for (int i = 0; i < n_row; ++i) {
        for (int j = 0; j < n_col; ++j) {
            size_t idx = i * n_col + j;
            double err = std::abs(M[idx] - N[idx]);
            if (err > max_err) max_err = err;
        }
    }
    return max_err;
}

void matmul(const double A[], const double B[], double C[], const int row_A, const int col_A, const int col_B) {
    memset(C, 0, row_A * col_B * sizeof(double));
    for (int i = 0; i < row_A; ++i) {
        for (int k = 0; k < col_A; ++k) {
            for (int j = 0; j < col_B; ++j) {
                C[i * col_B + j] += A[i * col_A + k] * B[k * col_B + j];
            }
        }
    }
}