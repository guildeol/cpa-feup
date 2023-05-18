#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

extern "C" {
    void dlacpy_(const char *uplo, const int *m, const int *n, const double *A, const int *lda, double *B, const int *ldb);
    void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha,
                const double *A, const int *lda, const double *B, const int *ldb, const double *beta, double *C, const int *ldc);
}

void pdgemm(int m, int n, int k, int nb, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc,
            int nprow, int npcol, int myrow, int mycol, MPI_Comm comm_row, MPI_Comm comm_col, double *work1, double *work2)
{
    int myrowcol = myrow * npcol + mycol;
    int iwrk, jwrk;
    int m_a, n_a, m_b, n_b, m_c, n_c;
    double one = 1.0;
    double zero = 0.0;
    MPI_Status status;

    m_a = ((m - 1) / nb + 1) * nb;
    n_a = k;
    m_b = k;
    n_b = ((n - 1) / nb + 1) * nb;
    m_c = ((m - 1) / nb + 1) * nb;
    n_c = ((n - 1) / nb + 1) * nb;

    for (int jj = 0; jj < n_c; jj += nb)
    {
        jwrk = n_c - jj;
        if (jwrk > nb)
            jwrk = nb;

        if (mycol == (jj / nb) % npcol)
        {
            dlacpy_("General", &m_a, &jwrk, &a[0], &lda, &work1[0], &m_a);
        }

        MPI_Bcast(&work1[0], m_a * jwrk, MPI_DOUBLE, (jj / nb) % npcol, comm_row);


        if (myrow == (jj / nb) % nprow)
        {
            dlacpy_("General", &jwrk, &n_b, &b[jj / nb * k], &ldb, &work2[0], &jwrk);
        }

        MPI_Bcast(&work1[0], m_a * jwrk, MPI_DOUBLE, (jj / nb) % npcol, comm_row);


        dgemm_("No transpose", "No transpose", &m_c, &n_c, &jwrk, &alpha, &work1[0], &m_a, &work2[0], &jwrk, &beta, &c[0], &ldc);
    }
}


int main(int argc, char *argv[])
{
    int m, n, k;
    int nprow, npcol, nb;
    int me, np;
    int myrow, mycol;
    int lda, ldb, ldc;
    double alpha, beta;
    double *a, *b, *c;
    double *work1, *work2;
    MPI_Comm comm, comm_row, comm_col;
    MPI_Group orig_group, row_group, col_group;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &me);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    m = 4;  // Number of rows in matrix A
    n = 4;  // Number of columns in matrix B
    k = 4;  // Number of columns in matrix A / rows in matrix B
    nb = 2; // Block size

    nprow = 2; // Number of rows in process grid
    npcol = 2; // Number of columns in process grid

    lda = ((m - 1) / nb + 1) * nb;
    ldb = ((k - 1) / nb + 1) * nb;
    ldc = ((m - 1) / nb + 1) * nb;

    if (me == 0)
    {
        a = (double *)malloc(lda * k * sizeof(double));
        b = (double *)malloc(ldb * n * sizeof(double));
        c = (double *)malloc(ldc * n * sizeof(double));

        // Initialize matrices A and B with some values
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < k; ++j)
            {
                a[i + j * lda] = i * k + j + 1;
            }
        }

        for (int i = 0; i < k; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                b[i + j * ldb] = i * n + j + 1;
            }
        }
    }

    // Create virtual 2D grid communicator
    int dims[2] = {nprow, npcol};
    int periods[2] = {1, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm);

    // Determine my row and column rank
    MPI_Comm_rank(comm, &me);
    MPI_Comm_rank(comm, &myrow);
    MPI_Comm_rank(comm, &mycol);

    // Create row and column communicators
    int remain_dims_row[2] = {0, 1};
    int remain_dims_col[2] = {1, 0};
    MPI_Cart_sub(comm, remain_dims_row, &comm_row);
    MPI_Cart_sub(comm, remain_dims_col, &comm_col);

    // Allocate workspace arrays
    work1 = (double *)malloc(lda * nb * sizeof(double));
    work2 = (double *)malloc(nb * ldb * sizeof(double));

    alpha = 1.0;
    beta = 0.0;

    // Perform parallel matrix multiplication
    pdgemm(m, n, k, nb, alpha, a, lda, b, ldb, beta, c, ldc, nprow, npcol, myrow, mycol, comm_row, comm_col, work1, work2);

    if (me == 0)
    {
        // Print the result matrix C
        printf("Result matrix C:\n");
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                printf("%6.2f ", c[i + j * ldc]);
            }
            printf("\n");
        }

        free(a);
        free(b);
        free(c);
    }

    free(work1);
    free(work2);

    MPI_Finalize();

    return 0;
}
