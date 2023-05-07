#include <stdint.h>
#include <iostream>
#include <cmath>
#include <mpi.h>
#include <omp.h>

using namespace std;

#define GRID_2D  (2)

double **A, **B, **C;

double **matrix_alloc(uint64_t n)
{
    double **M;

	// Allocate n * n contiguous items
	double *p = (double*)calloc((n * n), sizeof(double));
	if (!p)
		return NULL;

	// Allocate row pointers
	M = (double **)malloc(n * sizeof(double *));
	if (!M) {
		free(p);
		return NULL;
	}

	// Set up the pointers into the contiguous memory
	for (int i = 0; i < n; i++)
		M[i] = &(p[i * n]);

	return M;
}

void matrix_create(uint64_t n)
{
    A = matrix_alloc(n);
    B = matrix_alloc(n);
    C = matrix_alloc(n);

    for (uint64_t i = 0; i < n; i++)
        for (uint64_t j = 0; j < n; j++)
            A[i][j] = 1.0;

    for (uint64_t  i = 0; i < n; i++)
        for (uint64_t j = 0; j < n; j++)
            B[i][j] = (double)(i + 1);
}

void matrix_mult(double **a, double **b, double ***c, uint64_t block_size)
{
    uint64_t i, j, k;

    #pragma omp parallel num_threads(4)
    {
        for (i = 0; i < block_size; i++)
        {
            for (j = 0; j < block_size; j++)
            {
                double val = 0;
                for (k = 0; k < block_size; k++)
                {
                    val += a[i][k] * b[k][j];
                }

                (*c)[i][j] = val;
            }
        }
    }
}

void matrix_print(double **M, uint64_t n, uint64_t i, uint64_t j, uint64_t block_size)
{
    uint64_t row, col;

    for (row = 0; row < n; row++) {
        printf("|");
        for (col = 0; col < n; col++) {
            if (row >= i && row < i + block_size && col >= j && col < j + block_size)
            {
                printf(" %2d ", (int)M[row][col]);
            }
            else
            {
                printf("    ");
            }
        }
        printf("|\n");
    }
}

void matrix_destroy(double **M)
{
	free(&M[0][0]);
	free(M);
}

int main(int argc, char *argv[])
{
    double begin, end;

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = atoi(argv[1]);  // Matrix dimension (n x n)
    int q = sqrt(size);     // Number of processes should be a perfect square

    // Check if the number of processes is a perfect square
    if (q * q != size)
    {
        if (rank == 0)
            cerr << "Number of processes must be a perfect square!" << endl;
        
        MPI_Finalize();
        return 1;
    }

    int block_size = n / q;

    // Create Cartesian topology
    int dims[GRID_2D] = {q, q};
    int periods[GRID_2D] = {1, 1};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, GRID_2D, dims, periods, 0, &cart_comm);

    int coords[GRID_2D];
    MPI_Cart_coords(cart_comm, rank, GRID_2D, coords);

    if (rank == 0)
    {
        matrix_create(n);
    }

    double **local_A = matrix_alloc(block_size);
    double **local_B = matrix_alloc(block_size);
    double **local_C = matrix_alloc(block_size);

    // Create datatype to describe the subarrays of the global array
	int matrix_dimensions[GRID_2D] = { n, n };
	int block_dimensions[GRID_2D] = { block_size, block_size };
	int starts[GRID_2D] = { 0,0 };
	MPI_Datatype type, subarrtype;
	MPI_Type_create_subarray(GRID_2D, matrix_dimensions, block_dimensions, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
	MPI_Type_create_resized(type, 0, block_size * sizeof(double), &subarrtype);
	MPI_Type_commit(&subarrtype);

	// Scatter the array to all processors
	int* sendCounts = (int*)malloc(sizeof(int) * size);
	int* displacements = (int*)malloc(sizeof(int) * size);

	if (rank == 0)
    {
		for (int i = 0; i < size; i++)
			sendCounts[i] = 1;

        for (int i = 0; i < size; i++) 
            displacements[i] = (i / q) * q * block_size + i % q;
	}

	double *globalptrA = NULL;
	double *globalptrB = NULL;
	double *globalptrC = NULL;
	if (rank == 0)
    {
		globalptrA = &(A[0][0]);
		globalptrB = &(B[0][0]);
		globalptrC = &(C[0][0]);

        begin = MPI_Wtime();
	}

	MPI_Scatterv(globalptrA, sendCounts, displacements, subarrtype, &local_A[0][0],
		         (n * n) / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatterv(globalptrB, sendCounts, displacements, subarrtype, &local_B[0][0],
		         (n * n) / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    int left, right, up, down;
    int coord[2];

	double **intermediate_matrix = matrix_alloc(block_size);
	for (int k = 0; k < q; k++) {
		matrix_mult(local_A, local_B, &intermediate_matrix, block_size);

		for (int i = 0; i < block_size; i++) {
			for (int j = 0; j < block_size; j++)
            {
				local_C[i][j] += intermediate_matrix[i][j];
			}
		}
		// Shift A once (left) and B once (up)
		MPI_Cart_shift(cart_comm, 1, 1, &left, &right);
		MPI_Cart_shift(cart_comm, 0, 1, &up, &down);
		MPI_Sendrecv_replace(&(local_A[0][0]), block_size * block_size, MPI_DOUBLE,
                             left, 1, right, 1, cart_comm, MPI_STATUS_IGNORE);

		MPI_Sendrecv_replace(&(local_B[0][0]), block_size * block_size, MPI_DOUBLE,
                             up, 1, down, 1, cart_comm, MPI_STATUS_IGNORE);
	}
	
	// Gather results
	MPI_Gatherv(&(local_C[0][0]), (n * n) / size, MPI_DOUBLE,
		        globalptrC, sendCounts, displacements, subarrtype,
		        0, MPI_COMM_WORLD);
                
    if (rank == 0)
    {
        end = MPI_Wtime();
        printf("Ellapsed time: %.2f seconds\n", end - begin);
        printf("C[0][0] = %lf\n", C[0][0]);
        //matrix_print(C, n, 0, 0, n);
    }

    MPI_Finalize();

    return 0;
}