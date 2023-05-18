#include <mpi.h>
#include <cstdlib>
#include <cstdio>
#include <stdlib.h>

#include <lapack.h>
#undef LAPACK_FORTRAN_STRLEN_END


//Macros for column major indexing

#define A(i,j) (a[j*lda+i])
#define B(i,j) (b[j*ldb+i])
#define C(i,j) (c[j*ldc+i])

#define min(x,y) ((x)<(y) ? (x):(y))

//Constants passed to BLAS call
int i_one=1;
double d_one=1.0, d_zero=0.0;

void RING_Bcast(double *buf, int count, MPI_Datatype type, int root, MPI_Comm comm);

void RING_SUM( double *buf, int count, MPI_Datatype type, int root, MPI_Comm comm, double *work );

// m, n and k - global matrix dimensions
// nb - panel width
// alpha and beta - multiplication constants
// m_a[], n_a[], m_b[], n_b[], m_c[] and n_c[] - dimensions of blocks of A, B and C
// lda, ldb, ldc - leading dimension of local arrays that hold local portions of matrices A, B, C
//*a, *b, *c - arrays that hold local parts of A, B, C
// *work1, *work2 - work arrays
// comm_row - communicator for this row of nodes
// comm_col - communicator for this column of nodes
void pdgemm(int m, int n, int k, int nb, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc,
    int m_a[], int n_a[], int m_b[], int n_b[], int m_c[], int n_c[], MPI_Comm comm_row, MPI_Comm comm_col, double *work1, double *work2 )
{
    // row and column indexes
    int myrow, mycol;
    // number of node rows and columns
    int nprow, npcol;
    // index variables
    int i, j, kk, iwrk;
    // index of row and column that hold current row and column, for rank-1 update
    int icurrow, icurcol;
    // local index (on icurrow and icurcol) of row and column for rank-1 update
    int ii, jj;
    // temporary pointer used in pdgemm_abt
    double *temp;   
    double *p;

    // get myrow, mycol
    MPI_Comm_rank( comm_row, &mycol );
    MPI_Comm_rank( comm_col, &myrow );

    // scale local block of C
    for ( j=0; j<n_c[mycol]; j++ ){
        for ( i=0; i<m_c[myrow]; i++ ){
            C( i,j ) = beta * C( i,j );
        }
    }
    
    icurrow = 0;
    icurcol = 0;
    ii = jj = 0;

    // malloc temp space for summation
    temp = (double *) malloc(m_c[myrow]*nb*sizeof(double) );

    for ( kk=0; kk<k; kk+=iwrk) {
        iwrk = min( nb, m_b[icurrow]-ii );
        iwrk = min( iwrk, n_a[icurcol]-jj );

        // pack current iwrk columns of A into work1 
        if ( mycol == icurcol ) {
            dlacpy_( "General", &m_a[myrow], &iwrk, &A(0, jj), &lda, work1, &m_a[myrow] );
        }
        
        // pack current iwrk rows of B into work2 
        if ( myrow == icurrow ) {
            dlacpy_( "General", &iwrk, &n_b[mycol], &B(ii, 0), &ldb, work2, &iwrk );
        }
        
        // broadcast work1 and work2
        RING_Bcast( work1 , m_a[ myrow ]*iwrk, MPI_DOUBLE, icurcol, comm_row );
        RING_Bcast( work2 , n_b[ mycol ]*iwrk, MPI_DOUBLE, icurrow, comm_col );

        // update local block 
        dgemm_( "No transpose", "No transpose", &m_c[ myrow ], &n_c[ mycol ], &iwrk, &alpha, work1, &m_b[ myrow ], work2, &iwrk, &d_one, c, &ldc );

        // update icurcol, icurrow, ii, jj
        ii += iwrk; jj += iwrk;
        if ( jj>=n_a[ icurcol ] ) { 
            icurcol++; jj = 0;
        }
        if ( ii>=m_b[ icurrow ] ) { 
            icurrow++; ii = 0;
        }
    }

    free( temp );
}

void RING_Bcast( double *buf, int count, MPI_Datatype type, int root, MPI_Comm comm )
{
    int me, np;

    MPI_Status status;
    MPI_Comm_rank( comm, me ); MPI_Comm_size( comm, np );

    if ( me != root){
        MPI_Recv( buf, count, type, (me-1+np)%np, MPI_ANY_TAG, comm );
    }
    if ( ( me+1 )%np != root ) {
        MPI_Send( buf, count, type, (me+1)%np, 0, comm );
    }

    int icurrow = 0;
    int icurcol = 0;
    int ii = jj = 0;

    double *p;

    // malloc temp space for summation
    temp = (double *) malloc( m_c[myrow]*nb*sizeof(double) );

    // loop over all column panels of C 
    for ( kk=0; kk<k; kk+=iwrk) {
        iwrk = min( nb, m_b[ icurrow ]-ii );
        iwrk = min( iwrk, n_c[ icurcol ]-jj );

        // pack current iwrk rows of B into work2
        if ( myrow == icurrow ) {
            dlacpy_( "General", &iwrk, &n_b[ mycol ], &B( ii, 0 ), &ldb, work2, &iwrk );
        }
        
        // broadcast work2
        RING_Bcast( work2 , n_b[ mycol ]*iwrk, MPI_DOUBLE, icurrow, comm_col );

        // Multiply local block of A times incoming rows of B 
        dgemm_( "No transpose", "Transpose", &m_c[ myrow ], &iwrk, &n_a[ mycol ], &alpha, a, &lda, work2, &iwrk, &d_zero, work1, &m_c[ myrow ] );

        //Sum to node that holds current columns of C
        RING_SUM( work1, m_c[ myrow ]*iwrk, MPI_DOUBLE, icurcol, comm_row, temp );

        // Add to current columns of C
        if ( mycol == icurcol ) {
            p = work1;
            for (int j=jj; j<jj+iwrk; j++) {
                daxpy_( &m_c[ myrow ], &d_one, p, &i_one, &C( 0,j ), &i_one );
                p += m_c[ myrow ];
            }
        }

        // update icurcol, icurrow, ii, jj
        ii += iwrk; jj += iwrk;
        if ( jj>=n_c[ icurcol ] ) { 
            icurcol++;
            jj = 0;
        }
        if ( ii>=m_b[ icurrow ] ) {
            icurrow++;
            ii = 0; 
        }
    }
    free( temp );
}

void RING_SUM( double *buf, int count, MPI_Datatype type, int root, MPI_Comm comm, double *work )
{
    int me, np;

    MPI_Status status;
    MPI_Comm_rank( comm, &me ); MPI_Comm_size( comm, &np );

    if ( me != (root+1)%np ) {
        MPI_Recv( work, count, type, (me-1+np)%np, MPI_ANY_TAG, comm, &status );
        daxpy_( &count, &d_one, work, &i_one, buf, &i_one );
    }

    if ( me != root ) {
        MPI_Send( buf, count, type, (me+1)%np, 0, comm );
    }
}