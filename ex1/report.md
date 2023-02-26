# **Performance evaluation of a single core**

## Problem Description

The goal of this report is to analyse the impact that the memory hierarchy has in a processor's performance, especially when a significant amount of data needs to be read, memory access can be particularly challenging. Considering this, optimizing the code is crucial to minimize the number of processor cache failures.

The problem consists in analysing the matrix multiplication of square matrices. Given two matrices A and B of n x n, the result of the multiplication of these two should be a matrix C of size n x n, defined as: 

``C[i][j] = sum(A[i][k] * B[k][j]) for k = 0 to n-1``, where C[i][j] represents the element in the ith row and jth column of matrix C.

To calculate each element of C, we take the dot product of the ith row of A with the jth column of B. This is done by multiplying each corresponding element of the row and column together and then summing the results.

To solve this problem we tested the perfomance of two different approaches - Basic Multiplication and Line Multiplication - using C++ language. In C++, multi-dimensional arrays are stored in memory in row-major order, with all elements of a row being stored in contiguous memory positions. This means that the last index of the array varies fastest in memory, followed by the second-to-last index, and so on.


## Algorithms Explanation

### Basic Multiplication
 This implementation multiplies each row of a matrix A for each column of matrix B, by following this method the computation of each cell in matrix C is done individually.

During each iteration, the algorithm accesses row i of matrix A and every element of matrix B to compute a row of matrix C, which is then written to memory.

 **PseudoCode**

 ```
    for (i=1; i<n; i++)
        for (j=1; j<n; j++)
            for (k=1; k<n; k++)
                c[i,j] += a[i,k] * b[k,j]

````

### Line Multiplication

The second implementation is more efficient than the first one, because it does the multiplicaiton line by line. This algorithm deviates from the steps performed in algebraic multiplication and instead leverages the architecture and operational style of the processor.

 **PseudoCode**

 ``` 
    for (i=1; i<n; i++)
        for (k=1; k<n; k++)
            for (j=1; j<n; j++)
                c[i,j] += a[i,k] * b[k,j]
````

The algorithm assumes that matrix C's elements are initialized to 0 and works by multiplying elements of a row of matrix A with every row of matrix B. This approach fills in one row of matrix C at a time, instead of computing all values at once.


## Performance Metrics
For all experiments, the execution time was recorded, as well as cache performance using the PAPI library, namely Data Cache Misses (DCM) on L1 an L2 caches. The performance metrics chosen to compare different algorithms with different sized matrices were DCM per FLOP and FLOP per second.

## Results

## Conclusions

