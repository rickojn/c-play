#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define M 10000
#define N 256
#define K 784

size_t min(size_t a, size_t b){
    return a < b ? a : b;
}

void matmul(const float* A, const float *B, float * C, size_t rows_C, size_t cols_C, size_t rows_B_cols_A){
    for (size_t i = 0; i < rows_C; i++ ){
        for (size_t j = 0; j < cols_C; j++){
            for (size_t k = 0; k < rows_B_cols_A; k++){            
                C[i * cols_C + j] += A[i * rows_B_cols_A + k] * B[j * cols_C + k]; 
            }
        }
    }
}


void tiled_matmul_cp(const float* A, const float *B, float * C, size_t rows_C, size_t cols_C, size_t rows_B_cols_A, size_t tile_size){
    for (size_t i = 0; i < rows_C; i += tile_size){
        for (size_t j = 0; j < cols_C; j += tile_size){
            for (size_t k = 0; k < rows_B_cols_A; k += tile_size){
                for (size_t ii = i; ii < i + tile_size && ii < rows_C; ii++){
                    for (size_t jj = j; jj < j + tile_size && jj < cols_C; jj++){
                        for (size_t kk = k; kk < k + tile_size && kk < rows_B_cols_A; kk++){
                            C[ii * cols_C + jj] += A[ii * rows_B_cols_A + kk] * B[jj *cols_C + kk]; 
                        }
                    }
                }
            }
        }
    }
 
}


void tiled_matmul_cp2(const float* A, const float *B, float * C, size_t rows_C, size_t cols_C, size_t rows_B_cols_A, size_t tile_size){
    for (size_t i = 0; i < rows_C; i += tile_size){
        for (size_t j = 0; j < cols_C; j += tile_size){
            for (size_t k = 0; k < rows_B_cols_A; k += tile_size){
                for (size_t ii = i; ii < i + tile_size && ii < rows_C; ii++){
                    for (size_t jj = j; jj < j + tile_size && jj < cols_C; jj++){
                        float sum = 0;
                        for (size_t kk = k; kk < k + tile_size && kk < rows_B_cols_A; kk++){
                            sum += A[ii * rows_B_cols_A + kk] * B[jj *cols_C + kk]; 
                        }
                        C[ii * cols_C + jj] += sum;
                    }
                }
            }
        }
    }
}



/*

1 1    1 1
2 2    2 2

3 3
6 6
sum of outer products

first column by first row:
1  (+)   1 1
2

=

1 1
2 2

second column by second row:

1  (+)   2 2
2

2 2 
4 4 


sum of outer products:

2 2 + 1 1
4 4 + 2 2
=
3 3 
6 6


*/


void outer_product_sum_matmul(const float * A, const float * B, float * C, size_t m, size_t n, size_t k, size_t size_tile){
    for (size_t tile_start_m = 0; tile_start_m < m; tile_start_m += size_tile){
        for (size_t tile_start_n = 0; tile_start_n < n; tile_start_n += size_tile){
            for (size_t tile_start_k = 0; tile_start_k < k; tile_start_k += size_tile){


                for (size_t idx_m = tile_start_m; idx_m < tile_start_m + size_tile && idx_m < m; idx_m++){
                    for (size_t idx_n = tile_start_n; idx_n < tile_start_n + size_tile && idx_n < n; idx_n++){
                        float sum = 0.0;
                        // each k (col of A and row of B)
                        size_t offset_a = idx_m * k;
                        size_t offset_b = idx_n * n;
                        for (size_t idx_k = tile_start_k; idx_k < tile_start_k + size_tile && idx_k < k; idx_k++){
                            sum += A[offset_a + idx_k] * B[idx_k + offset_b];
                            }
                        C[idx_m * n + idx_n] += sum;
                    }
                }


            }
        }
    }
}




void initialise_large_matrices(float *A_large, float * B_large, float * C_large){
    srand(42);
    for (size_t i = 0; i < M * K; i++){
        A_large[i] = rand() % 100;
    }
    for (size_t i = 0; i < K * N; i++){
        B_large[i] = rand() % 100;
    }
    for (size_t i = 0; i < M * N; i++){
        C_large[i] = 0.0f;
    }

}

void check_result(const float * ref_result, const float * result, size_t m, size_t n){
    for (size_t idx = 0; idx < m * n; idx++){
        if (ref_result[idx] != result[idx]){
            printf("this does not aggree with naive matmul\n");
            printf("ref_result[%zu] = %f  result[%zu] = %f\n", idx, ref_result[idx], idx, result[idx]);
            return;
        }
    }
}



int main() {
    printf("matmul\n");

    /*
    1111  1234
    2222  1234
    3333  1234
    4444  1234

    4,8,12,16
    8,16,24,32
    12,24,36,48
    16,32,48,64
    */

    float A[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // row major
    float B[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // column major
    float C[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    matmul(A, B, C, 4, 4, 4);

    for (size_t i = 0; i < 4; i++)
    {
        printf("\n");
        for (size_t j = 0; j < 4; j++){
            printf("%f\t", C[i * 4 + j]);
        }
        printf("\n");
    }

    printf("\n");
    printf("tiled_matmul_cp2\n");
    
    float AT[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    float BT[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // column major
    float CT[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    tiled_matmul_cp2(AT, BT, CT, 4, 4, 4, 2);

    for (size_t i = 0; i < 4; i++)
    {
        printf("\n");
        for (size_t j = 0; j < 4; j++){
            printf("%f\t", CT[i * 4 + j]);
        }
        printf("\n");
    }


    printf("\n");
    printf("outer product matmul:\n");

    float ATP[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // row major
    float BTP[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // column major
    float CTP[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // row major
    
    outer_product_sum_matmul(ATP, BTP, CTP, 4, 4, 4, 2);

    printf("\n");
    for (size_t i = 0; i < 4; i++)
    {
        printf("\n");   
        for (size_t j = 0; j < 4; j++){
            printf("%f\t", CTP[i * 4 + j]);
        }
        printf("\n");
    }
    
    
    // exit(0);


    printf("executing matmul now ...\n");
    float * LA = malloc(M * K * sizeof(float));
    float * LB = malloc(K * N * sizeof(float));
    float * LC = malloc(M * N * sizeof(float));
    initialise_large_matrices(LA, LB, LC);
    float * ref_C = calloc(M * N, sizeof(float));

    clock_t start = clock();
    matmul(LA, LB, ref_C, M, N, K);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time spent on matmul: %f seconds\n", time_spent);


    printf("executing matmulcp now ...\n");
    initialise_large_matrices(LA, LB, LC);
    start = clock();
    tiled_matmul_cp(LA, LB, LC, M, N, K, 256);
    end = clock();
    time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time spent on tiled_matmul_cp: %f seconds\n", time_spent);
    check_result(ref_C, LC, 1024, 1024);
    
    printf("executing matmulcp2 now ...\n");
    initialise_large_matrices(LA, LB, LC);
    start = clock();
    tiled_matmul_cp2(LA, LB, LC, M, N, K, 256);
    end = clock();
    time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time spent on tiled_matmul_cp2: %f seconds\n", time_spent);
    check_result(ref_C, LC, 1024, 1024);



    printf("executing outer product matmul now ...\n");
    initialise_large_matrices(LA, LB, LC);
    start = clock();
    outer_product_sum_matmul(LA, LB, LC, M, N, K, 256);
    end = clock();
    time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time spent on outer product matmul2: %f seconds\n", time_spent);
    check_result(ref_C, LC, 1024, 1024);



    free(LA);
    free(LB);
    free(LC);

}


// gcc -o main main.c -lm
// ./main