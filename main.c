#include <stdio.h>
#include <stdlib.h>
#include <time.h>



size_t min(size_t a, size_t b){
    return a < b ? a : b;
}

void matmul(const float* A, const float *B, float * C, size_t rows_C, size_t cols_C, size_t rows_B_cols_A){
    for (size_t i = 0; i < rows_C; i++ ){
        for (size_t j = 0; j < cols_C; j++){
            for (size_t k = 0; k < rows_B_cols_A; k++){            
                C[i * rows_C + j] += A[i * rows_C + k] * B[j *cols_C + k]; 
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
                            C[ii * rows_C + jj] += A[ii * rows_C + kk] * B[jj *cols_C + kk]; 
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
                            sum += A[ii * rows_C + kk] * B[jj *cols_C + kk]; 
                        }
                        C[ii * rows_C + jj] += sum;
                    }
                }
            }
        }
    }
}

void tiled_matmul_me(const float* A, const float *B, float * C, size_t rows_C, size_t cols_C, size_t rows_B_cols_A, size_t tile_size){
    for (size_t idx_row_C = 0; idx_row_C < rows_C; idx_row_C += tile_size){
        for (size_t idx_col_C = 0; idx_col_C < cols_C; idx_col_C += tile_size){
            for (size_t idx_rowB_colA = 0; idx_rowB_colA < rows_B_cols_A; idx_rowB_colA += tile_size){
                for (size_t idx_tile_C_row = idx_row_C; idx_tile_C_row < idx_row_C + tile_size && idx_tile_C_row < rows_C; idx_tile_C_row++){
                    for (size_t idx_tile_C_col = idx_col_C; idx_tile_C_col <idx_col_C + tile_size && idx_tile_C_col < cols_C; idx_tile_C_col++){
                        float dot_product = 0;
                        for (size_t idx_tile_inner = idx_rowB_colA; idx_tile_inner < idx_rowB_colA + tile_size && idx_tile_inner < rows_B_cols_A; idx_tile_inner++){
                            size_t offset_A = idx_tile_C_row * cols_C + idx_tile_inner;
                            size_t offset_B = idx_tile_C_col * cols_C + idx_tile_inner;
                            dot_product += A[offset_A] * B[offset_B];
                            // dot_product += A[idx_tile_C_row * cols_C + idx_tile_inner] * B[idx_tile_C_col * cols_C + idx_tile_inner];
                        }
                        size_t offset_C = idx_tile_C_row * rows_C + idx_tile_C_col;
                        C[offset_C] += dot_product;
                        // C[idx_tile_C_row * rows_C + idx_tile_C_col] += dot_product;
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


void matmul_outer_product(const float * A, const float * B, float * C, size_t M, size_t N, size_t K, size_t size_tile){
    for (size_t tile_start_m = 0; tile_start_m < M; tile_start_m += size_tile){
        for (size_t tile_start_n = 0; tile_start_n < N; N += size_tile){
            for (size_t tile_start_k = 0; tile_start_k < K; tile_start_k += size_tile){
                
                for (size_t idx_m = tile_start_m; idx_m < tile_start_m + size_tile && idx_m < M; idx_m++){
                    for (size_t idx_n = tile_start_n; idx_n < tile_start_n + size_tile && idx_n < N; idx_n++){
                        
                    }
                }
            }
        }
    }
}


void initialise_large_matrices(float *A_large, float * B_large, float * C_large){
    srand(42);
    for (size_t i = 0; i < 1024 * 1024; i++){
        A_large[i] = rand() % 100;
        B_large[i] = rand() % 100;
        C_large[i] = 0.0;
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

    // float A[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    // float B[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    // float C[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    // matmul(A, B, C, 4, 4, 4);

    // for (size_t i = 0; i < 4; i++)
    // {
    //     printf("\n");
    //     for (size_t j = 0; j < 4; j++){
    //         printf("%f\t", C[i * 4 + j]);
    //     }
    //     printf("\n");
    // }

    // printf("\n");
    // printf("tiled_matmul_cp2\n");
    
    // float AT[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    // float BT[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    // float CT[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    // tiled_matmul_cp2(AT, BT, CT, 4, 4, 4, 2);

    // for (size_t i = 0; i < 4; i++)
    // {
    //     printf("\n");
    //     for (size_t j = 0; j < 4; j++){
    //         printf("%f\t", CT[i * 4 + j]);
    //     }
    //     printf("\n");
    // }

    // printf("\n");
    // printf("tiled_matmul_me\n");
    
    // float ATM[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    // float BTM[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    // // float BTM[] = {1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4};
    // float CTM[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    // tiled_matmul_me(ATM, BTM, CTM, 4, 4, 4, 2);

    // for (size_t i = 0; i < 4; i++)
    // {
    //     printf("\n");
    //     for (size_t j = 0; j < 4; j++){
    //         printf("%f\t", CTM[i * 4 + j]);
    //     }
    //     printf("\n");
    // }
    
    
    // exit(0);


    printf("executing matmul now ...\n");
    float * LA = malloc(1024 * 1024 * sizeof(float));
    float * LB = malloc(1024 * 1024 * sizeof(float));
    float * LC = malloc(1024 * 1024 * sizeof(float));
    initialise_large_matrices(LA, LB, LC);

    clock_t start = clock();
    matmul(LB, LB, LC, 1024, 1024, 1024);
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time spent on matmul: %f seconds\n", time_spent);


    printf("executing matmulcp now ...\n");
    initialise_large_matrices(LA, LB, LC);
    start = clock();
    tiled_matmul_cp(LA, LB, LC, 1024, 1024, 1024, 256);
    end = clock();
    time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time spent on tiled_matmul_cp: %f seconds\n", time_spent);
    
    printf("executing matmulcp2 now ...\n");
    initialise_large_matrices(LA, LB, LC);
    start = clock();
    tiled_matmul_cp2(LA, LB, LC, 1024, 1024, 1024, 256);
    end = clock();
    time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time spent on tiled_matmul_cp2: %f seconds\n", time_spent);


    printf("executing matmulme now ...\n");
    initialise_large_matrices(LA, LB, LC);
    start = clock();
    tiled_matmul_me(LA, LB, LC, 1024, 1024, 1024, 256);
    end = clock();
    time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Time spent on tiled_matmul_me: %f seconds\n", time_spent);



    free(LA);
    free(LB);
    free(LC);

}



// gcc -o main main.c -lm
// ./main