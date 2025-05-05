#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define M 10000
#define N 256
#define K 784
#define TILE 128
#define INNER_TILE 32
#define NAIVE 0
#define OUTER 0
#define TILED 0
#define L1 1
#define ONLY_LARGE 1



size_t min(size_t a, size_t b){
    return a < b ? a : b;
}

void naive_matmul(const float* A, const float *B, float * C, size_t rows_C, size_t cols_C, size_t rows_B_cols_A){
    for (size_t i = 0; i < rows_C; i++ ){
        for (size_t j = 0; j < cols_C; j++){
            for (size_t k = 0; k < rows_B_cols_A; k++){            
                C[i * cols_C + j] += A[i * rows_B_cols_A + k] * B[j * cols_C + k]; 
            }
        }
    }
}







/*

m = 3 n= 4 k = 2

1 10    1 1 1 1
2 20    2 2 2 2
3 30    

21 21 21 21
42 42 42 42
63 63 63 63

sum of outer products

first column by first row:
1  (x)   1 1 1 1
2
3

=

1 1 1 1
2 2 2 2
3 3 3 3



second column by second row:

10  (x)   2 2 2 2
20
30

20 20 20 20
40 40 40 40
60 60 60 60
 


sum of outer products:

1 1 1 1   20 20 20 20
2 2 2 2 + 40 40 40 40
3 3 3 3   60 60 60 60
=
21 21 21 21
42 42 42 42
63 63 63 63


*/

void outer_product_matmul(const float * a, const float * b, float * c, size_t m, size_t n, size_t k ){
    for (size_t idx_k = 0; idx_k < k; idx_k++){
        for (size_t idx_m = 0; idx_m < m; idx_m++){
            for (size_t idx_n = 0; idx_n < n; idx_n++){
                // size_t offset_a = idx_k * m  + idx_m; // a[m][k] col major 
                // size_t offset_b = idx_k * n + idx_n;  // b[k][n]  row major
                // size_t offset_c = idx_m * n + idx_n; // row major
                // c[offset_c] += a[offset_a] * b[offset_b];
                c[idx_m * n +idx_n] += a[idx_k * m +idx_m] * b[idx_k * n + idx_n];
            }
        }
    }
}

void tiled_matmul(const float * A, const float * B, float * C, size_t m, size_t n, size_t k, size_t size_tile){
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



void l1_tiled_matmul(const float * A, const float * B, float * C, size_t m, size_t n, size_t k, size_t size_outer_tile, size_t size_inner_tile){
    for (size_t idx_m = 0; idx_m < m; idx_m += size_outer_tile){
        for (size_t idx_n = 0; idx_n < n; idx_n += size_outer_tile){
            for (size_t idx_k = 0; idx_k < k; idx_k += size_outer_tile){
                
                for (size_t idx_mm = idx_m; idx_mm < idx_m + size_outer_tile && idx_mm < m; idx_mm += size_inner_tile){
                    for (size_t idx_nn = idx_n; idx_nn < idx_n + size_outer_tile && idx_nn < n; idx_nn += size_inner_tile){
                        for (size_t idx_kk = idx_k; idx_kk < idx_k + size_outer_tile && idx_kk < k; idx_kk += size_inner_tile){

                            for (size_t idx_mmm = idx_mm; idx_mmm < idx_mm + size_inner_tile && idx_mmm < m; idx_mmm++){
                                for (size_t idx_nnn = idx_nn; idx_nnn < idx_nn + size_inner_tile && idx_nnn < n; idx_nnn++){
                                    float sum = 0;
                                    size_t offset_a = idx_mmm * k; // A[idx_mmm][0] row major
                                    size_t offset_b = idx_nnn * k; // B[0][idx_nnn] col major
                                    size_t idx_kkk = idx_kk;
                                    for (;idx_kkk < idx_kk + size_inner_tile && idx_kkk < k; idx_kkk+= 8){
                                        sum += A[offset_a + idx_kkk] * B[offset_b + idx_kkk];
                                        sum += A[offset_a + idx_kkk + 1] * B[offset_b + idx_kkk + 1];
                                        sum += A[offset_a + idx_kkk + 2] * B[offset_b + idx_kkk + 2];
                                        sum += A[offset_a + idx_kkk + 3] * B[offset_b + idx_kkk + 3];
                                        sum += A[offset_a + idx_kkk + 4] * B[offset_b + idx_kkk + 4];
                                        sum += A[offset_a + idx_kkk + 5] * B[offset_b + idx_kkk + 5];
                                        sum += A[offset_a + idx_kkk + 6] * B[offset_b + idx_kkk + 6];
                                        sum += A[offset_a + idx_kkk + 7] * B[offset_b + idx_kkk + 7];
                                    }
                                    for (; idx_kkk < idx_kk + size_inner_tile && idx_kkk < k; idx_kkk++){
                                        sum += A[offset_a + idx_kkk] * B[offset_b +idx_kkk];
                                    }
                                    C[idx_mmm * n + idx_nnn] += sum;
                                }
                            }
                        }
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

    if (!ONLY_LARGE && NAIVE){
        printf("naive matmul ... \n");
        float A[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // row major
        float B[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // column major
        float C[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        naive_matmul(A, B, C, 4, 4, 4);
    
        for (size_t i = 0; i < 4; i++)
        {
            printf("\n");
            for (size_t j = 0; j < 4; j++){
                printf("%f\t", C[i * 4 + j]);
            }
            printf("\n");
        }
    
        printf("\n");
    }

    if (!ONLY_LARGE && OUTER){
        printf("outer product matmul ... \n");
        float AO[] = {1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4}; // column major
        float BO[] = {1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4}; // row major
        float CO[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        outer_product_matmul(AO, BO, CO, 4, 4, 4);
    
        for (size_t i = 0; i < 4; i++)
        {
            printf("\n");
            for (size_t j = 0; j < 4; j++){
                printf("%f\t", CO[i * 4 + j]);
            }
            printf("\n");
        }
    
        printf("\n");
    }
    



    if (!ONLY_LARGE && TILED){
        printf("\n");
        printf("tiled matmul:\n");
    
        float ATP[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // row major
        float BTP[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // column major
        float CTP[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // row major
        
        tiled_matmul(ATP, BTP, CTP, 4, 4, 4, 2);
    
        printf("\n");
        for (size_t i = 0; i < 4; i++)
        {
            printf("\n");   
            for (size_t j = 0; j < 4; j++){
                printf("%f\t", CTP[i * 4 + j]);
            }
            printf("\n");
        }
    }

    if (!ONLY_LARGE && L1){
        printf("\n");
        printf("l1 tiled matmul:\n");
    
        float ATL[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // row major
        float BTL[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4}; // column major
        float CTL[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}; // row major
        
        l1_tiled_matmul(ATL, BTL, CTL, 4, 4, 4, 2, 1);
    
        printf("\n");
        for (size_t i = 0; i < 4; i++)
        {
            printf("\n");   
            for (size_t j = 0; j < 4; j++){
                printf("%f\t", CTL[i * 4 + j]);
            }
            printf("\n");
        }
        
    }


    
    
    // exit(0);


    float * LA = malloc(M * K * sizeof(float));
    float * LB = malloc(K * N * sizeof(float));
    float * LC = malloc(M * N * sizeof(float));
    initialise_large_matrices(LA, LB, LC);
    float * ref_C = calloc(M * N, sizeof(float));
    clock_t start, end;
    double time_spent;

    if (NAIVE)
    {
        printf("Naive .. \n");
        start = clock();
        naive_matmul(LA, LB, ref_C, M, N, K);
        end = clock();
        time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Time spent on matmul: %f seconds\n", time_spent);
    }

    if (TILED){
        printf("executing dot product matmul now with tile %d ...\n", TILE);
        initialise_large_matrices(LA, LB, LC);
        start = clock();
        tiled_matmul(LA, LB, LC, M, N, K, TILE);
        end = clock();
        time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Time spent on tiled matmul: %f seconds\n", time_spent);
        if (NAIVE){
            check_result(ref_C, LC, 1024, 1024);
        }
    }

    if (OUTER){
        printf("executing outer product matmul now ...\n");
        initialise_large_matrices(LA, LB, LC);
        start = clock();
        outer_product_matmul(LA, LB, LC, M, N, K);
        end = clock();
        time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Time spent on outer matmul: %f seconds\n", time_spent);
        if (NAIVE){
            check_result(ref_C, LC, 1024, 1024);
        }
    }

    if (L1){
        printf("executing l1 matmul now with tile %d and inner tile %d ...\n", TILE, INNER_TILE);
        initialise_large_matrices(LA, LB, LC);
        start = clock();
        l1_tiled_matmul(LA, LB, LC, M, N, K, TILE, INNER_TILE);
        end = clock();
        time_spent = (double)(end - start) / CLOCKS_PER_SEC;
        printf("Time spent on l1 matmul: %f seconds\n", time_spent);
        if (NAIVE){
            check_result(ref_C, LC, 1024, 1024);
        }
    }






    free(LA);
    free(LB);
    free(LC);

}


// gcc -o main main.c -lm
// ./main
/*
perf stat -e L1-dcache-loads,L1-dcache-load-misses,l2_rqsts.references,l2_rqsts.miss ./main

matmul
executing l1 matmul now with tile 128 and inner tile 32 ...
Time spent on l1 matmul: 4.677510 seconds

 Performance counter stats for './main':

       21718147654      L1-dcache-loads:u
          18910847      L1-dcache-load-misses:u   #    0.09% of all L1-dcache accesses
          36205262      l2_rqsts.references:u
          14787220      l2_rqsts.miss:u

       4.815381895 seconds time elapsed

       4.805314000 seconds user
       0.009990000 seconds sys


*/