#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include  "CustardFlow.h"

#define M 12000
#define N 256
#define K 784
#define TILE 128
#define INNER_TILE 32
#define NAIVE 1
#define OUTER 0
#define TILED 0
#define L1 0
#define ONLY_LARGE 0



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
        naive_matmul(A, B, C, 4, 4, 4, 4, 4, 4);
    
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
        naive_matmul(LA, LB, ref_C, M, N, K, M, K, N);
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