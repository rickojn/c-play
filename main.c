#include <stdio.h>
#include <stdlib.h>
#include <time.h>



void matmul(const float* A, const float *B, float * C, size_t rows_C, size_t cols_C, size_t rows_B_cols_A){
    for (size_t i = 0; i < rows_C; i++ ){
        for (size_t j = 0; j < cols_C; j++){
            for (size_t k = 0; k < rows_B_cols_A; k++){            
                C[i * rows_C + j] += A[i * rows_C + k] * B[j *cols_C + k]; 
            }
        }
    }
}


void tiled_matmul(const float* A, const float *B, float * C, size_t rows_C, size_t cols_C, size_t rows_B_cols_A, size_t tile_size){
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

    float A[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    float B[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    float C[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    matmul(&A, &B, &C, 4, 4, 4);

    for (size_t i = 0; i < 4; i++)
    {
        printf("\n");
        for (size_t j = 0; j < 4; j++){
            printf("%f\t", C[i * 4 + j]);
        }
        printf("\n");
    }

    printf("\n");
    printf("tiled_matmul\n");
    
    float AT[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    float BT[] = {1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4};
    float CT[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    
    tiled_matmul(&AT, &BT, &CT, 4, 4, 4, 2);

    for (size_t i = 0; i < 4; i++)
    {
        printf("\n");
        for (size_t j = 0; j < 4; j++){
            printf("%f\t", C[i * 4 + j]);
        }
        printf("\n");
    }
}



// gcc -o main main.c -lm
// ./main