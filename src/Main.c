// main.c
#include <stdio.h>

// Declare the external function (same signature as in the header)
extern void simd_mul_add(float *a, float *b);

int main(void) {
    float a[8] = { 1, 2, 3, 4, 5, 6, 7, 8 };
    float b[8] = { 0.5f, 1.5f, -2.0f, 4.0f, 0.25f, -1.0f, 2.0f, 3.0f };

    simd_mul_add(a, b);

    printf("Result:\n");
    for (int i = 0; i < 8; i++)
        printf("  a[%d] = %f\n", i, a[i]);
    return 0;
}
// This program demonstrates the use of SIMD (Single Instruction, Multiple Data)
// operations to perform a multiply-add operation on two arrays of floats.
// The function simd_mul_add multiplies each element of array a by the corresponding
// element of array b and adds the result to the corresponding element of array a.
// The result is stored back in array a. The program initializes two arrays, a and b,
// calls the simd_mul_add function, and prints the resulting array a.
