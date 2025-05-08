#include <immintrin.h>
void simd_mul_add(float *a, float *b) {
    __m256 va = _mm256_loadu_ps(a);
    __m256 vb = _mm256_loadu_ps(b);
    va = _mm256_fmadd_ps(va, vb, va);
    _mm256_storeu_ps(a, va);
}
