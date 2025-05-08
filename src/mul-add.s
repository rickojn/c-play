    .text
    .global simd_mul_add
    .type   simd_mul_add, @function
simd_mul_add:
    # System V AMD64: pointer to a[] in %rdi, pointer to b[] in %rsi

    vmovups  (%rdi), %ymm0         # load 8 floats from a into ymm0
    vmovups  (%rsi), %ymm1         # load 8 floats from b into ymm1
    vfmadd231ps %ymm1, %ymm0, %ymm0 # ymm0 = ymm0 * ymm1 + ymm0
    vmovups  %ymm0, (%rdi)         # store result back into a
    vzeroupper                     # avoid AVX→SSE transition penalty
    ret
