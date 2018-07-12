#include <stdio.h>
#include <x86intrin.h>
#include <time.h>

#include "mul.h"

uint64_t mul(uint64_t a, uint64_t b) {
    const __m128i LOW_MODULO = _mm_set_epi64x(0, ((uint64_t)1 << 4) | (1 << 3) | (1 << 1) | 1);

    __m128i a1 = _mm_cvtsi64_si128(a);
    __m128i b1 = _mm_cvtsi64_si128(b);
    
    __m128i product1 = _mm_clmulepi64_si128(a1, b1, 0x00);
    __m128i result = product1;
    
    __m128i product2 = _mm_clmulepi64_si128(LOW_MODULO, product1, 0xF2);
    result = _mm_xor_si128(result, product2);
    __m128i product3 = _mm_clmulepi64_si128(LOW_MODULO, product2, 0xF2);
    
    result = _mm_xor_si128(result, product3);
    return _mm_cvtsi128_si64(result);
}

void measure_mul_time(int tries) {
    uint64_t a = (uint64_t)(1) << 33;
    uint64_t b = (uint64_t)(1) << 33;
    
    int n = 1000000;
    double total_times = 0;
    clock_t start, end;
    for (int try = 0; try < tries; try++) {
        start = clock() ;
        for (int i = 0; i < n; i++) {
            a = mul(a, b);
        }
        end = clock();
        double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
        total_times += n / elapsed_time;
    }
    printf("Average mul times per 1 second: %llu\n", (uint64_t)(total_times / tries));
}
