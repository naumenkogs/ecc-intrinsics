#include <x86intrin.h>
#include "mul.h"

uint64_t reduce(__m128i a) {
    const __m128i LOW_MODULO = _mm_set_epi64x(0, ((uint64_t)1 << 4) | (1 << 3) | (1 << 1) | 1);
    __m128i result = a;
    
    __m128i product2 = _mm_clmulepi64_si128(LOW_MODULO, a, 0xF2);
    result = _mm_xor_si128(result, product2);
    __m128i product3 = _mm_clmulepi64_si128(LOW_MODULO, product2, 0xF2);

    result = _mm_xor_si128(result, product3);
    return _mm_cvtsi128_si64(result);
}

uint64_t mul(uint64_t a, uint64_t b) {
    __m128i a1 = _mm_cvtsi64_si128(a);
    __m128i b1 = _mm_cvtsi64_si128(b);
    
    __m128i product = _mm_clmulepi64_si128(a1, b1, 0x00);
    return reduce(product);
}

uint64_t square(uint64_t a) {
    const uint64_t SQUARE_MASK = 0x0AAAAAAAAAAAAAAAAULL >> 1;
    uint64_t low_a = _pdep_u64(a, SQUARE_MASK);
    uint64_t high_a = _pdep_u64(a >> 32, SQUARE_MASK);
    __m128i a1 = _mm_set_epi64x(high_a, low_a);
    return reduce(a1);
}
