#include <stdio.h>
#include <x86intrin.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include "mul.h"
#include "inverse.h"
#include "syndrome.h"

// n — set length
uint64_t* find_odd_syndromes(uint64_t set[], int n, int syndromes_to_calc) {
    static uint64_t matr[2][MAX_DEGREE];
    uint64_t* result = calloc(syndromes_to_calc, sizeof(uint64_t));
    
    // init matrix with power 1 and 2, calc first syndrome
    for (int j = 0; j < n; j++) {
        matr[0][j] = set[j];
        matr[1][j] = mul(set[j], set[j]);
        result[0] ^= matr[0][j];
    }
    
    for (int i = 1; i < syndromes_to_calc; i++)
    {
        for (int j = 0; j < n; j++) {
            matr[0][j] = mul(matr[0][j], matr[1][j]);
            result[i] ^= matr[0][j];
        }
    }
    return result;
}

void reconstruct_all_syndromes(uint64_t odd_syndromes[], int n, uint64_t all_syndromes[]) {
    for (int i = 0; i < n; i++) {
        all_syndromes[i*2] = odd_syndromes[i];
        all_syndromes[i*2+1] = mul(all_syndromes[i], all_syndromes[i]);
    }
}

uint64_t* xor_sets(uint64_t basic_set[], uint64_t add_set[], int n) {
    uint64_t* res = malloc(sizeof(uint64_t) * n);
    for (int i = 0; i < n; i++)
        res[i] = basic_set[i] ^ add_set[i];
    return res;
}


void decode_syndromesPGZ(uint64_t syndromes[], int n, uint64_t error_loc_poly[]) {
    assert(~(n&1));
    assert(n <= MAX_DEGREE);
    static uint64_t A[MAX_DEGREE][MAX_DEGREE];
    int low_dim = n / 2;
    for (int i = 0; i < low_dim; i++) {
        for (int j = 0; j < low_dim + 1; j++) {
            A[i][j] = syndromes[i+j];
        }
    }
    
    static uint64_t echelon_elements[MAX_DEGREE];
    // {i,i} — echelon
    for (int i = 0; i < low_dim; i++) {
        // rows we want to get 0 in
        for (int j = i + 1; j < low_dim + i; j++) {
            int current_row = j % low_dim;
            if (A[current_row][i] == 0)
                continue;
            for (int k = i + 1; k < low_dim + 1; k++) {
                A[current_row][k] = mul(A[i][i], A[current_row][k]) ^ mul(A[i][k], A[current_row][i]);
            }
            A[current_row][i] = 0;
        }
        echelon_elements[i] = A[i][i];
    }
    
    // multiply echelons
    for (int i = 1; i < low_dim; i++) {
        for (int j = 0; j < i; j++)
            echelon_elements[j] = mul(echelon_elements[j], echelon_elements[i]);
    }
    
    static uint64_t echelon_inversions[MAX_DEGREE];
    inverses(echelon_elements, low_dim, echelon_inversions);
    
    for (int i = 0; i < low_dim; i++) {
        error_loc_poly[i + 1] = mul(echelon_inversions[low_dim - i - 1], A[low_dim - i - 1][low_dim]);
    }
    error_loc_poly[0] = 1;
}

// Horner
uint64_t eval_in_poly(uint64_t poly[], int size, uint64_t x0) {
    uint64_t res = poly[0];
    for (int i = 1 ; i < size; i++) {
        res = poly[i] ^ mul(res, x0);
    }
    return res;
}


uint64_t* find_diff(uint64_t error_loc_poly[], int size1, uint64_t candidates[], int size2, int* diffs_found){
    uint64_t* res = malloc(size1 * sizeof(uint64_t));
    int count = 0;
    for (int i = 0; i < size2; i++) {
        if (eval_in_poly(error_loc_poly, size1, candidates[i]) == 0) {
            res[count] = candidates[i];
            count++;
        }
    }
    *diffs_found = count;
    return realloc(res, count * sizeof(uint64_t));
}



static inline void pzero(int degree, uint64_t dst[]) {
    for (int i = 0; i < degree; i++) {
        dst[i] = 0;
    }
}

static inline void pcopy(int degree, uint64_t dst[], uint64_t src[]) {
    for (int i = 0; i < degree; i++) {
        dst[i] = src[i];
    }
}

static inline void pshift(int degree, uint64_t dst[], const uint64_t src[]) {
    for (int i = degree - 1; i > 0; i--) {
        dst[i] = src[i - 1];
    }
    dst[0] = 0;
}


static void pmul(int degree, uint64_t dst[], const uint64_t x[], const uint64_t y[]) {
    for (int i = 0; i < degree; i++) {
        dst[i] = 0;
    }
    for (int shift = 0; shift < degree; shift++) {
        for (int i = 0; i + shift < degree; i++) {
            dst[i + shift] ^= mul(x[i], y[shift]);
        }
    }
}

static void init_gamma(int degree, uint64_t gamma[], int codelen) {
    pzero(degree, gamma);
    gamma[0] = 1;
}

static inline uint64_t compute_discrepancy(const uint64_t lambda[], const uint64_t S[], int L, int n) {
    uint64_t sum = 0;
    for (int i = 0; i <= L; i++) {
        sum ^= mul(lambda[i], S[n - i]);
    }
    return sum;
}

void decode_syndromesBM(uint64_t syndromes[], int N, uint64_t lambda[], int codelen) {
    static uint64_t D[MAX_DEGREE];
    init_gamma(2 * N, lambda, codelen);
    pshift(2 * N, D, lambda);
    int k = -1, L = 0;
    for (int n = 0; n < N; n++) {
        uint64_t d = compute_discrepancy(lambda, syndromes, L, n);
        if (d != 0) {
            static uint64_t lambda2[MAX_DEGREE];
            for (int i = 0; i < 2 * N; i++) {
                lambda2[i] = lambda[i] ^ mul(D[i], d);
            }
            if (L < n - k) {
                int L2 = n - k;
                k = n - L;
                uint64_t d_inv = inverse(d);
                for (int i = 0; i < 2 * N; i++) {
                    D[i] = mul(lambda[i], d_inv);
                }
                L = L2;
            }
            pcopy(2 * N, lambda, lambda2);
        }
        pshift(2 * N, D, D);
    }
}
