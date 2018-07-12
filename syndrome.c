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
    uint64_t matr[2][n];
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
        all_syndromes[i*2+1] = mul(odd_syndromes[i], odd_syndromes[i]);
    }
}

void decode_syndromes(uint64_t syndromes[], int n, uint64_t lambdas[]) {
    assert(~(n&1));
    int low_dim = n / 2;
    uint64_t A[low_dim][low_dim + 1];
    for (int i = 0; i < low_dim; i++) {
        for (int j = 0; j < low_dim + 1; j++) {
            A[i][j] = syndromes[i+j];
        }
    }

    uint64_t echelon_elements[low_dim];
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

    uint64_t echelon_inversions[low_dim];
    inverses(echelon_elements, low_dim, echelon_inversions);

    for (int i = 0; i < low_dim; i++) {
        lambdas[i] = mul(echelon_inversions[i], A[i][low_dim]);
    }
}


void measure_syn_time(int set_size, int tries) {
    int errors = set_size;
    uint64_t* set = malloc(set_size);
    
    int n = 100000;
    double total_times = 0;
    clock_t start, end;
    for (int try = 0; try < tries; try++) {
        start = clock() ;
        for (int i = 0; i < n; i++) {
            set = find_odd_syndromes(set, set_size, errors);
            
        }
        end = clock();
        double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
        total_times += n / elapsed_time;
    }
    printf("Average syndromes times per 1 second: %llu\n", (uint64_t)(total_times / tries));
}
