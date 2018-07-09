#include <stdio.h>
#include <x86intrin.h>
#include <time.h>
#include <assert.h>
#include <string.h>
#include "mul.h"
#include "syndrome.h"

// n - set length
// syndromes_to_calc - only odd ones (1, 3, 5, 7, ...)
uint64_t* syndromes(uint64_t set[], int n, int syndromes_to_calc) {
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

void measure_syn_time(int set_size, int tries) {
    int errors = set_size;
    uint64_t* set = malloc(set_size);
    
    int n = 100000;
    double total_times = 0;
    clock_t start, end;
    for (int try = 0; try < tries; try++) {
        start = clock() ;
        for (int i = 0; i < n; i++) {
            set = syndromes(set, set_size, errors);
            
        }
        end = clock();
        double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
        total_times += n / elapsed_time;
    }
    printf("Average syndromes times per 1 second: %llu\n", (uint64_t)(total_times / tries));
}


