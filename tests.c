#include "tests.h"
#include "mul.h"
#include "syndrome.h"
#include "inverse.h"

#include <x86intrin.h>
#include <assert.h>
#include <time.h>
#include <stdio.h>

void test_inv_mul() {
    int tries = 100;
    for (int i = 0; i < tries; i++) {
        uint64_t poly = rand();
        uint64_t inv_poly = inverse(poly);
        uint64_t product = mul(inv_poly, poly);
        assert(product == 1);
    }
}

void test_inverses() {
    int tries = 100;
    int n = 5;
    uint64_t polys[n];
    uint64_t inversed_polys[n];
    for (int i = 0; i < tries; i++) {
        for (int j = 0; j < n; j++) {
            polys[j] = rand();
        }
        inverses(polys, n, inversed_polys);
        for (int j = 0; j < n; j++) {
            uint64_t product = mul(polys[j], inversed_polys[j]);
            assert(product == 1);
        }
    }
}

void test_square() {
    int tries = 100;
    for (int i = 0; i < tries; i++) {
        uint64_t val = rand() * rand();
        assert(square(val) == mul(val, val));
    }
}

uint64_t* generate_random_set(int n) {
    uint64_t* res = malloc(n * sizeof(uint64_t));
    for (int i = 0; i < n; i++)
        res[i] = rand() * rand();
    return res;
}

void test_full_bch() {
    int tries = 100;
    int set_size = 100;
    int syndromes = set_size * 2;
    
    uint64_t* odd_syndromes;
    uint64_t all_syndromes[syndromes];
    uint64_t error_loc_poly[syndromes / 2 + 1];
    
    for (int t = 0; t < tries; t++) {
        uint64_t* set = generate_random_set(set_size);
        odd_syndromes = find_odd_syndromes(set, set_size, syndromes / 2);
        reconstruct_all_syndromes(odd_syndromes, syndromes / 2, all_syndromes);
        decode_syndromes(all_syndromes, syndromes, error_loc_poly);
        
        for (int j = 0; j < set_size; j++) {
            uint64_t res = eval_in_poly(error_loc_poly, syndromes / 2 + 1, set[j]);
            assert(res == 0);
        }
    }
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

void measure_square_time(int tries) {
    uint64_t a = (uint64_t)(1) << 33;
    
    int n = 1000000;
    double total_times = 0;
    clock_t start, end;
    for (int try = 0; try < tries; try++) {
        start = clock() ;
        for (int i = 0; i < n; i++) {
            a = square(a);
        }
        end = clock();
        double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
        total_times += n / elapsed_time;
    }
    printf("Average squaring times per 1 second: %llu\n", (uint64_t)(total_times / tries));
}

void measure_calc_syn_time(int set_size, int tries) {
    int errors = set_size;
    uint64_t* set = generate_random_set(set_size);
    
    double total_times = 0;
    double elapsed_time = 0;
    clock_t start, end;
    for (int try = 0; try < tries; try++) {
        start = clock() ;
        set = find_odd_syndromes(set, set_size, errors);
        end = clock();
        elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
        total_times += elapsed_time;
        set = generate_random_set(set_size);
    }
    printf("Average execution time for calc syndromes (sender side): %f\n", total_times / tries);
}

void measure_decode_syn_time(int set_size, int tries) {
    int errors = set_size;
    int n_syndromes = errors * 2;
    uint64_t syndromes[n_syndromes];
    uint64_t error_loc_poly[errors + 1];
    
    uint64_t* sender_set = generate_random_set(set_size);
    uint64_t* client_set = generate_random_set(set_size);
    uint64_t* odd_syndromes = find_odd_syndromes(sender_set, set_size, errors);
    
    int diffs_found;
    double total_times = 0;
    clock_t start, end;
    double elapsed_time;
    
    for (int try = 0; try < tries; try++) {
        start = clock() ;
        
        reconstruct_all_syndromes(odd_syndromes, n_syndromes / 2, syndromes);
        decode_syndromes(syndromes, n_syndromes / 2, error_loc_poly);
        uint64_t* diff = find_diff(error_loc_poly, errors + 1, client_set, set_size, &diffs_found);
        end = clock();
        elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
        total_times += elapsed_time;
        
        assert(diff);
        client_set = generate_random_set(set_size);
        sender_set = generate_random_set(set_size);
        
    }
    printf("Average execution time for decoding syndromes (receiver side): %f\n", total_times / tries);
}
