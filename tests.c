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
    uint64_t* res = malloc(sizeof(uint64_t) * n);
    for (int i = 0; i < n; i++)
        res[i] = rand() * rand();
    return res;
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

void measure_calc_syn_time(int set_size, int max_errors, int tries) {
    assert(max_errors < set_size * 2);
    uint64_t* set = generate_random_set(set_size);
    
    double total_times = 0;
    double elapsed_time = 0;
    clock_t start, end;
    for (int try = 0; try < tries; try++) {
        start = clock() ;
        set = find_odd_syndromes(set, set_size, max_errors / 2);
        end = clock();
        elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
        total_times += elapsed_time;
        set = generate_random_set(set_size);
    }
    printf("Average execution time for calc syndromes (sender side): %f\n", total_times / tries);
}

void measure_decode_syn_time(int set_size, int max_errors, int errors, int tries) {
    int syndromes = max_errors * 2;
    assert(errors <= max_errors);
    assert(errors % 2 == 0);
    
    uint64_t sender_syndromes[syndromes];
    uint64_t receiver_syndromes[syndromes];
    static uint64_t error_loc_poly[MAX_DEGREE];
    
    int diffs_found = 0;
    double total_times = 0;
    clock_t start, after1, after2, end;
    double elapsed_time = 0;
    
    double total_time_reconstruct = 0;
    double total_time_decode = 0;
    double total_time_find_diff = 0;
    
    for (int try = 0; try < tries; try++) {
        uint64_t* receiver_set = generate_random_set(set_size);
        uint64_t* sender_set = generate_random_set(set_size);
        memcpy(receiver_set, sender_set, (set_size - errors / 2) * sizeof(uint64_t));
        uint64_t* sender_odd_syndromes = find_odd_syndromes(sender_set, set_size, syndromes / 2);
        uint64_t* receiver_odd_syndromes = find_odd_syndromes(receiver_set, set_size, syndromes / 2);

        reconstruct_all_syndromes(receiver_odd_syndromes, syndromes / 2, receiver_syndromes);
        start = clock();
        reconstruct_all_syndromes(sender_odd_syndromes, syndromes / 2, sender_syndromes);
        
        after1 = clock();
        total_time_reconstruct += (after1-start)/(double)CLOCKS_PER_SEC;

        uint64_t* diff_syndromes = xor_sets(receiver_syndromes, sender_syndromes, syndromes);
        
        decode_syndromesBM(diff_syndromes, syndromes, error_loc_poly, max_errors);
        after2 = clock();
        total_time_decode += (after2-after1)/(double)CLOCKS_PER_SEC;
        uint64_t* diff = find_diff(error_loc_poly, max_errors + 1, receiver_set, set_size, &diffs_found);
        end = clock();
        total_time_find_diff += (end-after2)/(double)CLOCKS_PER_SEC;
        total_times += elapsed_time;
        
        assert(error_loc_poly);
        assert(diff);
        assert(diffs_found * 2 == errors);
        for (int i = 0; i < diffs_found; i++) {
            assert(diff[i] == receiver_set[set_size - (errors / 2) + i]);
        }
        
    }
    printf("Average execution time for reconstructing syndromes (receiver side): %f\n", total_time_reconstruct / tries);
    printf("Average execution time for decoding syndromes (receiver side): %f\n", total_time_decode / tries);
    printf("Average execution time for finding diff syndromes (receiver side): %f\n", total_time_find_diff / tries);
}

void test_full_bch() {
    int tries = 100;
    int set_size = 100;
    int errors = set_size;
    int syndromes = errors * 2;

    uint64_t* odd_syndromes;
    uint64_t all_syndromes[syndromes];
    uint64_t error_loc_poly[MAX_DEGREE];
    uint64_t error_loc_poly2[MAX_DEGREE];
    
    for (int t = 0; t < tries; t++) {
        uint64_t* set = generate_random_set(set_size);
        odd_syndromes = find_odd_syndromes(set, set_size, syndromes / 2);
        reconstruct_all_syndromes(odd_syndromes, syndromes / 2, all_syndromes);
        decode_syndromesPGZ(all_syndromes, syndromes, error_loc_poly);
        decode_syndromesBM(all_syndromes, syndromes, error_loc_poly2, errors);
        
        for (int i = 0; i < errors + 1; i++) {
            assert(error_loc_poly2[i] == error_loc_poly[i]);
        }
        
        for (int j = 0; j < set_size; j++) {
            uint64_t res = eval_in_poly(error_loc_poly, errors + 1, set[j]);
            assert(res == 0);
        }
    }
}
