#include "tests.h"
#include "mul.h"
#include "syndrome.h"
#include "inverse.h"

#include <x86intrin.h>
#include <assert.h>
#include <time.h>
#include <stdio.h>
#include "minisketch.h"
#include <vector>


void TestInvMul() {
    int tries = 100;
    for (int i = 0; i < tries; i++) {
        uint64_t poly = rand();
        uint64_t inv_poly = inverse(poly);
        uint64_t product = mul(inv_poly, poly);
        assert(product == 1);
    }
}

void TestInverses() {
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

void TestSquare() {
    int tries = 100;
    for (int i = 0; i < tries; i++) {
        uint64_t val = rand() * rand();
        assert(square(val) == mul(val, val));
    }
}

uint64_t* GenerateRandomSet(int n) {
    uint64_t* res = new uint64_t[n];
    for (int i = 0; i < n; i++)
        res[i] = rand() * rand();
    return res;
}

void MeasureMulTime(int tries) {
    uint64_t a = (uint64_t)(1) << 33;
    uint64_t b = (uint64_t)(1) << 33;
    
    int n = 1000000;
    double total_times = 0;
    clock_t start, end;
    for (int tried = 0; tried < tries; tried++) {
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

void MeasureSquareTime(int tries) {
    uint64_t a = (uint64_t)(1) << 33;
    
    int n = 1000000;
    double total_times = 0;
    clock_t start, end;
    for (int tried = 0; tried < tries; tried++) {
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

void MeasureCalcSynTime(int set_size, int max_errors, int tries) {
    assert(max_errors <= set_size * 2);
    uint64_t* set = GenerateRandomSet(set_size);
    uint64_t odd_syndromes[max_errors / 2];
    
    double total_times = 0;
    double elapsed_time = 0;
    clock_t start, end;
    for (int tried = 0; tried < tries; tried++) {
        start = clock() ;
        find_odd_syndromes(set, set_size, max_errors / 2, odd_syndromes);
        end = clock();
        elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
        total_times += elapsed_time;
        set = GenerateRandomSet(set_size);
    }
    printf("Average execution time for calc syndromes (sender side): %f\n", total_times / tries);
}

void MeasureDecodeSynTime(int set_size, int max_errors, int errors, int tries) {
    int syndromes = max_errors * 2;
    assert(errors <= max_errors);
    assert(errors % 2 == 0);
    
    static uint64_t error_loc_poly[MAX_DEGREE];
    
    int diffs_found = 0;
    clock_t start, end;
    double elapsed_time = 0;
    
    double total_time_reconstruct = 0;
    double total_time_decode = 0;
    double total_time_find_diff = 0;
    
    for (int tried = 0; tried < tries; tried++) {
        uint64_t* receiver_set = GenerateRandomSet(set_size);
        uint64_t* sender_set = GenerateRandomSet(set_size);
        uint64_t* sender_odd_syndromes = new uint64_t[syndromes / 2]();
        uint64_t* receiver_odd_syndromes = new uint64_t[syndromes / 2]();
        uint64_t* diff_odd_syndromes = new uint64_t[syndromes / 2]();
        uint64_t* diff_all_syndromes = new uint64_t[syndromes]();

        for (int i = 0; i < set_size - errors / 2; i++)
            receiver_set[i] = sender_set[i];

        find_odd_syndromes(sender_set, set_size, syndromes / 2, sender_odd_syndromes);
        find_odd_syndromes(receiver_set, set_size, syndromes / 2, receiver_odd_syndromes);

        start = clock();
        xor_sets(receiver_odd_syndromes, sender_odd_syndromes, syndromes / 2, diff_odd_syndromes);
        reconstruct_all_syndromes(diff_odd_syndromes, syndromes / 2, diff_all_syndromes);
        end = clock();
        total_time_reconstruct += (end-start)/(double)CLOCKS_PER_SEC;

        start = clock();
        decode_syndromesBM(diff_all_syndromes, syndromes, error_loc_poly, max_errors);
        end = clock();
        total_time_decode += (end - start)/(double)CLOCKS_PER_SEC;

//        uint64_t* diff = find_diff(error_loc_poly, max_errors + 1, receiver_set, set_size, &diffs_found);
//        end = clock();
//        total_time_find_diff += (end-after2)/(double)CLOCKS_PER_SEC;
//
//        assert(diff);
//        assert(diffs_found * 2 == errors);
//        for (int i = 0; i < diffs_found; i++) {
//            assert(diff[i] == receiver_set[set_size - (errors / 2) + i]);
//        }
        std::vector<uint64_t> poly(error_loc_poly, error_loc_poly + errors + 1);
        std::reverse(poly.begin(), poly.end());

        start = clock();
        std::vector<uint64_t> roots;
        FindRoots(poly, roots);
        end = clock();
        total_time_find_diff += (end-start)/(double)CLOCKS_PER_SEC;

        for (int i = set_size - errors / 2; i < set_size; i++)
            assert(std::find(roots.begin(), roots.end(), receiver_set[i]) != roots.end());

        
    }
    printf("Average execution time for reconstructing diff syndromes (receiver side): %f\n", total_time_reconstruct / tries);
    printf("Average execution time for decoding syndromes (receiver side): %f\n", total_time_decode / tries);
    printf("Average execution time for finding diff syndromes (receiver side): %f\n", total_time_find_diff / tries);
}

void TestFullBCH() {
    int tries = 10;
    int set_size = 100;
    int errors = set_size;
    int syndromes = errors * 2;

    static uint64_t error_loc_poly[MAX_DEGREE];
    static uint64_t error_loc_poly2[MAX_DEGREE];
    for (int t = 0; t < tries; t++) {
        uint64_t all_syndromes[syndromes];
        uint64_t* odd_syndromes = new uint64_t[syndromes / 2]();

        uint64_t* set = GenerateRandomSet(set_size);
        set_size = 2;
        set[0] = 1;
        set[1] = 2;
        errors = set_size;
        syndromes = errors * 2;

        find_odd_syndromes(set, set_size, syndromes / 2, odd_syndromes);
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

void TestFindRoots() {
    int tries = 1;
    int set_size = 1000;
    int errors = set_size;
    int syndromes = errors * 2;

    
    uint64_t all_syndromes[syndromes];
    uint64_t error_loc_poly[MAX_DEGREE];
    uint64_t error_loc_poly2[MAX_DEGREE];
    
    
    for (int t = 0; t < tries; t++) {
        uint64_t* odd_syndromes = new uint64_t[syndromes / 2]();
        uint64_t* set = GenerateRandomSet(set_size);

    
        set_size = 2;
        set[0] = 1;
        set[1] = 2;
        errors = set_size;
        syndromes = errors*2;
        
    
        find_odd_syndromes(set, set_size, syndromes / 2, odd_syndromes);
        reconstruct_all_syndromes(odd_syndromes, syndromes / 2, all_syndromes);
        decode_syndromesBM(all_syndromes, syndromes, error_loc_poly, errors);


        std::vector<uint64_t> poly;

        for (int i = errors; i >= 0; i--)
            poly.push_back(error_loc_poly[i]);

        std::vector<uint64_t> roots;
        FindRoots(poly, roots);
        for (int i = 0; i < set_size; i++)
            assert(std::find(roots.begin(), roots.end(), set[i]) != roots.end());
        
    }
}


