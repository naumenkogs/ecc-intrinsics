#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include "mul.h"
#include "string.h"
#include <vector>

#include "syndrome.h"
#include "inverse.h"
#include "tests.h"
#include "assert.h"
#include "string.h"

#include "minisketch.h"

//void print_poly(uint64_t poly[], int size, std::string msg) {
//    assert(size >= 0);
//    printf("\n%s\n", msg);
//    for (int i = 0; i < size; i++)
//        printf("%llu, ", poly[i]);
//    printf("\n");
//}

int find_remainder(uint64_t* a, int size1, uint64_t* b, int size2, uint64_t* quotient) {
    // A is always bigger
    assert(size1 >= size2);
    uint64_t monic_b[size2];
    memcpy(monic_b, b, size1 * sizeof(uint64_t));
    uint64_t monic_coef = 1;
    if (b[0] != 1) {
        monic_coef = inverse(b[0]);
        monic_b[0] = 1;
        for (int i = 1; i < size2; i++) {
            monic_b[i] = mul(b[i], monic_coef);
        }
    }
    for (int i = 0; i < size1 - size2 + 1; i++) {
        quotient[i] = a[i];
        if (a[i] == 0)
            continue;
        for (int j = 1; j < size2; j++) {
            a[i+j] = mul(a[i], monic_b[j]) ^ a[i + j];
        }
    }

    int rem_starts_at = size1 - size2 + 1;
    int zeros = 0;
    for (int i = rem_starts_at; i < size1; i++) {
        if (a[i] == 0) {
            zeros++;
        } else {
            break;
        }
    }
    int rem_size = size2 - zeros - 1;
    if (monic_coef != 1) {
        for (int i = size1 - rem_size; i < size2; i++) {
            a[i] = mul(a[i], monic_coef);
        }
        for (int i = 0; i < size1 - size2 + 1; i++)
            quotient[i] = mul(quotient[i], monic_coef);
    }
    return rem_size;
}

uint64_t* gcd(uint64_t a[], int size1, uint64_t b[], int size2, int* gcd_size) {
    uint64_t* b1 = new uint64_t[size2];
    memcpy(b1, b, sizeof(uint64_t) * size2);
    uint64_t* a1 = new uint64_t[size1];
    memcpy(a1, a, sizeof(uint64_t) * size1);
    static uint64_t quotient[MAX_DEGREE];

    while (size2 > 1) {
        int old_size2 = size2;
        size2 = find_remainder(a1, size1, b1, size2, quotient);
        if (size2 == 0) {
            size2 = old_size2;
            break;
        } else {
            uint64_t* tmp_array = b1;
            b1 = &a1[size1 - size2];
            a1 = tmp_array;
            size1 = old_size2;
        }
    }
    *gcd_size = size2;
    return b1;
}


int main(int argc, const char * argv[]) {
    srand(time(NULL));
//    int size = 5;
//    uint64_t poly[5] = {1,1,1,1,0};
//    sff(poly, size);
//    uint64_t set[2] = {1,2};
//    int set_size = 2;
//    int syndromes = 4;
//    int errors = 2;
//    uint64_t all_syndromes[syndromes];
//    uint64_t error_loc_poly2[syndromes];
//    uint64_t* odd_syndromes = find_odd_syndromes(set, set_size, syndromes / 2);
//    reconstruct_all_syndromes(odd_syndromes, syndromes / 2, all_syndromes);
//    print_poly(all_syndromes, syndromes, "all syn");
//    decode_syndromesPGZ(all_syndromes, syndromes, error_loc_poly2);
//    print_poly(error_loc_poly2, errors+1, "poly");
//    sff(error_loc_poly2, errors + 1);

//    TestFindRoots();
    
    
    
    TestFullBCH();
    int transactions = atoi(argv[1]);
    int max_errors = atoi(argv[2]);
    int errors = atoi(argv[3]);
//    MeasureCalcSynTime(transactions, max_errors, 10);
//    MeasureDecodeSynTime(transactions, max_errors, errors, 10);

    return 0;
}
