#include <stdio.h>
#include "mul.h"
#include "syndrome.h"
#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>


int main(int argc, const char * argv[]) {
    assert(mul((uint64_t)(1) << 33, (uint64_t)(1) << 33) == 108);
    int set_size = 3;
    int syndromes = set_size * 2;
    uint64_t set[] = {0, 2, 3};
    uint64_t* odd_syndromes = find_odd_syndromes(set, set_size, syndromes / 2);
    
    uint64_t all_syndromes[syndromes];
    reconstruct_all_syndromes(odd_syndromes, syndromes / 2, all_syndromes);
    uint64_t lambdas[syndromes];
    decode_syndromes(all_syndromes, syndromes, lambdas);
    
    printf("Result: \n");
    for (int i = 0; i < syndromes; i++) {
        printf("%llu ", lambdas[i]);
    }
    printf("\n");
//    measure_mul_time(100);
//    measure_syn_time(10, 100);


    return 0;
}
