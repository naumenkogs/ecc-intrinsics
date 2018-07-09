#include <stdio.h>
#include "mul.h"
#include "syndrome.h"
#include <assert.h>
#include <inttypes.h>


int main(int argc, const char * argv[]) {
    assert(mul((uint64_t)(1) << 33, (uint64_t)(1) << 33) == 108);
    uint64_t set[] = {1, 2, 1};
    uint64_t* res = syndromes(set, 3, 4);
    assert(syndromes(set, 3, 4)[3] == 128);

    measure_mul_time(100);
    measure_syn_time(10, 100);


    return 0;
}
