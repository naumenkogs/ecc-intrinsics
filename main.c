#include <stdlib.h>
#include <time.h>
#include <stdio.h>

#include "mul.h"
#include "syndrome.h"
#include "inverse.h"
#include "tests.h"

int main(int argc, const char * argv[]) {
    srand(time(NULL));
    test_full_bch();
    
    int transactions = 1000;
//    measure_mul_time(100);
//    measure_square_time(100);
    measure_calc_syn_time(transactions, 10);
    measure_decode_syn_time(transactions, 10);


    return 0;
}
