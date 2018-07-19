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
    
    int transactions = atoi(argv[1]);
    int max_errors = atoi(argv[2]);
    int errors = atoi(argv[3]);
    measure_calc_syn_time(transactions, max_errors, 10);
    measure_decode_syn_time(transactions, max_errors, errors, 10);

    return 0;
}
