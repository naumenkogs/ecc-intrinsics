#include <inttypes.h>

uint64_t* find_odd_syndromes(uint64_t set[], int n, int syndromes_to_calc);
void reconstruct_all_syndromes(uint64_t* odd_syndromes, int n, uint64_t* all_syndromes);
void decode_syndromes(uint64_t* syndromes, int n, uint64_t* lambdas);
void measure_syn_time(int set_size, int tries);
