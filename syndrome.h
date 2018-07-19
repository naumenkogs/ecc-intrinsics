#include <inttypes.h>

uint64_t* find_odd_syndromes(uint64_t set[], int n, int syndromes_to_calc);
void reconstruct_all_syndromes(uint64_t* odd_syndromes, int n, uint64_t* all_syndromes);
uint64_t* xor_sets(uint64_t basic_set[], uint64_t add_set[], int n);

// Peterson–Gorenstein–Zierler algorithm (slower)
void decode_syndromesPGZ(uint64_t* syndromes, int n, uint64_t* error_loc_poly);
// Berlekamp-Massey (faster)
void decode_syndromesBM(uint64_t* syndromes, int n, uint64_t* error_loc_poly, int errors);

uint64_t eval_in_poly(uint64_t poly[], int size, uint64_t x0);
uint64_t* find_diff(uint64_t error_loc_poly[], int size1, uint64_t candidates[], int size2, int* diffs_found);
