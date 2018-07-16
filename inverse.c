#include "inverse.h"
#include <string.h>
#include "mul.h"

uint64_t inverse(uint64_t a) {
    const int squares1[] = {1, 2, 4, 8, 16};
    const int squares2[] = {4, 8, 16, 32};
    uint64_t powers[5];
    uint64_t cur_element = a;
    // iterations: p2, p4, p8, p16, p32
    for (int i = 0; i < 5; i++) {
        uint64_t prev_element = cur_element;
        for (int j = 0; j < squares1[i]; j++) {
            cur_element = mul(cur_element, cur_element);
        }
        cur_element = mul(cur_element, prev_element);
        powers[i] = cur_element;
    }
    uint64_t res = mul(a, powers[0]);
    for (int i = 0; i < 4; i++) {
        res = mul(res, powers[i]);
        for (int j = 0; j < squares2[i]; j++) {
            res = mul(res, res);
        }
    }
    res = mul(res, powers[4]);

    return mul(res, res);
}

void inverses(uint64_t values[], int n, uint64_t inverses[]) {
    // products: [a, ab, abc, abcd, ...]
    uint64_t products[n];
    products[0] = values[0];
    for (int i = 1; i < n; i++)
        products[i] = mul(products[i-1], values[i]);
    uint64_t product_inverse = inverse(products[n - 1]);
    // helper product: 1, e, ed, edc, edcb, ...
    uint64_t helper_product = 1;
    for (int i = n - 1; i >= 1; i--) {
        uint64_t tmp_res = mul(product_inverse, products[i - 1]);
        inverses[i] = mul(tmp_res, helper_product);
        helper_product = mul(helper_product, values[i]);
    }
    inverses[0] = mul(helper_product, product_inverse);
}
