#include "inverse.h"
#include <string.h>
#include "mul.h"

uint64_t inverse(uint64_t a) {
    uint64_t res = a;
    // 5 iterations: 3, 7, 15, 31, 62
    for (int i = 0; i < 5; i++) {
        uint64_t old_res = res;
        for (int j = 0; j < i + 1; j++) {
            res = mul(res, res);
        }
        res = mul(res, old_res);
    }
    return mul(res, a);
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
