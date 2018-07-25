#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <time.h>

#include <vector>
#include "mul.h"
#include "inverse.h"
#include "minisketch.h"

#include <iostream>
#include <random>
#include <limits>


/** Square a field element. */
static uint64_t sqr(uint64_t x) {
    return mul(x, x);
}


/** Compute the remainder of a polynomial division of val by mod, putting the result in mod. */
void Mod(const std::vector<uint64_t>& mod, std::vector<uint64_t>& val) {
    size_t modsize = mod.size();
    assert(modsize > 0 && mod.back() == 1);
    if (val.size() < modsize) return;
    assert(val.back() != 0);
    while (val.size() >= modsize) {
        uint64_t term = val.back();
        val.pop_back();
        if (term) {
            for (size_t x = 0; x < mod.size() - 1; ++x) {
                val[val.size() - modsize + 1 + x] ^= mul(term, mod[x]);
            }
        }
    }
    while (val.size() > 0 && val.back() == 0) val.pop_back();
}

/** Compute the quotient of a polynomial division of val by mod, putting the quotient in div and the remainder in val. */
void DivMod(const std::vector<uint64_t>& mod, std::vector<uint64_t>& val, std::vector<uint64_t>& div) {
    size_t modsize = mod.size();
    assert(mod.size() > 0 && mod.back() == 1);
    if (val.size() < mod.size()) {
        div.clear();
        return;
    }
    assert(val.back() != 0);
    div.resize(val.size() - mod.size() + 1);
    while (val.size() >= modsize) {
        uint64_t term = val.back();
        div[val.size() - modsize] = term;
        val.pop_back();
        if (term) {
            for (size_t x = 0; x < mod.size() - 1; ++x) {
                val[val.size() - modsize + 1 + x] ^= mul(term, mod[x]);
            }
        }
    }
}

/** Make a polynomial monic. */
void MakeMonic(std::vector<uint64_t>& a) {
    assert(a.back() != 0);
    if (a.back() == 1) return;
    uint64_t fac = inverse(a.back());
    a.back() = 1;
    for (size_t i = 0; i < a.size() - 1; ++i) {
        a[i] = mul(a[i], fac);
    }
}


/** Compute the GCD of two polynomials, putting the result in a. b will be cleared. */
void GCD(std::vector<uint64_t>& a, std::vector<uint64_t>& b) {
    if (a.size() < b.size()) std::swap(a, b);
    while (b.size() > 0) {
        if (b.size() == 1) {
            a.resize(1);
            a[0] = 1;
            return;
        }
        MakeMonic(b);
        Mod(b, a);
        std::swap(a, b);
    }
}

static void Derivative(const std::vector<uint64_t>& a, std::vector<uint64_t>& res) {
    int even_degree = a.size() % 2;
    for (int i = even_degree; i < a.size() - 1; i+=2) {
        res.push_back(a[i]);
        res.push_back(0);
    }
    res.pop_back();
}


/** Square a polynomial. */
void Sqr(std::vector<uint64_t>& poly) {
    if (poly.size() == 0) return;
    poly.resize(poly.size() * 2 - 1);
    for (int x = poly.size() - 1; x >= 0; --x) {
        poly[x] = (x & 1) ? 0 : sqr(poly[x / 2]);
    }
}

/** Compute the trace map of (param*x) modulo mod, putting the result in out. */
void TraceMod(const std::vector<uint64_t>& mod, std::vector<uint64_t>& out, uint64_t param) {
    out.resize(2);
    out[0] = 0;
    out[1] = param;
    out.reserve(2 * mod.size());
    
    for (int i = 0; i < 63; ++i) {
        Sqr(out);
        if (out.size() < 2) out.resize(2);
        out[1] = param;
        Mod(mod, out);
    }
}


/** One step of the root finding algorithm; finds roots of poly and puts them in roots. */
bool RecFindRoots(std::vector<uint64_t>& poly, std::vector<uint64_t>& roots, bool known_distinct) {
    assert(poly.size() > 0 && poly.back() == 1);
    if (poly.size() == 1) return true;
    if (poly.size() == 2) {
        roots.push_back(poly[0]);
        return true;
    }

    std::vector<uint64_t> trace;
    std::vector<uint64_t> tmp;
    for (int iter = 0;; ++iter) {
        std::random_device rd;
        std::mt19937_64 eng(rd());
        std::uniform_int_distribution<unsigned long long> distr;

        
        TraceMod(poly, trace, distr(eng));

        if (iter == 1 && !known_distinct) {
            // Only check for distinct roots after a failed iteration
            tmp = trace;
            Sqr(tmp);
            for (size_t i = 0; i < trace.size(); ++i) {
                tmp[i] ^= trace[i];
            }
            while (tmp.size() && tmp.back() == 0) tmp.pop_back();
            Mod(poly, tmp);
            if (tmp.size() != 0) return false;
            known_distinct = true;
        }
        tmp = poly;
        GCD(trace, tmp);
        if (trace.size() != poly.size() && trace.size() > 1) {
            MakeMonic(trace);
            DivMod(trace, poly, tmp);
            if (!RecFindRoots(trace, roots, known_distinct)) return false;
            if (!RecFindRoots(tmp, roots, known_distinct)) return false;
            break;
        }
    }
    return true;
}

bool IsSquareFree(const std::vector<uint64_t>& poly) {
    std::vector<uint64_t> der;
    Derivative(poly, der);
    auto copy = poly;
    GCD(copy, der);
    return copy.size() == 1 && copy[0] == 1;
}

/** Find roots of poly and put them in roots. Poly must be square free and only have 1st degree factors. */
void FindRoots(const std::vector<uint64_t>& poly, std::vector<uint64_t>& roots) {
    if (!IsSquareFree(poly)) {
        printf("Not square free \n");
        return;
    }
    roots.clear();
    roots.reserve(poly.size());
    auto copy = poly;
    RecFindRoots(copy, roots, false);
}
