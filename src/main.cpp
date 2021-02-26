#include "base.hpp"
#include "utils.hpp"
#include "strassen.hpp"

#include <cassert>
#include <chrono>
#include <cstdio>

using namespace std::chrono;
namespace ut = utils;

using Time = steady_clock::time_point;

int main(int argc, char** argv) {
    // A = m × l.
    // B = l × n.
    // C = m × n.
    size_t m = 256, n = 256, l = 256;

    if (argc > 3) {
        m = std::atoi(argv[1]);
        n = std::atoi(argv[2]);
        l = std::atoi(argv[3]);
    }

    auto A = ut::allocate(m, l);
    auto B = ut::allocate(l, n);

    ut::fill(m, l, A, 8);
    ut::fill(l, n, B, 8);

    auto C = ut::allocate(m, n);
    auto D = ut::allocate(m, n);
    auto E = ut::allocate(m, n);

    auto b1 = steady_clock::now();

    base::multiply(m, n, l, A, B, C);

    auto b2 = steady_clock::now();

    base::multiply_cache_friendly(m, n, l, A, B, D);

    auto b3 = steady_clock::now();

    strassen::multiply(m, n, l, A, B, E);

    auto end = steady_clock::now();

    assert(ut::compare(m, n, C, D));
    assert(ut::compare(m, n, D, E));

    auto delta = [](const Time& x, const Time& y) {
        return duration_cast<microseconds>(x - y).count();
    };

    printf("Base time (µs): %li\n", delta(b2, b1));
    printf("Cache friendly time (µs): %li\n", delta(b3, b2));
    printf("Strassen time (µs): %li\n", delta(end, b3));
    
    ut::free(A);
    ut::free(B);
    ut::free(C);
    ut::free(D);
    ut::free(E);

    return 0;
}
