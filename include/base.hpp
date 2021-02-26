#ifndef MM_BASIC_HPP
#define MM_BASIC_HPP

#include <cstddef>
#include <cmath>
#include <type_traits>

namespace base {
template<
    typename T, 
    typename = std::enable_if_t<std::is_arithmetic_v<T>>
>
void multiply(size_t m, size_t n, size_t l, 
              T** A, T** B, T** C) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C[i][j] = 0;
            for (size_t k = 0; k < l; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// Uses a cache line size.
template<
    typename T, 
    size_t CacheSize = 64,
    typename = std::enable_if_t<std::is_arithmetic_v<T>>
>
void multiply_cache_friendly(size_t m, size_t n, size_t l, 
                             T** A, T** B, T** C) {
    auto t = static_cast<size_t>(std::sqrt(CacheSize));

    for (size_t i = 0; i < m; i += t) {
        auto ip = std::min(i + t, m);
        for (size_t j = 0; j < n; j += t) {
            auto jq = std::min(j + t, n);
            C[i][j] = 0;
            for (size_t k = 0; k < l; k += t) {
                auto kr = std::min(k + t, l);
                for (size_t p = i; p < ip; ++p) {
                    for (size_t q = j; q < jq; ++q) {
                        for (size_t r = k; r < kr; ++r) {
                            C[p][q] += A[p][r] * B[r][q];
                        }
                    }
                }
            }
        }
    }
}
}

#endif
