#ifndef MM_UTILS_HPP
#define MM_UTILS_HPP

#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <type_traits>

namespace utils {
// Returns a 2D memory-contiguous array which supports the standard subscript operator.
template<typename T = double>
T** allocate(size_t m, size_t n) {
    T** array = new T* [m];
    *array = new T [m * n];

    for (size_t i = 1; i < m; ++i) {
        *(array + i) = *(array + i - 1) + n;
    }

    return array;
}

template<typename T = double>
void free(T** array) {
    delete[] *array;
    delete[] array;
}

template<typename T = double>
void fill(size_t m, size_t n, T** A, size_t bound) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            A[i][j] = static_cast<T>(rand() % bound);
        }
    }
}

template<typename T = double>
void print(size_t m, size_t n, T** A) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            printf("%f\t", A[i][j]);
        }
        printf("\n");
    }
}

template<typename T = double>
bool compare(size_t m, size_t n, T** A, T** B) {
    auto result = true;
    for (size_t i = 0; i < m; ++i) {
        auto c = std::memcmp(A[i], B[i], n);
        result &= c == 0;
    }
    return result;
}

template<typename T>
T min(T x, T y) {
    return x < y ? x : y;
}

template<typename T, typename... Ts>
T min(T x, T y, Ts... xs) {
    return min(min(x, y), xs...);
}
}

#endif
