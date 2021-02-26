#ifndef MM_STRASSEN_HPP
#define MM_STRASSEN_HPP

#include "base.hpp"
#include "utils.hpp"

#include <cassert>

#include "omp.h"

namespace strassen {
namespace detail {
namespace ut = utils;

template<typename T>
void add(size_t m, size_t n, T** A, T** B, T** C) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

template<typename T>
void sub(size_t m, size_t n, T** A, T** B, T** C) {
    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}

template<typename T>
void copy_block(size_t m, size_t ma, size_t na, T** A, T** B) {
    for (size_t i = 0; i < m; ++i) {
        A[i] = &B[ma + i][na];
    }
}

template<typename T, size_t N>
void recursive_multiply(size_t mi, size_t ni, size_t li,
                        T** A, T** B, T** C) {
    // auto s = min(ml, nl, ll);

    if (mi * ni * li < N) {
        base::multiply(mi, ni, li, A, B, C);
        return;
    }

    auto m = mi / 2;
    auto n = ni / 2;
    auto l = li / 2;

    // Store results of add/sub.
    auto __AP1 = ut::allocate<T>(m, l);
    auto __BP1 = ut::allocate<T>(l, n);
    auto __AP2 = ut::allocate<T>(m, l);
    auto __BP3 = ut::allocate<T>(l, n);
    auto __BP4 = ut::allocate<T>(l, n);
    auto __AP5 = ut::allocate<T>(m, l);
    auto __AP6 = ut::allocate<T>(m, l);
    auto __BP6 = ut::allocate<T>(l, n);
    auto __AP7 = ut::allocate<T>(m, l);
    auto __BP7 = ut::allocate<T>(l, n);

    auto P1 = ut::allocate<T>(m, n);
    auto P2 = ut::allocate<T>(m, n);
    auto P3 = ut::allocate<T>(m, n);
    auto P4 = ut::allocate<T>(m, n);
    auto P5 = ut::allocate<T>(m, n);
    auto P6 = ut::allocate<T>(m, n);
    auto P7 = ut::allocate<T>(m, n);

    auto A11 = new T* [m];
    auto A12 = new T* [m];
    auto A21 = new T* [m];
    auto A22 = new T* [m];
    
    // Upper left.
    copy_block(m, 0, 0, A11, A);
    // Upper right.
    copy_block(m, 0, l, A12, A);
    // Lower left.
    copy_block(m, m, 0, A21, A);
    // Lower right.
    copy_block(m, m, l, A22, A);

    auto B11 = new T* [l];
    auto B12 = new T* [l];
    auto B21 = new T* [l];
    auto B22 = new T* [l];

    copy_block(l, 0, 0, B11, B);
    copy_block(l, 0, n, B12, B);
    copy_block(l, l, 0, B21, B);
    copy_block(l, l, n, B22, B);

    auto C11 = new T* [m];
    auto C12 = new T* [m];
    auto C21 = new T* [m];
    auto C22 = new T* [m];

    copy_block(m, 0, 0, C11, C);
    copy_block(m, 0, n, C12, C);
    copy_block(m, m, 0, C21, C);
    copy_block(m, m, n, C22, C);

#pragma omp task
    {
        // P1 := (A11 + A22) · (B11 + B22).
        add(m, l, A11, A22, __AP1);
        add(l, n, B11, B22, __BP1);
        recursive_multiply<T, N>(m, n, l, __AP1, __BP1, P1);
    }

#pragma omp task
    {
        // P2 := (A21 + A22) · B11.
        add(m, l, A21, A22, __AP2);
        recursive_multiply<T, N>(m, n, l, __AP2, B11, P2);
    }

#pragma omp task
    {
        // P3 := A11 · (B12 - B22).
        sub(l, n, B12, B22, __BP3);
        recursive_multiply<T, N>(m, n, l, A11, __BP3, P3);
    }

#pragma omp task
    {
        // P4 := A22 · (B21 - B11).
        sub(l, n, B21, B11, __BP4);
        recursive_multiply<T, N>(m, n, l, A22, __BP4, P4);
    }

#pragma omp task
    {
        // P5 := (A11 + A12) · B22.
        add(m, l, A11, A12, __AP5);
        recursive_multiply<T, N>(m, n, l, __AP5, B22, P5);
    }

#pragma omp task
    {
        // P6 := (A21 - A11) · (B11 + B12).
        sub(m, l, A21, A11, __AP6);
        add(l, n, B11, B12, __BP6);
        recursive_multiply<T, N>(m, n, l, __AP6, __BP6, P6);
    }

#pragma omp task
    {
        // P7 := (A12 - A22) · (B21 + B22).
        sub(m, l, A12, A22, __AP7);
        add(l, n, B21, B22, __BP7);
        recursive_multiply<T, N>(m, n, l, __AP7, __BP7, P7);
    }

#pragma omp taskwait

    for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C11[i][j] = P1[i][j] + P4[i][j] - P5[i][j] + P7[i][j];
            C12[i][j] = P3[i][j] + P5[i][j];
            C21[i][j] = P2[i][j] + P4[i][j];
            C22[i][j] = P1[i][j] - P2[i][j] + P3[i][j] + P6[i][j];
        }
    }

    ut::free(__AP1);
    ut::free(__BP1);
    ut::free(__AP2);
    ut::free(__BP3);
    ut::free(__BP4);
    ut::free(__AP5);
    ut::free(__AP6);
    ut::free(__BP6);
    ut::free(__AP7);
    ut::free(__BP7);

    ut::free(P1);
    ut::free(P2);
    ut::free(P3);
    ut::free(P4);
    ut::free(P5);
    ut::free(P6);
    ut::free(P7);

    // Shallow delete.
    delete[] A11;
    delete[] A12;
    delete[] A21;
    delete[] A22;

    delete[] B11;
    delete[] B12;
    delete[] B21;
    delete[] B22;

    delete[] C11;
    delete[] C12;
    delete[] C21;
    delete[] C22;
}
}

template<
    typename T, 
    size_t Threshold = 64 * 64 * 64,
    size_t Threads = 4,
    typename = std::enable_if_t<std::is_arithmetic_v<T>>
>
void multiply(size_t m, size_t n, size_t l, 
              T** A, T** B, T** C) {
    // assert((m & (m - 1)) == 0);
    // assert((n & (n - 1)) == 0);
    // assert((l & (l - 1)) == 0);
    
    omp_set_dynamic(0);
    omp_set_num_threads(Threads);

#pragma omp parallel
    {
#pragma omp single
        {
            detail::recursive_multiply<T, Threshold>(m, n, l, A, B, C);
        }
    }
}
}

#endif
