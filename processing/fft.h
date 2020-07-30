//
// Created by ydrml on 2020/5/26.
//

#ifndef FFT_FFT_H
#define FFT_FFT_H

#include <utility>
#include "../signal/complex_t.hpp"

/// 计算所需的 FFT 点数
template<class num_t, num_t n0, num_t n1>
constexpr num_t MIN_2_POW(num_t i = 1u) {
    return i < n0 + n1 ? MIN_2_POW<num_t, n0, n1>(i << 1u) : i;
}

/// 用于正变换的 ω<n,k>
template<auto _n>
complex_t omega(decltype(_n) k) {
    constexpr static auto t0 = 2 * M_PI / _n;
    
    float theta = t0 * k;
    return {std::cosf(theta), std::sinf(theta)};
}

/// 用于反变换的 ω<n,k>
template<auto _n>
complex_t omega_inverse(decltype(_n) k) {
    constexpr static auto t0 = 2 * M_PI / _n;
    
    float theta = t0 * k;
    return {std::cosf(theta), -std::sinf(theta)};
}

/// 基 2 快速傅里叶正变换
template<auto _n>
void fft(
    complex_t memory[_n],
    complex_t _omega(decltype(_n)) = omega<_n>
) {
    // 错序
    for (size_t i = 0, j = 0; i < _n; ++i) {
        if (i > j) std::swap(memory[i], memory[j]);
        for (size_t l = _n >> 1u; (j ^= l) < l; l >>= 1u);
    }
    // 变换
    for (size_t m = 1; m < _n; m <<= 1u) {
        auto        a = memory, b = a + m;
        const auto  s = _n / m / 2;
        for (size_t i = 0; i < s; ++i) {
            for (size_t j = 0; j < m; ++j, ++a, ++b)
                if (b->is_zero())
                    *b = *a;
                else {
                    const auto t = *b * _omega(s * j);
                    *b = *a - t;
                    *a += t;
                }
            a = b;
            b += m;
        }
    }
}

/// 基 2 快速傅里叶反变换
template<auto _n>
void ifft(complex_t memory[_n]) {
    fft<_n>(memory, omega_inverse<_n>);
    for (auto p = memory; p < memory + _n; ++p)
        *p /= _n;
}

#endif //FFT_FFT_H
