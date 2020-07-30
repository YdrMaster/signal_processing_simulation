//
// Created by ydrml on 2020/6/15.
//

#ifndef SIMULATION_PAM_H
#define SIMULATION_PAM_H

#include <tuple>

#include "../signal/complex_t.hpp"

inline std::pair<float, float> pam(
    float fs,
    float f1,
    float f2,
    const float *data,
    size_t length
) {
    auto      omega1 = 2 * M_PI * f1 / fs,
              omega2 = 2 * M_PI * f2 / fs;
    auto      a      = complex_t::zero,
              b      = complex_t::zero;
    for (auto i      = 0; i < length; ++i) {
        auto k = *data++;
        a += {k * std::cosf(omega1 * i), k * std::sinf(omega1 * i)};
        b += {k * std::cosf(omega2 * i), k * std::sinf(omega2 * i)};
    }
    return {std::fmax(a.norm(), b.norm()) * 2 / length, (a / b).arg()};
}

template<class sample_t>
std::pair<float, float> pam(
    float fs,
    float f1,
    float f2,
    const sample_t *data,
    size_t length,
    float block(sample_t)
) {
    auto      omega1 = 2 * M_PI * f1 / fs,
              omega2 = 2 * M_PI * f2 / fs;
    auto      a      = complex_t::zero,
              b      = complex_t::zero;
    for (auto i      = 0; i < length; ++i) {
        auto k = block(*data++);
        a += {k * std::cosf(omega1 * i), k * std::sinf(omega1 * i)};
        b += {k * std::cosf(omega2 * i), k * std::sinf(omega2 * i)};
    }
    return {std::fmax(a.norm(), b.norm()) * 2 / length, (a / b).arg()};
}

#endif //SIMULATION_PAM_H
