//
// Created by ydrml on 2020/7/23.
//

#ifndef SIMULATION_NOISE_H
#define SIMULATION_NOISE_H

#include <vector>
#include <cmath>
#include <numeric>
#include <random>

#include "../signal/complex_t.hpp"

struct db_t {
    float value;
    
    [[nodiscard]]
    inline float to_float() const {
        return std::powf(10, value / 10);
    }
    
    [[nodiscard]]
    inline db_t operator-() const {
        return {-value};
    }
};

/// 构造信噪比值的分贝字面值
constexpr db_t operator ""_db(long double db) {
    return {static_cast<float>(db)};
}

/// 构造信噪比值的分贝字面值
constexpr db_t operator ""_db(unsigned long long db) {
    return {static_cast<float>(db)};
}

/// 计算实信号能量
/// \tparam sample_t 信号类型
/// \param signal 信号
/// \return 能量值
template<class sample_t>
float energy(std::vector<sample_t> const &signal) {
    return std::accumulate(
        signal.begin(), signal.end(), .0f,
        [](float sum, float x) { return sum + x * x; }
    ) / signal.size();
}

/// 为信号加噪
/// \tparam snr_t 信噪比类型
/// \param signal 信号
/// \param snr 信噪比
template<class snr_t>
void add_noise(std::vector<float> &signal, snr_t snr) {
    float sigma = std::sqrtf(energy(signal) / snr);
    if (sigma == 0) return;
    
    std::random_device         rd{};
    std::mt19937               gen{rd()};
    std::normal_distribution<> d{0, sigma};
    
    for (auto &x:signal) x += d(gen);
}

void add_noise(std::vector<float> &signal, db_t snr) {
    add_noise(signal, snr.to_float());
}


#endif // SIMULATION_NOISE_H
