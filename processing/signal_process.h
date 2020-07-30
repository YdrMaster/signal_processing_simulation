//
// Created by ydrml on 2020/7/9.
//

#ifndef SIMULATION_SIGNAL_PROCESS_H
#define SIMULATION_SIGNAL_PROCESS_H

#include <vector>
#include <array>
#include <algorithm>

#include "../signal/complex_t.hpp"
#include "static_check.h"
#include "fft.h"

/**
 * 切片并复制向量
 * @tparam sample_t 采样点类型
 * @param vec 原向量
 * @param begin 起点下标
 * @param length 长度
 * @return 切片的副本
 */
template<class sample_t>
std::vector<sample_t> slice(
    const std::vector<sample_t> &vec,
    size_t begin,
    size_t length = -1
) {
    auto iterator = std::begin(vec) + begin;
    auto end      = length < 0 ? std::end(vec) : iterator + length;
    return std::vector<sample_t>(iterator, end);
}

/**
 * 标准化
 * @tparam sample_t 采样类型
 * @param vec 信号存储
 * @param target 目标最大值（绝对值）
 */
template<class sample_t>
void normalize(
    std::vector<sample_t> &vec,
    sample_t target = 1
) {
    sample_t  max = vec.front();
    for (auto p   = vec.begin() + 1; p < vec.end(); ++p)
        if (std::abs(*p) > max) max = std::abs(*p);
    for (auto &n:vec)
        n = n * target / max;
}

/**
 * 重采样
 * @remarks 重采样用于把某一采样率的信号用新的采样率重新采样，可以进行升采样，也可以进行降采样。
 *          重采样的原理是先大倍数升采样，再在近似新采样率下抽取，
 *          因此，仅当新采样率与原采样率有整倍数关系，重采样才保证准确性。
 *          否则，倍率越大，采样越准。
 * @tparam times 处理倍率
 * @tparam size0 原信号长度（确保 `signal.size() < size0`）
 * @tparam size1 新信号长度（点数不够将补 0）
 * @param signal 原信号
 * @param f0 原采样率
 * @param f1 新采样率
 * @return 重采样信号
 */
template<auto times, auto size0, auto size1>
std::vector<float> resample(
    std::vector<float> const &signal,
    float f0,
    float f1
) {
    auto n_downsampling = std::lroundf(f0 * times / f1);
    auto enlarged       = fft_real<size0>(signal);
    enlarged.resize(times * size0);
    
    for (auto p = enlarged.begin() + size0 / 2, q = enlarged.end() - size0 / 2; q < enlarged.end(); ++p, ++q) {
        auto temp = *p;
        *p = *q;
        *q = temp;
    }
    ifft<times * size0>(enlarged.data());
    
    auto      target = std::vector<float>(size1, 0);
    for (auto i      = 0; i < size1; ++i) {
        auto j = n_downsampling * i;
        if (j >= enlarged.size()) break;
        target[i] = enlarged[j].re;
    }
    return target;
}

/**
 * 用 FFT 变换实信号
 * @tparam _size_per_group 每组 FFT 长度，必须是 2 的整数次幂
 * @tparam _group_count FFT 分组数量
 * @param signal 原信号
 * @return 变换
 */
template<auto _size_per_group, auto _group_count = 1>
std::vector<complex_t> fft_real(std::vector<float> const &signal) {
    constexpr static auto _size = _group_count * _size_per_group;
    static_assert(_group_count > 0);
    static_assert(check_power_2<_size>(), "size is not power of 2");
    
    std::vector<complex_t> spectrum(_size, complex_t::zero);
    
    if constexpr (_group_count == 1) {
        std::transform(signal.begin(), signal.end(), spectrum.begin(),
                       [](float z) -> complex_t { return {z, 0}; });
        fft<_size_per_group>(spectrum.data());
    } else {
        complex_t parts[_group_count][_size_per_group]{{}};
        
        { // 分组
            complex_t *iterators[_group_count];
            
            for (size_t i = 0; i < _group_count; ++i)
                iterators[i] = parts[i];
            
            for (auto ptr = signal.begin(); ptr < signal.end();)
                for (auto &it : iterators) *it++ = {*ptr++, 0};
        }
        
        // 变换
        for (auto v : parts) fft<decltype(_size_per_group), _size_per_group>(v);
        
        // 合并
        for (size_t i = 0; i < _group_count; ++i)
            for (size_t j = 0; j < _size_per_group; ++j) {
                auto n = i * _size_per_group + j;
                spectrum[n] = parts[0][j];
                for (size_t k = 1; k < _group_count; ++k)
                    spectrum[n] += omega<_size>(n * k) * parts[k][j];
            }
    }
    
    return spectrum;
}

/**
 * 快速卷积
 * @tparam _size 计算长度
 * @param a 信号a
 * @param b 信号b
 * @return 卷积信号
 */
template<auto _size>
std::vector<float> convolve(
    std::vector<float> const &a,
    std::vector<float> const &b
) {
    static_assert(check_power_2<_size>(), "size is not power of 2");
    
    auto      fa = fft_real<_size>(a),
              fb = fft_real<_size>(b);
    for (auto p  = fa.begin(),
              q  = fb.begin();
         p < fa.end(); ++p, ++q)
        *p *= *q;
    ifft<_size>(fa.data());
    
    std::vector<float> result(_size);
    std::transform(fa.begin(), fa.end(), result.begin(),
                   [](complex_t z) { return z.re; });
    return result;
}

/// 希尔伯特变换
template<auto _size>
std::vector<complex_t> hilbert(std::vector<float> const &x) {
    static_assert(check_power_2<_size>(), "size is not power of 2");
    
    // 生成超前 90° 的信号（虚部）
    auto result = fft_real<_size>(x);
    {
        auto p = result.begin();
        ++p; // 避开 0 频率点，前一半，正频率部分，超前 90°
        while (p < result.begin() + _size / 2)
            *p++ = {p->im, -p->re};
        ++p; // 避开 0 频率点，后一半，负频率部分，滞后 90°
        while (p < result.end())
            *p++ = {-p->im, p->re};
    }
    ifft<_size>(result.data());
    
    // 与原信号合并为复信号
    result.resize(x.size());
    auto p = x.begin();
    auto q = result.begin();
    while (p < x.end()) *q++ = {*p++, q->re};
    return result;
}

/**
 * 互相关（静态部分）=== 傅里叶变换并取共轭
 * @tparam _size 变换长度，必须是 2 的整数次幂
 * @param signal 原信号
 * @return 相关滤波器谱
 */
template<auto _size>
std::vector<complex_t> xcorr_init(std::vector<float> const &signal) {
    static_assert(check_power_2<_size>(), "size is not power of 2");
    
    auto spectrum = fft_real<_size>(signal);
    
    for (auto &it : spectrum) it.im = -it.im;
    
    return spectrum;
}

/**
 * 互相关（动态部分）=== 白化并乘以滤波器谱
 * @tparam _size 变换长度，必须是 2 的整数次幂
 * @param filter 相关滤波器谱，来自 `xcorr_init`
 * @param signal 原信号
 */
template<auto _size>
void xcorr(std::vector<complex_t> const &filter, std::vector<float> &signal) {
    static_assert(check_power_2<_size>(), "size is not power of 2");
    
    auto spectrum = fft_real<_size>(signal);
    auto p        = filter.begin();
    auto q        = spectrum.begin();
    while (p < filter.end()) *q++ = q->normalize() * *p++;
    
    ifft<_size>(spectrum.data());
    std::transform(spectrum.begin(), spectrum.end(), signal.begin(), [](complex_t it) { return it.re; });
}

#endif // SIMULATION_SIGNAL_PROCESS_H
