//
// Created by ydrml on 2020/5/26.
//

#ifndef FFT_BANDPASS_FILTER_T_HPP
#define FFT_BANDPASS_FILTER_T_HPP

#include <algorithm>
#include "../signal/complex_t.hpp"

/// 频域带通滤波器模板
template<class num_t, num_t _n, int _fs, int _f0, int _bw>
struct bandpass_filter_t {
    constexpr static num_t
        _i0 = _n * (_f0 - _bw / 2) / _fs,
        _i1 = _n * (_f0 + _bw / 2) / _fs,
        _i2 = _n - _i1 - 1,
        _i3 = _n - _i0 - 1;
    
    static_assert(_bw >= 2, "bandwidth < 2");
    static_assert(_i0 < _n / 2, "Nyquist!");
    static_assert(_i1 < _n / 2, "Nyquist!");
    
    /// 滤波
    static void filter(complex_t data[_n]) {
        std::fill(data, data + _i0, complex_t::zero);
        std::fill(data + _i1 + 1, data + _i2, complex_t::zero);
        std::fill(data + _i3 + 1, data + _n, complex_t::zero);
    }
};

template<class num_t, num_t n, int fs, int f0, int bw>
using bandpass = bandpass_filter_t<num_t, n, fs, f0, bw>;

#endif //FFT_BANDPASS_FILTER_T_HPP
