#include <iostream>
#include <random>
#include <sstream>
#include <numeric>

#include "../signal/complex_t.hpp"
#include "fft.h"
#include "bandpass_filter_t.hpp"
#include "simulation.h"

#define F(X) static_cast<float>(X)

template<class t>
using vec = std::vector<t>;

template<class iterator_t, class separator_t>
std::string join_to_string(iterator_t, iterator_t, separator_t);

template<size_t _size>
struct spectrum_t {
    float     fs;
    complex_t *data;
    
    [[nodiscard]]
    inline complex_t operator[](float f) const {
        auto i = static_cast<int>(std::round(f / fs * _size));
        return data[i];
    }
};

std::vector<float> send_signal(std::vector<float> const &x0);

void simulate();

int main() {
    simulate();
    return 0;
}


template<class iterator_t, class separator_t>
std::string join_to_string(iterator_t begin, iterator_t end, separator_t separator) {
    std::stringstream builder;
    builder << *begin++;
    while (begin < end) builder << separator << *begin++;
    return builder.str();
}

void simulate() {
    // 构造激励信号
    auto x0 = build_signal<float, (long) (1e-3 * 600e3)>(
        600e3,
        [](float t) {
            return F(std::sin(2 * M_PI * 40e3f * t));
        });
    
    // 计算发射信号
    auto x1 = send_signal(x0);
    
    // 构造多径信道描述
    auto path_info = vec<path_info_t>();
    {
        auto        src = complex_t{0, 0};
        auto        tgt = complex_t{1, 0};
        auto        ob0 = complex_t{0, -.2};
        auto        ob1 = complex_t{1, -.2};
        auto        s0  = (tgt - src).norm();
        for (size_t i   = 0; i <= 100; ++i) {
            auto ob = ob0 * (1 - F(i) / 100) + ob1 * F(i) / 100;
            path_info.push_back({.3, (src - ob).norm() + (tgt - ob).norm() - s0, 1});
        }
    }
    
    // 构造多径信道响应
    auto multi_path = build_multi_path_response<4096>(
        600e3f, 20.048f * F(std::sqrt(15 + 273.15)), path_info);
    
    // 通过多径信道
    auto y0 = convolve<float, 4096>(x1, {1});
    
    // 加噪
    auto y1 = vec<float>(y0);
    add_noise(y1, 100_db);
    
    // 接收（高通滤波）
    {
        auto      temp = y1[0];
        for (auto &y:y1) {
            std::swap(y, temp);
            y = temp - y;
        }
    }
    
    auto buffer0 = fft<size_t, float, 8192>(x1);
    auto buffer1 = fft<size_t, float, 8192>(y1);
    
    bandpass<size_t, 8192, 600, 50, 24>::filter(buffer0.data());
    
    // 频谱乘
    for (auto p = buffer1.begin(),
              q = buffer0.begin();
         p < buffer1.end(); ++p, ++q)
        if (!p->is_zero()) *p *= q->normalize();
    // 原地反 fft
    ifft<size_t, 8192>(buffer1.data());
    
    SAVE_SIGNAL("../data/x0.txt", float, x0, x);
    SAVE_SIGNAL("../data/x1.txt", float, x1, x);
    SAVE_SIGNAL("../data/mp.txt", float, multi_path, x);
    SAVE_SIGNAL("../data/y0.txt", float, y0, x);
    SAVE_SIGNAL("../data/y1.txt", float, y1, x);
    SAVE_SIGNAL("../data/xcorr.txt", complex_t, buffer1, x.re);
}
