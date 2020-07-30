//
// Created by ydrml on 2020/6/10.
//

#ifndef FFT_SIMULATION_H
#define FFT_SIMULATION_H

#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>
#include <functional>
#include <algorithm>
#include <numeric>

#include "../signal/complex_t.hpp"
#include "signal_process.h"

#define SAVE_SIGNAL(PATH, S) \
save_signal(PATH, S, [](std::ofstream &file, typename decltype(S)::value_type x) { file << x << std::endl; })

#define SAVE_SIGNAL_TF(PATH, S, TF) \
save_signal(PATH, S, [](std::ofstream &file, typename decltype(S)::value_type x) { file << (TF) << std::endl; })

#define SAVE_SIGNAL_FORMAT(PATH, S, TF) \
save_signal(PATH, S, [](std::ofstream &file, typename decltype(S)::value_type x) { file << TF; })

template<auto _length, class sample_t = float>
std::vector<sample_t> build_signal(float fs, std::function<sample_t(float)> const &function) {
    auto signal = std::vector<sample_t>(_length);
    auto ts     = 1 / fs;
    auto t      = .0f;
    
    for (auto p = signal.begin(); p < signal.end(); ++p)
        *p = function(t), t += ts;
    
    return signal;
}

template<class code_t, class sample_t>
std::vector<sample_t> encode(
    std::vector<code_t> const &code,
    std::unordered_map<code_t, std::vector<sample_t>> const &map
) {
    std::vector<sample_t> result;
    size_t                length = 0;
    
    for (auto &c:code) {
        auto slice = map.at(c);
        result.resize(length + slice.size());
        std::copy(slice.begin(), slice.end(), result.begin() + length);
        length += slice.size();
    }
    
    return result;
}

/**
 * 从文件读取数字信号
 * @tparam sample_t 采样点类型
 * @tparam _length 信号长度
 * @param file_name 文件名
 * @param function 解码函数
 * @return 数字信号
 */
template<class sample_t, long _length>
std::vector<sample_t> load_signal(
    const std::string &file_name,
    bool function(std::ifstream &, sample_t &)
) {
    std::vector<sample_t> signal(_length, sample_t{});
    std::ifstream         file(file_name);
    sample_t              value;
    
    for (auto p = signal.begin(); p < signal.end() && function(file, value); ++p)
        *p = value;
    return signal;
}

template<auto _size>
std::vector<float> send_signal(std::vector<float> const &x0) {
    // 加载发射端冲激响应原始数据
    auto t0 = load_signal<float, 2048>(
        "C:\\Users\\ydrml\\Desktop\\数据\\2048_1M.txt",
        [](std::ifstream &file, float &value) {
            return (bool) (file >> value);
        });
    
    // 处理发射端冲激响应（去除时延，去除直流分量）
    auto      mean = std::accumulate(t0.begin(), t0.end(), .0f) / t0.size();
    for (auto &x:t0) x -= mean;
    
    // 计算发射信号
    return convolve<_size>(x0, t0);
}

/**
 * 保存数字信号到文件
 * @tparam sample_t 采样点类型
 * @param file_name 文件名
 * @param signal 信号
 * @param function 编码函数
 */
template<class collector_t>
void save_signal(
    std::string const &file_name,
    collector_t const &signal,
    std::function<void(std::ofstream &, typename collector_t::value_type)>
    const &function
) {
    std::ofstream file(file_name);
    for (auto     x:signal) function(file, x);
    file.close();
}

#endif //FFT_SIMULATION_H
