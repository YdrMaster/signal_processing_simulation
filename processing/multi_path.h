//
// Created by ydrml on 2020/7/9.
//

#ifndef SIMULATION_MULTI_PATH_H
#define SIMULATION_MULTI_PATH_H

#include <vector>

/**
 * 多径声径信息
 * @param a 路径增益
 * @param ds 额外路程（声径比第一径长的路程，米）
 * @param reflect_times 反射次数
 */
struct path_info_t {
    float a;
    float ds;
    int   reflect_times;
};

/**
 * 从多径信道描述构造时域响应
 * @tparam _length 长度
 * @param fs 采样率
 * @param c 声速
 * @param path_info 信道描述
 * @return 时域响应（单位脉冲响应信号）
 */
template<size_t _length>
std::vector<float> build_multi_path_response(
    float fs, float c,
    const std::vector<path_info_t> &path_info
) {
    auto signal = std::vector<float>(_length, 0);
    signal[0] = 1;
    for (auto info:path_info)
        signal[static_cast<int>(info.ds / c * fs )]
            = info.reflect_times % 2 ? -info.a : info.a;
    return signal;
}

#endif //SIMULATION_MULTI_PATH_H
