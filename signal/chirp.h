//
// Created by ydrml on 2020/7/9.
//

#ifndef SIMULATION_CHIRP_H
#define SIMULATION_CHIRP_H

#include <functional>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

std::function<float(float)>
chirp_linear(float f0_hz, float f1_hz, float t_s, float phi0 = 0) {
    return [k = (f1_hz - f0_hz) / t_s / 2, f0_hz, phi0](float t) {
        return std::sinf(2 * M_PI * (k * t + f0_hz) * t + phi0);
    };
}

#endif // SIMULATION_CHIRP_H
