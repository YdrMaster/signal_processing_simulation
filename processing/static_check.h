//
// Created by ydrml on 2020/7/20.
//

#ifndef SIMULATION_STATIC_CHECK_H
#define SIMULATION_STATIC_CHECK_H

template<auto value>
constexpr bool check_power_2() {
    unsigned char      k    = 0;
    unsigned long long copy = value;
    while (copy) {
        k += copy & 1u;
        copy >>= 1u;
    }
    return !value || k == 1;
}

#endif //SIMULATION_STATIC_CHECK_H
