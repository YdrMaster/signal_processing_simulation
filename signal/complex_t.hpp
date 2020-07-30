//
// Created by ydrml on 2020/5/21.
//

#ifndef FFT_COMPLEX_T_HPP
#define FFT_COMPLEX_T_HPP

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

struct complex_t {
    float re, im;
    
    const static complex_t zero;
    
    [[nodiscard]]
    inline float norm() const {
        return std::hypot(re, im);
    }
    
    [[nodiscard]]
    inline float arg() const {
        return std::atan2(im, re);
    }
    
    [[nodiscard]]
    inline complex_t conjugate() const {
        return {re, -im};
    }
    
    [[nodiscard]]
    inline complex_t normalize() const {
        auto l = norm();
        return l == 0 ? zero : complex_t{re / l, im / l};
    }
    
    [[nodiscard]]
    inline bool is_zero() const {
        return re == 0 && im == 0;
    }
    
    [[nodiscard]]
    inline complex_t operator+() const {
        return {re, im};
    }
    
    [[nodiscard]]
    inline complex_t operator-() const {
        return {-re, -im};
    }
    
    [[nodiscard]]
    inline complex_t operator+(const complex_t &others) const {
        return {re + others.re, im + others.im};
    }
    
    [[nodiscard]]
    inline complex_t operator-(const complex_t &others) const {
        return {re - others.re, im - others.im};
    }
    
    [[nodiscard]]
    inline complex_t operator*(const complex_t &others) const {
        return {re * others.re - im * others.im, re * others.im + im * others.re};
    }
    
    [[nodiscard]]
    inline complex_t operator/(const complex_t &others) const {
        float k = 1 / (others.re * others.re + others.im * others.im);
        return {(re * others.re + im * others.im) * k, (im * others.re - re * others.im) * k};
    }
    
    template<class num_t>
    [[nodiscard]]
    inline complex_t operator+(const num_t &others) const {
        return {re + others, im};
    }
    
    template<class num_t>
    [[nodiscard]]
    inline complex_t operator-(const num_t &others) const {
        return {re - others, im};
    }
    
    template<class num_t>
    [[nodiscard]]
    inline complex_t operator*(const num_t &others) const {
        return {re * others, im * others};
    }
    
    template<class num_t>
    [[nodiscard]]
    inline complex_t operator/(const num_t &others) const {
        return {re / others, im / others};
    }
    
    inline complex_t operator+=(const complex_t &others) {
        return *this = {re + others.re, im + others.im};
    }
    
    inline complex_t operator-=(const complex_t &others) {
        return *this = {re - others.re, im - others.im};
    }
    
    inline complex_t operator*=(const complex_t &others) {
        return *this = {re * others.re - im * others.im, re * others.im + im * others.re};
    }
    
    inline complex_t operator/=(const complex_t &others) {
        float k = 1 / (others.re * others.re + others.im * others.im);
        return *this = {(re * others.re + im * others.im) * k, (im * others.re - re * others.im) * k};
    }
    
    template<class num_t>
    inline complex_t &operator+=(const num_t &others) {
        return *this = {re + others, im};
    }
    
    template<class num_t>
    inline complex_t &operator-=(const num_t &others) {
        return *this = {re - others, im};
    }
    
    template<class num_t>
    inline complex_t &operator*=(const num_t &others) {
        return *this = {re * others, im * others};
    }
    
    template<class num_t>
    inline complex_t &operator/=(const num_t &others) {
        return *this = {re / others, im / others};
    }
};

const complex_t complex_t::zero = {0, 0};

#endif //FFT_COMPLEX_T_HPP
