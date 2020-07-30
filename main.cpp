#include <iostream>
#include <chrono>

#include "processing/noise.h"
#include "processing/simulation.h"
#include "signal/chirp.h"
#include "signal/walsh.hpp"

template<class t>
using vec = std::vector<t>;
using w3 = walsh_t<3>;

int main() {
    auto x0 = build_signal<1000>(1e6, chirp_linear(39e3f, 60.5e3f, 1e-3f));
    
    for (auto i = 0; i < x0.size(); ++i) x0[i] *= 1 + 5.0f * i / x0.size();
    
    auto normalized = std::vector<float>(x0);
    normalize(normalized, 1000.0f);
    SAVE_SIGNAL_TF("../data/shorted_to_send.txt", normalized, static_cast<unsigned short>(x + 2048));
    
    // 重采样到 808 周期对应的参考信号
    auto resampled = resample<64, 8192, 256>(send_signal<8192>(x0), 1e6f, 1e8f / 808);
    normalize(resampled, 1024.0f);
    SAVE_SIGNAL_FORMAT("../data/shorted_for_reference.txt", resampled, static_cast<short>(x) << ',');
    return 0;
}

void encode() {
    w3::fill();
    
    std::vector<uint8_t>                    code(w3::dim);
    std::unordered_map<uint8_t, vec<float>> map;
    map[1]   = build_signal<200>(200e3, chirp_linear(39e3f, 61e3f, 1e-3f));
    map[255] = build_signal<200>(200e3, chirp_linear(61e3f, 39e3f, 1e-3f));
    
    std::copy(w3::memory[0], w3::memory[1], code.begin());
    auto x0 = encode(code, map);
    std::copy(w3::memory[2], w3::memory[3], code.begin());
    auto x1 = encode(code, map);
    
    auto y0 = send_signal<1024>(x0);
    auto y1 = send_signal<1024>(x1);
    
    std::vector<float> yy(4096, 0);
    
    {
        decltype(yy.begin()) p;
        p = yy.begin() + 400;
        for (auto q = y0.begin(); q < y0.end(); ++p, ++q) *p += *q;
        p = yy.begin() + 800;
        for (auto q = y1.begin(); q < y1.end(); ++p, ++q) *p += *q;
        normalize<float>(yy, 4096);
    }
    
    auto Y0 = std::vector(yy);
    auto Y1 = std::vector(yy);
    
    xcorr<4096>(xcorr_init<4096>(y0), Y0);
    xcorr<4096>(xcorr_init<4096>(y1), Y1);
    
    SAVE_SIGNAL("../data/yy.txt", yy);
    SAVE_SIGNAL("../data/y0.txt", y0);
    SAVE_SIGNAL("../data/y1.txt", y1);
    SAVE_SIGNAL("../data/Y0X.txt", Y0);
    SAVE_SIGNAL("../data/Y1X.txt", Y1);
}
