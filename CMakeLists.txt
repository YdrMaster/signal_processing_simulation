cmake_minimum_required(VERSION 3.10)
project(simulation)

set(CMAKE_CXX_STANDARD 20)

add_executable(simulation main.cpp
        signal/complex_t.hpp

        processing/fft.h
        processing/bandpass_filter_t.hpp

        processing/multi_path.h
        signal/chirp.h
        signal/walsh.hpp

        processing/signal_process.h
        processing/simulation.h processing/static_check.h processing/noise.h)
