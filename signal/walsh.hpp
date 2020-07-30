//
// Created by ydrml on 2020/7/10.
//

#ifndef SIMULATION_WALSH_HPP
#define SIMULATION_WALSH_HPP

template<unsigned n>
struct walsh_t {
    constexpr static unsigned
        dim = 1u << (n - 1);
    
    static unsigned char memory[dim][dim];
    
    static void fill() {
        for (size_t i = 0; i < dim; ++i) {
            for (size_t j = 0; j < dim; ++j)
                get(i, j);
        }
    }
    
    static unsigned char get(unsigned i, unsigned j) {
        if (memory[i][j])
            return memory[i][j];
        else {
            auto div   = dim / 2;
            auto value = walsh_t<n - 1>::get(i % div, j % div);
            return memory[i][j] = i < div || j < div ? value : -value;
        }
    }
};

template<unsigned n>
unsigned char walsh_t<n>::memory[walsh_t<n>::dim][walsh_t<n>::dim]{{0}};

template<>
struct walsh_t<1> {
    constexpr static unsigned dim = 1;
    
    static unsigned char get(unsigned i, unsigned j) {
        return 1;
    }
};

#endif //SIMULATION_WALSH_HPP
