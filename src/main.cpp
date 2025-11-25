#include <iostream>
#include <array>

#include "utils.hpp"

// #define N_total 1000
#define N_local 100

int main() {

    std::array<double, N_TOTAL * 3> data = load_data();

    for (int i = 0; i < N_TOTAL; i++) {
        std::cout << data[i * 3] << "  " << data[i * 3 + 1] << "  " << data[i * 3 + 2] << std::endl;
    }

    return 0;
}