#include <iostream>
#include "simu.hpp"

int main() {

    Simu simu;

    simu.compute_energy();
    // std::cout << simu.get_U() << std::endl;
    std::cout << simu.compute_u_ij(0, 1) << std::endl;
    std::cout << simu.dist_squared(0, 1) << std::endl;

    return 0;
}