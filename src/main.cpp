#include <iostream>
#include "simu.hpp"

int main() {

    Kokkos::initialize();

    {
        Simu simu;

        simu.print();
        simu.save();
        for (int i = 0; i < 100; i++) {
            simu.tick();
            simu.print();
            simu.save();
        }
    }

    Kokkos::finalize();
    return 0;
}