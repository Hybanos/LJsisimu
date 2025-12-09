#include <iostream>
#include "simu.hpp"

int main() {

    Kokkos::initialize();

    {
        Simu simu;

        simu.print();
        simu.save();
        for (int i = 0; i < 1000000; i++) {
            simu.tick();
            if (i % 50 == 0) {
                simu.print();
                simu.save();
            }
        }
    }

    Kokkos::finalize();
    return 0;
}