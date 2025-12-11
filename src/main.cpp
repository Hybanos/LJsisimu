#include <iostream>
#include "simu.hpp"

int main() {

    Kokkos::initialize();

    {
        Simu simu;

        simu.print();
        simu.save();
        for (int i = 1; i < 10000; i++) {
            simu.tick();
            if (i % 10 == 0) {
                simu.print();
                simu.save();
            }
        }
    }

    Kokkos::finalize();
    return 0;
}