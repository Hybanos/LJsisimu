#include <iostream>
#include "simu.hpp"

int main() {

    Kokkos::initialize();

    {
        Simu simu;

        simu.print();
        simu.save();
        for (int i = 1; i < 100000; i++) {
            simu.step();
            if (i % 100 == 0) {
                simu.print();
                simu.save();
            }
        }
    }

    Kokkos::finalize();
    return 0;
}