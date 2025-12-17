#include <iostream>
#include "simu.hpp"

int main() {

    Kokkos::initialize();

    {
        Simu simu;
        simu.set_save_cond([](int i){ return (i % 10) == 0; });

        simu.print();
        simu.save();
        for (int i = 1; i < 10000; i++) {
            simu.step();
            if (i % 10 == 0) simu.print();
        }
    }

    Kokkos::finalize();
    return 0;
}