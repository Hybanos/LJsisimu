#include <vector>
#include <cmath>
#include <iostream>

#include <Kokkos_Core.hpp>

#include "utils.hpp"

using array = Kokkos::View<double *, Kokkos::LayoutRight>;
using mat = Kokkos::View<double **, Kokkos::LayoutRight>;

class Simu {
    private:
        array x = array("x", N_LOCAL);
        array y = array("y", N_LOCAL);
        array z = array("z", N_LOCAL);

        double U = 0;

        double r_star = 3.0;
        double r_star_squared = std::pow(r_star, 2);
        double epsilon_star = 0.2;

    public:
        double dist_squared(int p1, int p2);
        double compute_u_ij(int i, int j);
        Simu();

        void compute_energy();
        double get_U() { return U; }
};