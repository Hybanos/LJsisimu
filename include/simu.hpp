#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <Kokkos_Core.hpp>

#include "traits/Kokkos_IterationPatternTrait.hpp"
#include "utils.hpp"

struct vec3 {
    double x = 0;
    double y = 0;
    double z = 0;

    vec3() {};
    vec3(double _x, double _y, double _z) : x{_x}, y{_y}, z{_z} {};
};

using array = Kokkos::View<double *, Kokkos::LayoutRight>;
using mat = Kokkos::View<double **, Kokkos::LayoutRight>;
using images = Kokkos::View<vec3 ***, Kokkos::LayoutRight>;

class Simu {
    private:
        size_t ticks = 0;
        double timestep = 0.005;

        array x = array("x", N_LOCAL);
        array y = array("y", N_LOCAL);
        array z = array("z", N_LOCAL);

        array fx = array("fx", N_LOCAL);
        array fy = array("fy", N_LOCAL);
        array fz = array("fz", N_LOCAL);

        images imgs = images("images", 3, 3, 3);


        double U = 0;

        double r_star = 3.0;
        double r_star_squared = std::pow(r_star, 2);
        double epsilon_star = 0.2;

    public:
        double dist_squared(int p1, int p2);
        Simu();

        void tick();
        void print();
        void save();
};
