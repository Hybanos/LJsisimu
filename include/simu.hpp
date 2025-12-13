#pragma once 

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <Kokkos_Core.hpp>

#include "utils.hpp"

#define N_SYM 27

struct vec3 {
    double x = 0;
    double y = 0;
    double z = 0;

    vec3() {};
    vec3(double _x, double _y, double _z) : x{_x}, y{_y}, z{_z} {};
    void print() { std::cout << "vec3    x:" << x << " y:" << y << " z:" << z << std::endl; }
};

using array = Kokkos::View<double *, Kokkos::LayoutRight>;
using mat = Kokkos::View<double **, Kokkos::LayoutRight>;
using images = Kokkos::View<vec3 *, Kokkos::LayoutRight>;

class Simu {
    private:
        size_t steps = 0;

        // global positions
        array x = array("x", N_LOCAL);
        array y = array("y", N_LOCAL);
        array z = array("z", N_LOCAL);

        // local positions (global positions modulo'd to the center image)
        array x_loc = array("xj_loc", N_LOCAL);
        array y_loc = array("yj_loc", N_LOCAL);
        array z_loc = array("zj_loc", N_LOCAL);

        // Center of mass position
        double Px = 0.0;
        double Py = 0.0;
        double Pz = 0.0;

        // forces
        array fx = array("fx", N_LOCAL);
        array fy = array("fy", N_LOCAL);
        array fz = array("fz", N_LOCAL);

        // momentum
        array px = array("mx", N_LOCAL);
        array py = array("my", N_LOCAL);
        array pz = array("mz", N_LOCAL);

        // images offsets
        images imgs = images("images", N_SYM);

        double U = 0.0;
        double T = 0.0;
        double E_k = 0.0;

        double timestep = 1;
        double m = 18;
        double T_0 = 300.0;
        double R_const = 0.00199;
        double force_conversion_factor = 0.0001 * 4.186;
        double R_cut_squared = std::pow(10, 2);
        double r_star = 3.0;
        double r_star_squared = std::pow(r_star, 2);
        double epsilon_star = 0.2;

        void compute_kinetic_temp();
        void compute_center_of_mass();
        void calibrate_momentums();
        void calibrate_center_of_mass();
        void lennard_jones();
        void velocity_verlet();
    public:
        Simu();

        void step();
        void print();
        void save();
};
