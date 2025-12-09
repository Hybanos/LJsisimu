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
        size_t ticks = 0;
        double timestep = 0.01;
        double R_cut_squared = std::pow(10, 2);

        array x = array("x", N_LOCAL);
        array y = array("y", N_LOCAL);
        array z = array("z", N_LOCAL);

        array x_loc = array("xj_loc", N_LOCAL);
        array y_loc = array("yj_loc", N_LOCAL);
        array z_loc = array("zj_loc", N_LOCAL);

        array fx = array("fx", N_LOCAL);
        array fy = array("fy", N_LOCAL);
        array fz = array("fz", N_LOCAL);

        images imgs = images("images", N_SYM);


        double U = 0;

        double r_star = 3.0;
        double r_star_squared = std::pow(r_star, 2);
        double epsilon_star = 0.2;

    public:
        // double dist_squared(int p1, int p2);
        Simu();

        void tick();
        void print();
        void save();
};
