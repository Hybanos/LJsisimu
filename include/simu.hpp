#include <vector>
#include <cmath>
#include <iostream>

#include "utils.hpp"

class Simu {
    private:
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;

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