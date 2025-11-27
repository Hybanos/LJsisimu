#include "simu.hpp"

Simu::Simu() {
    auto data = load_data();

    x.reserve(N_LOCAL);
    y.reserve(N_LOCAL);
    z.reserve(N_LOCAL);

    for (int i = 0; i < N_LOCAL; i++) {
        x[i] = data[i * 3];
        y[i] = data[i * 3 + 1];
        z[i] = data[i * 3 + 2];
    }
}

void Simu::compute_energy() {
    U = 0;

    for (int i = 0; i < N_LOCAL; i++) {
        for (int j = 0; j < i; j++) {
            U += compute_u_ij(i, j);
        }

        for (int j = i+1; j < N_LOCAL; j++) {
            U += compute_u_ij(i, j);
        }
    }

    U = U * 4;
}

double Simu::dist_squared(int p1, int p2) {
    return std::pow(x[p1] - x[p2], 2) + std::pow(y[p1] - y[p2], 2) + std::pow(z[p1] - z[p2], 2);
}

double Simu::compute_u_ij(int i, int j) {
    double r_ij_squared = dist_squared(i, j);
    double r_r_squared = r_star_squared / r_ij_squared;
    double u_ij = epsilon_star * 
        (std::pow(r_r_squared, 6) - 2 * std::pow(r_r_squared, 3));

    return u_ij;
}