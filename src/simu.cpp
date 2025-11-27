#include "simu.hpp"

Simu::Simu() {
    auto data = load_data();

    for (int i = 0; i < N_LOCAL; i++) {
        x[i] = data[i * 3];
        y[i] = data[i * 3 + 1];
        z[i] = data[i * 3 + 2];
    }
}

void Simu::compute_energy() {
    U = 0;
    for (int i = 0; i < N_LOCAL; i++) {
        for (int j = 0; j < i; j++) U += compute_u_ij(i, j);
        for (int j = i+1; j < N_LOCAL; j++) U += compute_u_ij(i, j);
    }
    U = U * 4;
}

double Simu::dist_squared(int i, int j) {
    return std::pow(x[i] - x[j], 2) + std::pow(y[i] - y[j], 2) + std::pow(z[i] - z[j], 2);
}

double Simu::compute_u_ij(int i, int j) {
    double r_ij_squared = dist_squared(i, j);
    double r_r_squared = r_star_squared / r_ij_squared;
    double u_ij = epsilon_star * 
        (std::pow(r_r_squared, 6) - 2 * std::pow(r_r_squared, 3));

    return u_ij;
}