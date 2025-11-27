#include "simu.hpp"

Simu::Simu() {
    auto data = load_data();

    for (int i = 0; i < N_LOCAL; i++) {
        x[i] = data[i * 3];
        y[i] = data[i * 3 + 1];
        z[i] = data[i * 3 + 2];
    }
}

void Simu::tick() {
    ticks++;

    U = 0;
    // for (int i = 0; i < N_LOCAL; i++) {
    Kokkos::parallel_for(
        "tick_outer",
        N_LOCAL,
        KOKKOS_LAMBDA(int i) {

            fx[i] = 0;
            fy[i] = 0;
            fz[i] = 0;
            for (int j = 0; j < i; j++) {
                if (i == j) continue;

                // Potential
                double r_ij_squared = dist_squared(i, j);
                double r_r_squared = r_star_squared / r_ij_squared;
                double u_ij = epsilon_star * 
                    (std::pow(r_r_squared, 6) - 2 * std::pow(r_r_squared, 3));

                U += u_ij;

                // Forces
                double tmp = -48 * epsilon_star * 
                    (std::pow(r_r_squared, 6) - std::pow(r_r_squared, 3));
                fx[i] += tmp * ((x[i] - x[j]) / r_ij_squared);
                fy[i] += tmp * ((y[i] - y[j]) / r_ij_squared);
                fz[i] += tmp * ((z[i] - z[j]) / r_ij_squared);
            }
        }
    );
    U = U * 4;

    // apply forces
    for (int i = 0; i < N_TOTAL; i++) {
        x[i] += fx[i];
        y[i] += fy[i];
        z[i] += fz[i];
    }
}

double Simu::dist_squared(int i, int j) {
    return std::pow(x[i] - x[j], 2) + std::pow(y[i] - y[j], 2) + std::pow(z[i] - z[j], 2);
}

void Simu::print() {
    std::cout << "iter: " << ticks << ", total energy: " << U << std::endl;
    for (int i = 0; i < 10; i++) {
        std::cout << "x: " << x[i] << " y: " << y[i] << " z: " << z[i] << std::endl;
    }
}

void Simu::save() {
    std::ofstream f; 
    f.open("saved/" + std::to_string(ticks) + ".data", std::ios::out | std::ios::binary);

    for (int i = 0; i < N_LOCAL; i++) {
        f.write(reinterpret_cast<const char *>(&x.data()[i]), sizeof(double));
        f.write(reinterpret_cast<const char *>(&y.data()[i]), sizeof(double));
        f.write(reinterpret_cast<const char *>(&z.data()[i]), sizeof(double));
    }

    f.close();
}