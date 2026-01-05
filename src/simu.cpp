#include "simu.hpp"

Simu::Simu() {
    auto data = load_data();
    f.open("out.data", std::ios::out | std::ios::binary);
    
    double n = N_LOCAL;
    double l = L;

    f.write(reinterpret_cast<const char *>(&n), sizeof(double));
    f.write(reinterpret_cast<const char *>(&l), sizeof(double));

    std::srand(0);
    Kokkos::parallel_for("init", policy,
        KOKKOS_LAMBDA(int i) {
            x[i] = data[i * 3]    ;
            y[i] = data[i * 3 + 1];
            z[i] = data[i * 3 + 2];

            // Init local pos
            x_loc[i] = std::fmod(x[i], L);
            y_loc[i] = std::fmod(y[i], L);
            z_loc[i] = std::fmod(z[i], L);
            x_loc[i] += L;
            y_loc[i] += L;
            z_loc[i] += L;
            x_loc[i] = std::fmod(x_loc[i], L);
            y_loc[i] = std::fmod(y_loc[i], L);
            z_loc[i] = std::fmod(z_loc[i], L);

            // init random momentums
            px[i] = ((double) std::rand() / RAND_MAX) * 2 - 1;
            py[i] = ((double) std::rand() / RAND_MAX) * 2 - 1;
            pz[i] = ((double) std::rand() / RAND_MAX) * 2 - 1;
        }
    );

    calibrate_momentums();
    calibrate_center_of_mass();
    calibrate_momentums();

    // init images coords
    int len = std::pow(N_SYM, 1.0/3.0);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            for (int k = 0; k < len; k++) {
                int img = i * 9 + j * 3 + k;
                int offset = len / 3;
                imgs[img] = vec3((i-offset)*L, (j-offset)*L, (k-offset)*L);
            }
        }
    }

    lennard_jones();
    compute_kinetic_temp();
    compute_center_of_mass();
}

void Simu::step() {
    steps++;

    auto t1 = std::chrono::high_resolution_clock().now();

    velocity_verlet_speed();
    velocity_verlet_position();
    lennard_jones();
    velocity_verlet_speed();
    compute_kinetic_temp();
    compute_center_of_mass();
    if (steps % m_step == 0) berendsen_thermostat();
    if (save_cond(steps)) save();

    auto t2 = std::chrono::high_resolution_clock().now();

    sps = 1e9 / (t2 - t1).count();
}

void Simu::lennard_jones() {

    U = 0;

    Kokkos::parallel_reduce("LJ", policy,
        KOKKOS_LAMBDA(int i, double &U_loc) {
            fx[i] = 0;
            fy[i] = 0;
            fz[i] = 0;
            for (int i_sym = 0; i_sym < N_SYM; i_sym++) {
                for (int j = 0; j < N_LOCAL; j++) {
                    if (i == j) continue;

                    // i to j_image distance
                    double xj_loc = x_loc[j] + imgs[i_sym].x;
                    double yj_loc = y_loc[j] + imgs[i_sym].y;
                    double zj_loc = z_loc[j] + imgs[i_sym].z;

                    double r_ij_squared = 
                        std::pow(x_loc[i] - xj_loc, 2) +
                        std::pow(y_loc[i] - yj_loc, 2) +
                        std::pow(z_loc[i] - zj_loc, 2);
                    if (r_ij_squared > R_cut_squared) continue; 

                    double r_r_squared = r_star_squared / r_ij_squared;

                    // Potential
                    if (i < j) {
                        double u_ij = std::pow(r_r_squared, 6) - 2 * std::pow(r_r_squared, 3);
                        U_loc += u_ij;
                    }

                    // Forces
                    double tmp = -48 * epsilon_star * 
                        (std::pow(r_r_squared, 6) - std::pow(r_r_squared, 3));
                    fx[i] += tmp * ((x_loc[i] - xj_loc) / r_ij_squared);
                    fy[i] += tmp * ((y_loc[i] - yj_loc) / r_ij_squared);
                    fz[i] += tmp * ((z_loc[i] - zj_loc) / r_ij_squared);
                }
            }
        }, U
    );

    U = U * epsilon_star * 2;
}

void Simu::velocity_verlet_speed() {
    Kokkos::parallel_for("velocity_verlet", policy,
        KOKKOS_LAMBDA(int i) {
            px[i] -= fx[i] * timestep * 0.5 * force_conversion_factor * magic_timestep_universe_fix;
            py[i] -= fy[i] * timestep * 0.5 * force_conversion_factor * magic_timestep_universe_fix;
            pz[i] -= fz[i] * timestep * 0.5 * force_conversion_factor * magic_timestep_universe_fix;
        }
    );
}

void Simu::velocity_verlet_position() {
    Kokkos::parallel_for("velocity_verlet", policy,
        KOKKOS_LAMBDA(int i) {
            x[i] += (px[i] * timestep) / m;
            y[i] += (py[i] * timestep) / m;
            z[i] += (pz[i] * timestep) / m;

            // mod positions in main image
            x_loc[i] = std::fmod(x[i], L);
            y_loc[i] = std::fmod(y[i], L);
            z_loc[i] = std::fmod(z[i], L);
            x_loc[i] += L;
            y_loc[i] += L;
            z_loc[i] += L;
            x_loc[i] = std::fmod(x_loc[i], L);
            y_loc[i] = std::fmod(y_loc[i], L);
            z_loc[i] = std::fmod(z_loc[i], L);
        }
    );
}

void Simu::compute_kinetic_temp() {
    double sum = 0.0;
    Kokkos::parallel_reduce("compute_kinetic_temp", policy,
            KOKKOS_LAMBDA(int i, double &lsum) {
            lsum += (px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i]);
        }, sum
    );
    sum /= m;
    E_k = (1.0 / (2.0 * force_conversion_factor)) * sum;

    T = 1 / (N_dl * R_const) * E_k;
}

void Simu::calibrate_momentums() {
    compute_kinetic_temp();
    double ratio = N_dl * R_const * T_0 / E_k;
    // idk if this square root is physically accurate but i'm going insane without it
    ratio = std::sqrt(ratio);
    Kokkos::parallel_for("velocity_verlet_2", policy,
        KOKKOS_LAMBDA(int i) {
            px[i] *= ratio;
            py[i] *= ratio;
            pz[i] *= ratio;
        }
    );
    compute_kinetic_temp();
}

void Simu::compute_center_of_mass() {
    Px = 0.0;
    Py = 0.0;
    Pz = 0.0;

    for (int i = 0; i < N_LOCAL; i++) {
        Px += px[i];
        Py += py[i];
        Pz += pz[i];
    }

    Px = Px / N_TOTAL;
    Py = Py / N_TOTAL;
    Pz = Pz / N_TOTAL;
}

void Simu::berendsen_thermostat() {
    Kokkos::parallel_for("berendsen_thermostat", policy,
        KOKKOS_LAMBDA(int i) {
            px[i] = px[i] + gamma * (T_0 / T - 1) * px[i];
            py[i] = py[i] + gamma * (T_0 / T - 1) * py[i];
            pz[i] = pz[i] + gamma * (T_0 / T - 1) * pz[i];
        }
    );
}

void Simu::calibrate_center_of_mass() {
    compute_center_of_mass();
    Kokkos::parallel_for("calibrate_center_of_mass", policy,
        KOKKOS_LAMBDA(int i) {
            px[i] -= Px;
            py[i] -= Py;
            pz[i] -= Pz;
        }
    );
    compute_center_of_mass();
}

void Simu::print() {
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "iter: " << steps << ", sps: " << sps << std::endl;
    std::cout << ", total energy: " << U + E_k << ", U: " << U << ", E_k: " << E_k << ", temp: " << T <<std::endl;
    double xx = 0, yy = 0, zz = 0;
    for (int i = 0; i < N_LOCAL; i++) {
        xx += fx[i];
        yy += fy[i];
        zz += fz[i];
    }
    std::cout << "f total: " << xx << ", " << yy << ", " << zz << std::endl;
    for (int i = 0; i < std::min(5, N_LOCAL); i++) {
        std::cout << "x: " << x[i] << " y: " << y[i] << " z: " << z[i] << std::endl;
        std::cout << "fx: " << fx[i] << " fy: " << fy[i] << " fz: " << fz[i] << std::endl;
        std::cout << "px: " << px[i] << " py: " << py[i] << " pz: " << pz[i] << std::endl;
    }
}

void Simu::save() {

    double iter = steps;

    f.write(reinterpret_cast<const char *>(&iter), sizeof(double));
    f.write(reinterpret_cast<const char *>(&U), sizeof(double));
    f.write(reinterpret_cast<const char *>(&T), sizeof(double));
    f.write(reinterpret_cast<const char *>(&E_k), sizeof(double));

    f.write(reinterpret_cast<const char *>(&Px), sizeof(double));
    f.write(reinterpret_cast<const char *>(&Py), sizeof(double));
    f.write(reinterpret_cast<const char *>(&Pz), sizeof(double));

    for (int i = 0; i < N_LOCAL; i++) {
        f.write(reinterpret_cast<const char *>(&x.data()[i]), sizeof(double));
        f.write(reinterpret_cast<const char *>(&y.data()[i]), sizeof(double));
        f.write(reinterpret_cast<const char *>(&z.data()[i]), sizeof(double));
    }

    f.flush();
}
