#include "simu.hpp"

Simu::Simu() {
    auto data = load_data();

    std::srand(std::time({}));
    for (int i = 0; i < N_LOCAL; i++) {
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
    compute_kinetic_temp();
    normalize_momentums();
    compute_kinetic_temp();

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
}

void Simu::tick() {
    ticks++;

    U = 0;

    for (int i = 0; i < N_LOCAL; i++) {
        fx_tmp[i] = fx[i];
        fy_tmp[i] = fy[i];
        fz_tmp[i] = fz[i];
        fx[i] = 0;
        fy[i] = 0;
        fz[i] = 0;
    }

    // LJ
    for (int i_sym = 0; i_sym < N_SYM; i_sym++) {
        for (int i = 0; i < N_LOCAL; i++) {

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
                double u_ij = std::pow(r_r_squared, 6) - 2 * std::pow(r_r_squared, 3);
                U += u_ij;

                // Forces
                double tmp = -48 * epsilon_star * 
                    (std::pow(r_r_squared, 6) - std::pow(r_r_squared, 3));
                fx[i] -= tmp * ((x_loc[i] - xj_loc) / r_ij_squared);
                fy[i] -= tmp * ((y_loc[i] - yj_loc) / r_ij_squared);
                fz[i] -= tmp * ((z_loc[i] - zj_loc) / r_ij_squared);
            }
        }
    }
    U = U * epsilon_star * 2;

    // Velocity-verlet
    for (int i = 0; i < N_LOCAL; i++) {
        double x_tp1 = x[i] + px[i] / m * timestep + fx[i] * timestep * timestep * 0.5 * force_conversion_factor;
        double y_tp1 = y[i] + py[i] / m * timestep + fy[i] * timestep * timestep * 0.5 * force_conversion_factor;
        double z_tp1 = z[i] + pz[i] / m * timestep + fz[i] * timestep * timestep * 0.5 * force_conversion_factor;

        double vx_tp1 = px[i] / m + (fx[i] + fx_tmp[i]) * timestep * 0.5 * force_conversion_factor;
        double vy_tp1 = py[i] / m + (fy[i] + fy_tmp[i]) * timestep * 0.5 * force_conversion_factor;
        double vz_tp1 = pz[i] / m + (fz[i] + fz_tmp[i]) * timestep * 0.5 * force_conversion_factor;

        x[i] = x_tp1;
        y[i] = y_tp1;
        z[i] = z_tp1;

        px[i] = vx_tp1 * m;
        py[i] = vy_tp1 * m;
        pz[i] = vz_tp1 * m;

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

    // Barycenter momentum conservation
    double Px = 0;
    double Py = 0;
    double Pz = 0;

    for (int i = 0; i < N_LOCAL; i++) {
        Px += px[i];
        Py += py[i];
        Pz += pz[i];
    }

    for (int i = 0; i < N_TOTAL; i++) {
        px[i] = px[i] - Px / N_TOTAL;
        py[i] = py[i] - Py / N_TOTAL;
        pz[i] = pz[i] - Pz / N_TOTAL;
    }

    compute_kinetic_temp();
    normalize_momentums();
    compute_kinetic_temp();
}

void Simu::compute_kinetic_temp() {
    E_k = 0.0;
    for (int i = 0; i < N_LOCAL; i++) {
        E_k += px[i] * px[i] + py[i] * py[i] + pz[i] * pz[i];
    }
    E_k = std::sqrt(E_k);
    E_k = E_k / m / (2 * force_conversion_factor);

    T = E_k / ((3 * N_LOCAL - 3) * R_const);
}

void Simu::normalize_momentums() {
    double ratio = (3 * N_LOCAL - 3) * R_const * T_0 / E_k;
    for (int i = 0; i < N_LOCAL; i++) {
        px[i] *= ratio;
        py[i] *= ratio;
        pz[i] *= ratio;
    }
}

void Simu::print() {
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "iter: " << ticks << ", total energy: " << U << ", kinetic energy: " << E_k << ", temp: " << T <<std::endl;
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
    std::ofstream f; 
    f.open("saved/" + std::to_string(ticks) + ".data", std::ios::out | std::ios::binary);

    double n = N_LOCAL;
    double l = L;
    double iter = ticks;

    f.write(reinterpret_cast<const char *>(&n), sizeof(double));
    f.write(reinterpret_cast<const char *>(&l), sizeof(double));

    f.write(reinterpret_cast<const char *>(&iter), sizeof(double));
    f.write(reinterpret_cast<const char *>(&U), sizeof(double));
    f.write(reinterpret_cast<const char *>(&T), sizeof(double));
    f.write(reinterpret_cast<const char *>(&E_k), sizeof(double));

    for (int i = 0; i < N_LOCAL; i++) {
        f.write(reinterpret_cast<const char *>(&x.data()[i]), sizeof(double));
        f.write(reinterpret_cast<const char *>(&y.data()[i]), sizeof(double));
        f.write(reinterpret_cast<const char *>(&z.data()[i]), sizeof(double));
    }


    f.close();
}
