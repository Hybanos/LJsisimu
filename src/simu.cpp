#include "simu.hpp"

Simu::Simu() {
    auto data = load_data();

    for (int i = 0; i < N_LOCAL; i++) {
        x[i] = data[i * 3]    ;
        y[i] = data[i * 3 + 1];
        z[i] = data[i * 3 + 2];

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

    int len = std::pow(N_SYM, 1.0/3.0);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            for (int k = 0; k < len; k++) {
                int img = i * 9 + j * 3 + k;
                int offset = len / 3;
                imgs[img] = vec3((i-offset)*L, (j-offset)*L, (k-offset)*L);
                vec3((i-offset)*L, (j-offset)*L, (k-offset)*L).print();
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

    // apply forces
    for (int i = 0; i < N_LOCAL; i++) {
        // Velocity-verlet
        double x_tp1 = x[i] + px[i] / m * timestep + fx[i] * timestep * timestep * 0.5 * force_conversion_factor;
        double y_tp1 = y[i] + py[i] / m * timestep + fy[i] * timestep * timestep * 0.5 * force_conversion_factor;
        double z_tp1 = z[i] + pz[i] / m * timestep + fz[i] * timestep * timestep * 0.5 * force_conversion_factor;

        double vx_tp1 = px[i] / m + (fx[i] + fx_tmp[i]) * timestep * 0.5 * force_conversion_factor * force_conversion_factor;
        double vy_tp1 = py[i] / m + (fy[i] + fy_tmp[i]) * timestep * 0.5 * force_conversion_factor * force_conversion_factor;
        double vz_tp1 = pz[i] / m + (fz[i] + fz_tmp[i]) * timestep * 0.5 * force_conversion_factor * force_conversion_factor;

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
}

// double Simu::dist_squared(int i, int j) {
//     return std::pow(x[i] - x[j], 2) + std::pow(y[i] - y[j], 2) + std::pow(z[i] - z[j], 2);
// }

void Simu::print() {
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "iter: " << ticks << ", total energy: " << U << std::endl;
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
