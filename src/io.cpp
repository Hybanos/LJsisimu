#include "simu.hpp"

void Simu::load() {
    std::array<double, N_TOTAL * 3> data;

    char *tmp = new char[N_TOTAL*3*sizeof(double)];

    std::fstream f;
    f.open("exam_xyz", std::ios::in | std::ios::binary);
    f.read(tmp, N_TOTAL * 3 * sizeof(double));
    f.close();

    for (int i = 0; i < N_TOTAL; i++) {
        x[i] = ((double *)tmp)[i * 3 + 0];
        y[i] = ((double *)tmp)[i * 3 + 1];
        z[i] = ((double *)tmp)[i * 3 + 2];
    }

    f.open("exam_mci", std::ios::in | std::ios::binary);
    f.read(tmp, N_TOTAL * 3 * sizeof(double));
    f.close();

    for (int i = 0; i < N_TOTAL; i++) {
        px[i] = ((double *)tmp)[i * 3 + 0];
        py[i] = ((double *)tmp)[i * 3 + 1];
        pz[i] = ((double *)tmp)[i * 3 + 2];
    }

    delete tmp;
}

void Simu::print() {
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "iter: " << steps << ", sps: " << sps << std::endl;
    std::cout << "M_rc: " << M_rc[0] << ", " << M_rc[1] << ", " << M_rc[2] << ", " << M_rc[3] << ", " << std::endl;
    std::cout << "total energy: " << U + E_k << ", U: " << U << ", E_k: " << E_k << ", temp: " << T <<std::endl;
    double xx = 0, yy = 0, zz = 0;
    double vx = 0, vy = 0, vz = 0;
    for (int i = 0; i < N_LOCAL; i++) {
        xx += fx[i];
        yy += fy[i];
        zz += fz[i];
    }
    for (int i = 0; i < N_LOCAL; i++) {
        vx += px[i] / m;
        vy += py[i] / m;
        vz += pz[i] / m;
    }
    std::cout << "f total: " << xx << ", " << yy << ", " << zz << std::endl;
    std::cout << "v total: " << vx << ", " << vy << ", " << vz << std::endl;
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

void Simu::print_dist() {
    std::cout << std::endl;
    int i = 0;
    for (auto a : temps) {
        std::cout << i * temps_step + temps_min << "\t"; 
        std::cout << a << std::endl;
        i++;
    }
}