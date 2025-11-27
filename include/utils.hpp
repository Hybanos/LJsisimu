#include <array>
#include <fstream>

static std::array<double, N_TOTAL * 3> load_data() {
    std::array<double, N_TOTAL * 3> data;

    char *tmp = new char[N_TOTAL*3*sizeof(double)];

    std::fstream f;
    f.open("particule.data", std::ios::in | std::ios::binary);
    f.read(tmp, N_TOTAL * 3 * sizeof(double));
    f.close();

    for (int i = 0; i < N_TOTAL * 3; i++) {
        data[i] = ((double *)tmp)[i];
    }

    delete tmp;

    return data;
}