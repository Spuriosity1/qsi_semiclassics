/* Executable to draw SF neutron scattering pattern
 *
 * Usage: spin_plot "system size" "input filename" "output filename"
 *
 * Created on 01/08/2018
 * Copyright (C) 2018 Attila Szabó <as2372@cam.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Creative Commons Attribution License (CC-BY),
 * version 4.0 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Please find a copy of the Creative Commons CC-BY License on
 * <https://creativecommons.org/licenses/by/4.0/>.
 */

#include <spin_manager.hh>
#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char** argv) {
    int n = atoi(argv[1]);
    spin_manager corr(n);

    FILE* in = fopen(argv[2], "rb");
    corr.load(in);
    fclose(in);

    double* sf = new double[(4*n+1)*(6*n+1)];
    int i = 0;

    for (int x = -2*n; x <= 2*n; ++x)
        for (int z = -3*n; z <= 3*n; ++z) {
            corr.tensor(x,z,spin_manager::ISING);
            sf[i] = corr.neutron_sf();
            ++i;
        }

    FILE* out = fopen(argv[3], "wb");
    fwrite(sf, sizeof(double), (4*n+1)*(6*n+1), out);
    fclose(out);

    delete[] sf;

    return 0;
}
