/* Executable to draw vison charge correlators
 *
 * Usage: vison_plot "system size" "input filename" "output filename"  [block]
 **   block specifies which block of infile is to be used (default: 0)
 *
 * Created on 11/09/2018
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

#include <magnetic_manager.hh>
#include <cstdio>
#include <cstdlib>

using namespace std;

int main(int argc, char** argv) {
    unsigned block = 0u;
    if (argc > 4)
        block = atoi(argv[4]);
    
    int n = atoi(argv[1]);
    magnetic_manager corr(n);

    FILE* in = fopen(argv[2], "rb");
    corr.load(in, block);
    fclose(in);

    double* qq = new double[(4*n+1)*(6*n+1)];
    int i = 0;
    
    for (int x = -2*n; x <= 2*n; ++x)
        for (int z = -3*n; z <= 3*n; ++z) {
            qq[i] = corr.vison(x,z);
            ++i;
        }

    FILE* out = fopen(argv[3], "wb");
    fwrite(qq, sizeof(double), (4*n+1)*(6*n+1), out);
    fclose(out);

    delete[] qq;

    return 0;
}
