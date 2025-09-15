/* Executable to find mean frequency and total structure factor
 *
 * Input file is assumed to contain "lines" of a fixed number of double
 * precision floats. These lines are summed and their first moments taken to
 * find the total structure factor and mean frequency, resp. These are printed
 * in the output file.
 *
 * Usage: mean_freq line_length infile outfile
 *
 * Created on 17/08/2018
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

#include <cstdio>
#include <cstdlib>

using namespace std;

int main (int argc, const char** argv) {
    int L = atoi(argv[1]); // Line length
    double* slurp = new double[L];

    FILE* f = fopen(argv[2],"rb");
    fseek(f, 0L, SEEK_END);
    size_t N = ftell(f); // length of file
    fprintf(stderr, "Found %lu bytes %lu doubles, %lu lines", N, N/sizeof(double), N/sizeof(double)/L);
    if (N % (3*L*sizeof(double)))
        throw std::runtime_error("File cannot be divided into lines");
    N /= (3*L*sizeof(double));
    rewind(f);

    for (int run = 0; run < 3; ++run) {
        for (unsigned i = 0; i < N; ++i) {
            fread(slurp,sizeof(double),L,f);
            double sum = 0, mom = 0;
            for (int w = 0; w < L; ++w) {
                sum += slurp[w];
                mom += slurp[w] * w;
            }
            // sanitise very weak signals at Gamma points
            if (sum < 1e-9)
                mom = 0;
            printf("%3u %10.6f %.6e\n",i,mom/sum,sum);
        }
        printf("\n");
    }

    fclose(f);
    delete[] slurp;
    return 0;
}
