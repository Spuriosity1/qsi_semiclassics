/* Takes the Fourier transform of the raw spin output data
 * Not sure what the driver program is here...
 *
 * Usage: autocorr_process L N infile > outfile
 **    L:       linear size of sample
 **    N:       number of samples
 **    infile:  binary file containing values of S^z
 **    outfile: text file listing of autocorrelation in time and freq. domain
 *
 * Created on 16/09/2018
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

#include <complex>
#include <fftw3.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

int main (int argc, char** argv) {
    size_t L = atoi(argv[1]);
    // size of a full spin configuration in bytes
    size_t line = 16ul*L*L*L * sizeof(double);
    // proceed by blocks of 128 spins (some may be lost at the end if L odd)
    size_t Nb = L*L*L/8;
    size_t block = 128 * sizeof(double);

    int N = atoi(argv[2]);
    if (N%2)
        throw "FFTW requires an even number of samples.";

    double* input  = fftw_alloc_real(128ul * N);
    fftw_complex* output_w = fftw_alloc_complex(128ul * (N/2+1));
    complex<double>* output = (complex<double>*) output_w;
    fftw_plan plan = fftw_plan_many_dft_r2c(1, &N, 128,
                                            input, NULL, 128, 1,
                                            output_w, NULL, 128, 1,
                                            FFTW_MEASURE);
    
    double* freq = fftw_alloc_real(N/2+1);
    double* time = fftw_alloc_real(N/2+1);
    fftw_plan result = fftw_plan_r2r_1d(N/2+1, freq, time, FFTW_REDFT00,
                                        FFTW_ESTIMATE);
    memset(freq,0,(N/2+1)*sizeof(double));
    
    FILE* f = fopen(argv[3], "rb");

    for (size_t i = 0; i < Nb; ++i) {
        for (int j = 0; j < N; ++j) {
            fseek(f, j*line+i*block, SEEK_SET);
            fread(input+128*j, sizeof(double), 128, f);
        }
        fftw_execute(plan);
        for (int j = 0; j < N/2+1; ++j)
            for (int k = 0; k < 128; ++k)
                freq[j] += norm(output[128*j+k]);
    }

    fclose(f);
    
    fftw_execute(result);

    for (int i = 0; i < N/2+1; ++i)
        printf("%.14e %.14e\n", freq[i], time[i]);

    fftw_free(input);
    fftw_free(output_w);
    fftw_free(freq);
    fftw_free(time);
    fftw_destroy_plan(plan);
    fftw_destroy_plan(result);

    return 0;
}
