/* Implementation of class fourier
 *
 * Created on 03/08/2018
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

#include "fourier.hh"
#include <vec3.hh>
#include <misc.hh>
#include <complex>
#include <fftw3.h>
#include <cstring>

using namespace std;

fourier::fourier (size_t N, bool direction):
    n(2*N)
{
    // Creating Fourier transform array using FFTW methods
    // and casting it onto complex<double>
    fft_w = fftw_alloc_complex(n*n*n);
    fft = (complex<double>*) fft_w;

    // Creating FFTW plan
    plan = fftw_plan_dft_3d(n, n, n, fft_w, fft_w,
                            (direction ? FFTW_FORWARD : FFTW_BACKWARD),
                            FFTW_MEASURE);
}

fourier::~fourier() {
    // Removing Fourier transform array
    // NB complex<double> pointer is just a typecast, no need to delete[] it
    fftw_free(fft_w);

    // Removing FFTW plan
    fftw_destroy_plan(plan);
}


// Indexing fft array in a row-major order
size_t fourier::index(int x, int y, int z) const {
    // Sanitise input, PBC
    x = mod(x, n);
    y = mod(y, n);
    z = mod(z, n);
    return (x*n + y)*n + z;
}

void fourier::clear_input() {
    // memset(fft, 0, n*n*n*sizeof(complex<double>));
    std::fill_n(fft, n*n*n, complex<double>(0));
}

void fourier::transform() {
    fftw_execute(plan);
}
