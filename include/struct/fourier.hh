/* Class fourier
 * Collects spin layouts and performs real space Fourier transforms on them 
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

#ifndef fourier_hh
#define fourier_hh

#include <vec3.hh>
#include <complex>
#include <fftw3.h>

class fourier {
private:
    // Number of q-points along each axis, twice the number of unit cells
    const size_t n;

    // FFTW workhorse array
    std::complex<double>* fft;
    fftw_complex* fft_w;

    // FFTW plan
    fftw_plan plan;

    // Indexing of full-size arrays in a row-major order as required by FFTW
    size_t index(int x, int y, int z) const;

public:
    /* Constructor
    **  n: number of unit cells per direction
    **  direction: true for forward FFT, false for backward FFT */
    fourier(size_t n, bool direction = true);
    // Destructor
    ~fourier();

    // Set real space input, indexing by dimension
    inline std::complex<double>& operator()(int x, int y, int z) {
        return fft[index(x,y,z)];
    }
    inline std::complex<double>& operator()(const vec3_int& v) {
        return fft[index(v[0],v[1],v[2])];
    }
    // Retrieve contents of fft array, indexing by dimension
    inline const
    std::complex<double>& operator()(int x, int y, int z) const {
        return fft[index(x,y,z)];
    }
    inline const
    std::complex<double>& operator()(const vec3_int& v) const {
        return fft[index(v[0],v[1],v[2])];
    }

    // Zero out input array
    void clear_input();

    // Perform FT on contents of fft
    void transform();
};

#endif
