/* Class static_correlator
 * Static spin correlators for pyrochlore Heisenberg models, 
 * restricted to the (hhk) plane.
 *
 * Created on 30/07/2018
 * Copyright (C) 2018,2019 Attila Szabó <as2372@cam.ac.uk>
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

#ifndef static_correlator_hh
#define static_correlator_hh

#include <complex>
#include <fftw3.h>
#include <vec3.hh>

class fourier;

class static_correlator {
public:
    // Indicator of what correlators are handled
    const unsigned what;

private:
    enum class bb_correlator_type{
        COS = 3,
        SIN = 2,
        ARG = 1,
        IMAG = 0
    };
    // Number of unit cells along each axis
    const size_t n;
    
    // Number of correlators stored
    size_t n_corr;

    // Real space Fourier transform generator
    const fourier& m_fourier;
    

    /* ---------- Arrays of Fourier transforms ----------
     * Since we look at functions defined on an fcc lattice, using FFTW naively
     * creates two Brillouin zones, with all q-point included.
     * We are only interested in the (hhk) plane, and also, there is no need to
     * store but half of each array.
     * This is accomplished by keeping only lines whose q_x and q_y coordinates
     * are equal and cutting them in half along the q_z axis.
     * The stored q-points then come from one (rather unconventional) reciprocal
     * space unit cell.
     *
     * Furthermore, the only correlators that are non-zero under our Hamiltonian
     * are <S^zS^z> and <S^+S^->. Accordingly, only these spin FT's and
     * correlators are stored, as well 3 proxies for Hermele's magnetic field:
     **  <bb> where b = Im(S^+S^-...)  (in S_b[0]);
     **  <BB> where B = arg(S^+S^-...) (in S_b[1]);
     **  <sin(B) sin(B)>               (in S_b[2]). */

    // Arrays to keep the spin FT, "reduced size", indexed by sublattice
    std::complex<double>* S_x[4];
    std::complex<double>* S_y[4];
    std::complex<double>* S_z[4];
    std::complex<double>* S_p[4];
    // First index: 0 = IMAG 1 = ARG 2 = SIN 3 = COS
    // TODO change this to an enum class
    std::complex<double>* S_b[4][4];

    // Arrays to collect spin correlators, "reduced size" 
    std::complex<double>* S_zz[4][4];
    std::complex<double>* S_xx[4][4];
    std::complex<double>* S_yy[4][4];
    std::complex<double>* S_pm[4][4];
    std::complex<double>* S_bb[4][4][4];
    // First index: 0 = IMAG 1 = ARG 2 = SIN 3 = COS
    // TODO change this to an enum class

    // Indexing reduced-size arrays in a row-major order
    size_t index(int x, int z) const;

public:
    // Constructor
    static_correlator(size_t n, unsigned what, const fourier* f);
    // Destructor
    ~static_correlator();

    /* Set spin data from m_fourier
    **   sublat: 0..3, index of pyrochlore sublattice
    **   what:   a single bit specifying what to store, if multiple given,
    **           the smallest one is used */
    void set (int sublat, unsigned what);

    // Add current spin set-up to correlator tally
    void add_correlator();

    /* Save correlators into binary files
     * Format: one file for each type
     **    B-field related: single multidimensional complex<double> array
     * Indexing (row major order): [i][j][qx][qz][re/im] */
    void save_corr(const char* root) const;
};

#endif
