/* Class dynamic_correlator
 * Dynamic spin correlators for pyrochlore Heisenberg models,
 * restricted to high-symmetry axes
 *
 * Created on 03/08/2018
 * Copyright (C) 2018,2019 Attila Szabó <as2372@cam.ac.uk>
 *cat .
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

#ifndef dynamic_correlator_hh
#define dynamic_correlator_hh

#include <complex>
#include <fftw3.h>
#include <vec3.hh>
#include <vector>
#include <cstdio>

class fourier;

class dynamic_correlator {
public:
    // Indicator of what correlators are kept
    const unsigned what;

private:
    // Number of unit cells along each axis 
    const size_t n;

    // Default unit of wave vector, 2pi/n
    const double Q;

    // Number of time steps performed
    const size_t n_time;

    // Number of frequencies saved in output
    const size_t n_freq;

    // Number of correlators stored
    size_t n_corr;

    // std::vector storing the path for correlators and for noise
    std::vector<vec3_int> kspace_path;
    std::vector<vec3_int> path_n;

    // Real space Fourier transform generator
    const fourier& m_fourier;

    /* ---------- Arrays for Fourier transforms & correlators ----------
     * S^z(q,t), S^+(q,t), B(q,t), b(q,t)
     * are stored in one contiguous array, in the above order,
     * each split by sublattice, each being q-major order.
     *
     * B = arg(S^+S^-...) and b = Im(S^+S^-...) ~ sin B are proxies for
     * Hermele's magnetic field B
     * 
     * The whole array is equipped with an FFTW plan to generate ***(q,omega)
     * which are then mod-squared and summed into correlator arrays */

    // Big array, big FFTW plan
    fftw_complex* corr_data = NULL;
    fftw_plan corr_plan = NULL;

    // Aux pointers for FTs of various components
    std::complex<double>* S_x[4];
    std::complex<double>* S_y[4];
    std::complex<double>* S_z[4];
    std::complex<double>* S_p[4];
    std::complex<double>* S_B[4];
    std::complex<double>* S_b[4];
    std::complex<double>* S_sin[4];
    std::complex<double>* S_cos[4];

    // Correlator arrays
    double* S_xx = NULL;
    double* S_yy = NULL;
    double* S_zz = NULL;
    double* S_pm = NULL;
    
    double* S_BB = NULL;
    double* S_bb = NULL;
    double* S_sinsin = NULL;
    double* S_coscos = NULL;

    /* ---------- Arrays for noise data & correlators ----------
     * Svec(q,t) is stored along the high-symmetry directions
     *     [100] (h = 1..2n), [110] (h = 1..2n), [111] (h = 1..n)
     * in one contiguous array, indexed as [t/omega][q][alpha]
     *
     * Svec(q,t) = sum_mu S_mu^z(q,t) e_mu is the Fourier transform of the
     * physical magnetisation of the system
     *
     * Each Svec(q,omega) thus forms a vec3_cplx in memory, which will also be
     * referenced as such.
     *
     * The whole array is equipped with an FFTW plan to generate Svec(q,omega)
     * which are mod-squared and summed using the appropriate kernel into
     * a noise correlator vector (index = frequency) */

    // Big array, big FFTW plan
    fftw_complex* noise_data = NULL;
    fftw_plan noise_plan = NULL;

    // Typecast into vec3_cplx
    vec3_cplx* noise_vec;

    // Noise correlator array [3 x n_freq]
    double* noise = NULL;

    /* ---------- Various indexing routines ----------
     * To facilitate FFTW Fourier transforms,
     ** spin/field arrays are q-major,
     ** spin vector arrays (for noise) are (t/omega)-major
     * For simplicity, all correlator arrays are omega-minor */

    // Indexing spin component arrays
    inline size_t index_spin(size_t q, size_t t) const {
        return q * n_time + t;
    }
    
    // Indexing spin vector arrays
    inline size_t index_noise(size_t q, size_t t) const {
        return t * path_n.size() + q;
    }

    // Indexing correlator arrays
    inline size_t index_corr(size_t q, size_t w) const {
        return q * n_freq + w;
    }

    // Generating the path for n = 4m 
    void path4();
    // Generating alternative path for n = 2m (to be used if n%4==2)
    void path2();
    // Generating the path for noise estimation
    void path_noise();
    
public:
    // Constructor
    dynamic_correlator(size_t n, size_t n_time, size_t n_freq, unsigned what,
                       const fourier* f);
    // Destructor
    ~dynamic_correlator();

    /* Set spin data at time step t from m_fourier
    **   sublat: 0..3, index of pyrochlore sublattice
    **   what:   a single bit specifying what to store, if multiple given,
    **           the smallest one is used
    **           NOISE is invalid, as that is calculated from Sz correlators */
    void set (size_t t, int sublat, unsigned what);

    /* Clear noise input arrays
     * has to be called before starting to save any new time evolution */
    void clear_input();

    // Add current spin set-up to correlator tally
    void add_correlator();

    /* Save correlators into binary file
     * Format: one file for each type, containing single double array
     * Indexing (row major order): [q/dir][omega] */
    void save_corr(const char* root) const;
};

#endif
