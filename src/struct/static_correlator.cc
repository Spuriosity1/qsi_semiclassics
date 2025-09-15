/* Implementation of class static_correlator.
 *
 * Created on 31/07/2018
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

#include "static_correlator.hh"
#include "corrtypes.hh"
#include <misc.hh>
#include <complex>
#include <fftw3.h>
#include "fourier.hh"
#include <cstdio>
#include <cstring>
#include <string>

using namespace std;
using namespace corrtypes;

static_correlator::static_correlator (size_t n, unsigned what,
                                      const fourier* f):
    what(what),
    n(n),
    n_corr(0),
    m_fourier(*f)
{
    // Creating storage arrays, all of size 2n^2
    // Zeroing out correlator arrays
    for (int i = 0; i < 4; ++i) {
        if (what & SZZ) 
            S_z[i] = new complex<double>[2*n*n];
        if (what & SXX) 
            S_x[i] = new complex<double>[2*n*n];
        if (what & SYY) 
            S_y[i] = new complex<double>[2*n*n];
        if (what & SPM)
            S_p[i] = new complex<double>[2*n*n];
        if (what & IMAG) 
            S_b[0][i] = new complex<double>[2*n*n];
        if (what & ARG)
            S_b[1][i] = new complex<double>[2*n*n];
        if (what & SIN)
            S_b[2][i] = new complex<double>[2*n*n];
        if (what & COS)
            S_b[3][i] = new complex<double>[2*n*n];
        
        for (int j = 0; j < 4; ++j) {
            if (what & SZZ) {
                S_zz[i][j] = new complex<double>[2*n*n];
                std::fill_n(S_zz[i][j], 2*n*n, complex<double>(0));
            }

            if (what & SXX) {
                S_xx[i][j] = new complex<double>[2*n*n];
                std::fill_n(S_xx[i][j], 2*n*n, complex<double>(0));
            }

            if (what & SYY) {
                S_yy[i][j] = new complex<double>[2*n*n];
                std::fill_n(S_yy[i][j], 2*n*n, complex<double>(0));
            }
            
            if (what & SPM) {
                S_pm[i][j] = new complex<double>[2*n*n];
                std::fill_n(S_pm[i][j], 2*n*n, complex<double>(0));
            }

            if (what & IMAG) {
                S_bb[0][i][j] = new complex<double>[2*n*n];
                std::fill_n(S_bb[0][i][j], 2*n*n, complex<double>(0));
            }

            if (what & ARG) {
                S_bb[1][i][j] = new complex<double>[2*n*n];
                std::fill_n(S_bb[1][i][j], 2*n*n, complex<double>(0));
            }

            if (what & SIN) {
                S_bb[2][i][j] = new complex<double>[2*n*n];
                std::fill_n(S_bb[2][i][j], 2*n*n, complex<double>(0));
            }

            if (what & COS) {
                S_bb[3][i][j] = new complex<double>[2*n*n];
                std::fill_n(S_bb[3][i][j], 2*n*n, complex<double>(0));
            }
        }
    }
}

static_correlator::~static_correlator() {
    // Removing various arrays
    for (int i = 0; i < 4; ++i) {
        if (what & SZZ) 
            delete[] S_z[i];
        if (what & SXX) 
            delete[] S_x[i];
        if (what & SYY) 
            delete[] S_y[i];
        if (what & SPM)
            delete[] S_p[i];
        if (what & IMAG)
            delete[] S_b[0][i];
        if (what & ARG)
            delete[] S_b[1][i];
        if (what & SIN)
            delete[] S_b[2][i];
        if (what & COS)
            delete[] S_b[3][i];
        
        for (int j = 0; j < 4; ++j) {
            if (what & SZZ)
                delete[] S_zz[i][j];
            if (what & SXX) 
                delete[] S_xx[i][j];
            if (what & SYY) 
                delete[] S_yy[i][j];
            if (what & SPM)
                delete[] S_pm[i][j];
            if (what & IMAG)
                delete[] S_bb[0][i][j];
            if (what & ARG)
                delete[] S_bb[1][i][j];
            if (what & SIN)
                delete[] S_bb[2][i][j];
            if (what & COS)
                delete[] S_bb[3][i][j];
        }
    }
}

// Indexing other arrays in a row-major order
size_t static_correlator::index(int x, int z) const {
    // Sanitise input, PBC
    x = mod(x, 2*n);
    z = mod(z, n);
    return x*n + z;
}

void static_correlator::set(int sublat, unsigned type) {
    // Always stands at the front of the spin array line being filled out
    complex<double>* run = NULL;
    // Find which correlator has been requested and if it is stored here
    if (type & IMAG)
        if (what & IMAG)
            run = S_b[0][sublat];
        else return;
    else if (type & ARG)
        if (what & ARG)
            run = S_b[1][sublat];
        else return;
    else if (type & SIN)
        if (what & SIN)
            run = S_b[2][sublat];
        else return;
    else if (type & COS)
        if (what & COS)
            run = S_b[3][sublat];
        else return;
    else if (type & SZZ)
        if (what & SZZ)
            run = S_z[sublat];
        else return;
    else if (type & SXX)
        if (what & SXX)
            run = S_x[sublat];
        else return;
    else if (type & SYY)
        if (what & SYY)
            run = S_y[sublat];
        else return;
    else if (type & SPM)
        if (what & SPM)
            run = S_p[sublat];
        else return;
    else return;

    // Do the actual copying
    for(size_t x = 0; x < 2*n; ++x) {
        memcpy(run, &m_fourier(x,x,0), n*sizeof(complex<double>));
        run += n;
    }
}
    
/* We simply need to calculate X_i(q) [X_j(q)]*
 *   for all q,i,j, and types of stored observables X,
 * and add them to whatever already is in correlator */
void static_correlator::add_correlator() {
    for (size_t q = 0; q < 2*n*n; ++q) // k-space point enumerations
        for (int i = 0; i < 4; ++i) // source sublattice
            for (int j = 0; j < 4; ++j) { // dest sublattice
                if (what & SZZ)
                    S_zz[i][j][q] += S_z[i][q] * conj(S_z[j][q]);
                if (what & SXX)
                    S_xx[i][j][q] += S_x[i][q] * conj(S_x[j][q]);
                if (what & SYY)
                    S_yy[i][j][q] += S_y[i][q] * conj(S_y[j][q]);
                if (what & SPM)
                    S_pm[i][j][q] += S_p[i][q] * conj(S_p[j][q]);
                if (what & IMAG)
                    S_bb[0][i][j][q] += S_b[0][i][q] * conj(S_b[0][j][q]);
                if (what & ARG)
                    S_bb[1][i][j][q] += S_b[1][i][q] * conj(S_b[1][j][q]);
                if (what & SIN)
                    S_bb[2][i][j][q] += S_b[2][i][q] * conj(S_b[2][j][q]);
                if (what & COS)
                    S_bb[3][i][j][q] += S_b[3][i][q] * conj(S_b[3][j][q]);
            }
    ++n_corr;
}

// Filename extensions
const string x_szz = ".szz.stat";
const string x_sxx = ".sxx.stat";
const string x_syy = ".syy.stat";
const string x_spm = ".spm.stat";
const string x_imag = ".imag.stat";
const string x_arg  = ".arg.stat";
const string x_sin = ".sin.stat";
const string x_cos = ".cos.stat";

void static_correlator::save_corr(const char* root) const {
    // array in which to divide total correlators by n_corr * (#spins)
    complex<double>* div = new complex<double>[2*n*n];
    double rec = 1.0 / n_corr / (16*n*n*n);

    if (what & SZZ) {
        FILE* f = fopen((root+x_szz).c_str(), "wb");
        // for each array, divide each q-point by n_corr*(#spins), then output
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                for (size_t q = 0; q < 2*n*n; ++q)
                    div[q] = S_zz[i][j][q] * rec;
                fwrite(div, sizeof(complex<double>), 2*n*n, f);
            }
        fclose(f);
        // Final binary structure: 4x4 blocks of 2n^2 
        // [0][0] blocks give S^z_0 S^z_0 corr
    }

    if ( what & SXX ){
        FILE* f = fopen((root+x_sxx).c_str(), "wb");
        // for each array, divide each q-point by n_corr*(#spins), then output
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                for (size_t q = 0; q < 2*n*n; ++q)
                    div[q] = S_xx[i][j][q] * rec;
                fwrite(div, sizeof(complex<double>), 2*n*n, f);
            }
        fclose(f);
    }

    if ( what & SYY ){
        FILE* f = fopen((root+x_syy).c_str(), "wb");
        // for each array, divide each q-point by n_corr*(#spins), then output
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                for (size_t q = 0; q < 2*n*n; ++q)
                    div[q] = S_yy[i][j][q] * rec;
                fwrite(div, sizeof(complex<double>), 2*n*n, f);
            }
        fclose(f);
    }
    
    // same for +/-
    if (what & SPM) {
        FILE* f = fopen((root+x_spm).c_str(), "wb");
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                for (size_t q = 0; q < 2*n*n; ++q)
                    div[q] = S_pm[i][j][q] * rec;
                fwrite(div, sizeof(complex<double>), 2*n*n, f);
            }
        // close output file
        fclose(f);
    }
    
    // same for b
    if (what & IMAG) {
        FILE* f = fopen((root+x_imag).c_str(), "wb");
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                for (size_t q = 0; q < 2*n*n; ++q)
                    div[q] = S_bb[0][i][j][q] * rec;
                fwrite(div, sizeof(complex<double>), 2*n*n, f);
            }
        fclose(f);
    }
    
    // same for B
    if (what & ARG) {
        FILE* f = fopen((root+x_arg).c_str(), "wb");
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                for (size_t q = 0; q < 2*n*n; ++q)
                    div[q] = S_bb[1][i][j][q] * rec;
                fwrite(div, sizeof(complex<double>), 2*n*n, f);
            }
        fclose(f);
    }
    
    // same for sin B
    if (what & SIN) {
        FILE* f = fopen((root+x_sin).c_str(), "wb");
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                for (size_t q = 0; q < 2*n*n; ++q)
                    div[q] = S_bb[2][i][j][q] * rec;
                fwrite(div, sizeof(complex<double>), 2*n*n, f);
            }
        fclose(f);
    }
       
    // same for cos B
    if (what & COS) {
        FILE* f = fopen((root+x_cos).c_str(), "wb");
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j){
                for (size_t q = 0; q < 2*n*n; ++q)
                    div[q] = S_bb[3][i][j][q] * rec;
                fwrite(div, sizeof(complex<double>), 2*n*n, f);
            }
        fclose(f);
    }
    
    // delete tmp array
    delete[] div;
}
