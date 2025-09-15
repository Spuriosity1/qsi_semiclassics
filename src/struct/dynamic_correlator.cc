/* Implementation of class dynamic_correlator
 *
 * Created on 03/08/2018
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

#include "dynamic_correlator.hh"
#include "corrtypes.hh"
#include "fourier.hh"
#include <complex>
#include <fftw3.h>
#include <vec3.hh>
#include <direction.hh>
#include <cstdio>
#include <vector>
#include <cstring>

using namespace std;
using namespace corrtypes;

#define TWOPI 6.283185307179586476925286766559005768394338798750211641949
#define PI 3.141592653589793238462643383279502884197169399375105820974

dynamic_correlator::dynamic_correlator (size_t n, size_t n_time, size_t n_freq,
                                        unsigned what, const fourier* f) :
    what(what),
    n(n),
    Q(TWOPI/n),
    n_time(n_time),
    n_freq(n_freq),
    n_corr(0),
    m_fourier(*f)
{
    // ---------- CONSTRUCTION FOR SPIN/FIELD CORRELATORS ----------
    if (what & ALL) {
        if (!(n%4)) // 4|n, use path4
            path4();
        else if (!(n%2)) // n=4m+2, use path2
            path2();
        else // n is odd, we're not prepared for this
            throw std::runtime_error("Can only take dynamic correlators for even sized boxes");

        // Count how many correlator arrays are needed
        int n_corrtypes = 0;
    
        // Correlator storage arrays, zeroed out
        if (what & SXX) {
            S_xx = new double[kspace_path.size()*n_freq];
            memset(S_xx, 0, kspace_path.size()*n_freq*sizeof(double));
            ++n_corrtypes;
        }
        if (what & SYY) {
            S_yy = new double[kspace_path.size()*n_freq];
            memset(S_yy, 0, kspace_path.size()*n_freq*sizeof(double));
            ++n_corrtypes;
        }
        if (what & SZZ) {
            S_zz = new double[kspace_path.size()*n_freq];
            memset(S_zz, 0, kspace_path.size()*n_freq*sizeof(double));
            ++n_corrtypes;
        }
        if (what & SPM) {
            S_pm = new double[kspace_path.size()*n_freq];
            memset(S_pm, 0, kspace_path.size()*n_freq*sizeof(double));
            ++n_corrtypes;
        }
        if (what & ARG) {
            S_BB = new double[kspace_path.size()*n_freq];
            memset(S_BB, 0, kspace_path.size()*n_freq*sizeof(double));
            ++n_corrtypes;
        }
        if (what & IMAG) {
            S_bb = new double[kspace_path.size()*n_freq];
            memset(S_bb, 0, kspace_path.size()*n_freq*sizeof(double));
            ++n_corrtypes;
        }
        if (what & SIN){
            S_sinsin = new double[kspace_path.size()*n_freq];
            memset(S_sinsin, 0, kspace_path.size()*n_freq*sizeof(double));
            ++n_corrtypes;
        }
        if (what & COS){
            S_coscos = new double[kspace_path.size()*n_freq];
            memset(S_coscos, 0, kspace_path.size()*n_freq*sizeof(double));
            ++n_corrtypes;
        }
        
        // FFTW array 
        corr_data = fftw_alloc_complex(n_corrtypes*4*kspace_path.size()*n_time);
    
        // FFTW plan: 1d forward FT in multiple adjacent arrays
        int time = n_time; // ffs
        corr_plan = fftw_plan_many_dft
            (1, &time, n_corrtypes*4*kspace_path.size(), //rank, n, howmany
             corr_data, NULL, 1, time, //in,  embed, stride, dist
             corr_data, NULL, 1, time, //out, embed, stride, dist
             FFTW_FORWARD, FFTW_MEASURE); //sign, flags
        
        // Distribute to components
        complex<double>* run = (complex<double>*) corr_data;
        if (what & SXX)
            for (int i = 0; i < 4; ++i) {
                S_x[i] = run;
                run += kspace_path.size()*n_time;
            }
        if (what & SYY)
            for (int i = 0; i < 4; ++i) {
                S_y[i] = run;
                run += kspace_path.size()*n_time;
            }
        if (what & SZZ)
            for (int i = 0; i < 4; ++i) {
                S_z[i] = run;
                run += kspace_path.size()*n_time;
            }
        if (what & SPM)
            for (int i = 0; i < 4; ++i) {
                S_p[i] = run;
                run += kspace_path.size()*n_time;
            }
        if (what & ARG)
            for (int i = 0; i < 4; ++i) {
                S_B[i] = run;
                run += kspace_path.size()*n_time;
            }
        if (what & IMAG)
            for (int i = 0; i < 4; ++i) {
                S_b[i] = run;
                run += kspace_path.size()*n_time;
            }
        if (what & SIN)
            for (int i = 0; i < 4; ++i) {
                S_sin[i] = run;
                run += kspace_path.size()*n_time;
            }
        if (what & COS)
            for (int i = 0; i < 4; ++i) {
                S_cos[i] = run;
                run += kspace_path.size()*n_time;
            }
    } // END SPIN/FIELD CORRELATORS

    // ---------- CONSTRUCTION FOR NOISE CORRELATORS ----------
    if (what & NOISE) {
        path_noise();

        // Correlator storage array; 3* due to 3 high-symmetry directions
        noise = new double[3*n_freq]; 
        memset(noise, 0, 3*n_freq*sizeof(double));

        // FFTW array; 3* due to 3 vector components
        noise_data = fftw_alloc_complex(3*path_n.size()*n_time);
        noise_vec = (vec3_cplx*) noise_data;

        // FFTW plan: 1d forward FT in multiple interwoven arrays
        int time = n_time; // ffs
        int howmany = 3 * path_n.size();
        noise_plan = fftw_plan_many_dft
            (1, &time, howmany, //rank, n, howmany
             noise_data, NULL, howmany, 1, //in,  embed, stride, dist
             noise_data, NULL, howmany, 1, //out, embed, stride, dist
             FFTW_FORWARD, FFTW_MEASURE); //sign, flags
    }
}

dynamic_correlator::~dynamic_correlator() {
    // Delete storage arrays
    if (S_zz)  delete[] S_zz;
    if (S_pm)  delete[] S_pm;
    if (S_BB)  delete[] S_BB;
    if (S_bb)  delete[] S_bb;
    if (noise) delete[] noise;
    if (S_sinsin) delete[] S_sinsin;
    if (S_coscos) delete[] S_coscos;

    // Delete FFTW arrays (complex<double> pointers are mere typecasts)
    if  (corr_data) fftw_free(corr_data);
    if (noise_data) fftw_free(noise_data);

    // Remove FFTW plans
    if  (corr_plan) fftw_destroy_plan(corr_plan);
    if (noise_plan) fftw_destroy_plan(noise_plan);
}

void dynamic_correlator::path4() {
    size_t m = n/4;
    // Gamma -> X
    for (size_t i = 0; i < 4*m; ++i)
        kspace_path.push_back(vec3_int(i,0,0));
    // X -> W
    for (size_t i = 0; i < 2*m; ++i)
        kspace_path.push_back(vec3_int(4*m,i,0));
    // W -> K
    for (size_t i = 0; i < m; ++i)
        kspace_path.push_back(vec3_int(4*m-i,2*m+i,0));
    // K -> Gamma
    for (size_t i = 3*m; i > 0; --i)
        kspace_path.push_back(vec3_int(i,i,0));
    // Gamma -> L
    for (size_t i = 0; i < 2*m; ++i)
        kspace_path.push_back(vec3_int(i,i,i));
    // L -> U
    for (size_t i = 0; i < m; ++i)
        kspace_path.push_back(vec3_int(2*m+2*i,2*m-i,2*m-i));
    // U -> X, note >= in the loop
    for (int i = m; i >= 0; --i)
        kspace_path.push_back(vec3_int(4*m,i,i));
}

void dynamic_correlator::path2() {
    size_t m = n/2;
    // Gamma -> X
    for (size_t i = 0; i < 2*m; ++i)
        kspace_path.push_back(vec3_int(i,0,0));
    // X -> W
    for (size_t i = 0; i < m; ++i)
        kspace_path.push_back(vec3_int(2*m,i,0));
    // W -> L
    for (size_t i = 0; i < m; ++i)
        kspace_path.push_back(vec3_int(2*m-i,m,i));
    // L -> Gamma
    for (size_t i = m; i > 0; --i)
        kspace_path.push_back(vec3_int(i,i,i));
    // Gamma -> K|U -> X' (note <= in loop)
    for (size_t i = 0; i <= 2*m; ++i)
        kspace_path.push_back(vec3_int(i,i,0));
}

void dynamic_correlator::path_noise() {
    // [100] direction
    for (size_t i = 1; i <= 2*n; ++i)
        path_n.push_back(vec3_int(i,0,0));
    // [110] direction
    for (size_t i = 1; i <= 2*n; ++i)
        path_n.push_back(vec3_int(i,i,0));
    // [111] direction
    for (size_t i = 1; i <= n; ++i)
        path_n.push_back(vec3_int(i,i,i));
}

void dynamic_correlator::set(size_t t, int sublat, unsigned type) {
    if (type & SXX) {
        if (what & SXX)
            for (size_t q = 0; q < kspace_path.size(); ++q) 
                S_x[sublat][index_spin(q,t)] = m_fourier(kspace_path[q]);
        return;
    }
    if (type & SYY) {
        if (what & SYY)
            for (size_t q = 0; q < kspace_path.size(); ++q) 
                S_y[sublat][index_spin(q,t)] = m_fourier(kspace_path[q]);
        return;
    }
    if (type & SZZ) {
        if (what & SZZ)
            for (size_t q = 0; q < kspace_path.size(); ++q) 
                S_z[sublat][index_spin(q,t)] = m_fourier(kspace_path[q]);
        if (what & NOISE)
            for (size_t q = 0; q < path_n.size(); ++q)
                noise_vec[index_noise(q,t)] += // note that this is incremental
                    m_fourier(path_n[q]) * direction::axis[sublat][2] *
                    polar(1.0,  -Q * path_n[q] % direction::r[sublat]);
        return;
    }
    if (type & SPM) {
        if (what & SPM)
            for (size_t q = 0; q < kspace_path.size(); ++q) 
                S_p[sublat][index_spin(q,t)] = m_fourier(kspace_path[q]);
        return;
    }
    if (type & ARG) {
        if (what & ARG)
            for (size_t q = 0; q < kspace_path.size(); ++q) 
                S_B[sublat][index_spin(q,t)] = m_fourier(kspace_path[q]);
        return;
    }
    if (type & IMAG) {
        if (what & IMAG)
            for (size_t q = 0; q < kspace_path.size(); ++q) 
                S_b[sublat][index_spin(q,t)] = m_fourier(kspace_path[q]);
        return;
    }
    if (type & SIN) {
        if (what & SIN)
            for (size_t q = 0; q < kspace_path.size(); ++q) 
                S_sin[sublat][index_spin(q,t)] = m_fourier(kspace_path[q]);
        return;
    }
    if (type & COS) {
        if (what & COS)
            for (size_t q = 0; q < kspace_path.size(); ++q) 
                S_cos[sublat][index_spin(q,t)] = m_fourier(kspace_path[q]);
        return;
    }
}

void dynamic_correlator::clear_input() {
    if (noise_data)
        memset(noise_data, 0, 3*path_n.size()*n_time*sizeof(complex<double>));
}

inline double sqr(double x) {return x*x;}

// G in units of 2pi/a_0
double noise_kernel(const vec3_cplx& spin, const vec3_int& q, double G, int n) {
    return sqr(PI / G) * (spin.len2() + norm(spin % normalise(q))) /
        sqr(sin(PI * q.len() / n / G));
}

#define S2 1.414213562373095048801688724209698078569671875376948073176
#define S3 1.732050807568877293527446341505872366942805253810380628055
void dynamic_correlator::add_correlator() {
    // Perform Fourier transforms
    if  (corr_plan) fftw_execute(corr_plan);
    if (noise_plan) fftw_execute(noise_plan);

    // Add mod-squared of appropriate FT's to correlator array
    for (int i = 0; i < 4; ++i)
        for (size_t q = 0; q < kspace_path.size(); ++q)
            for (size_t w = 0; w < n_freq; ++w) {
                if (what & SXX)
                    S_xx[index_corr(q,w)] += norm(S_x[i][index_spin(q,w)]);
                if (what & SYY)
                    S_yy[index_corr(q,w)] += norm(S_y[i][index_spin(q,w)]);
                if (what & SZZ)
                    S_zz[index_corr(q,w)] += norm(S_z[i][index_spin(q,w)]);
                if (what & SPM)
                    S_pm[index_corr(q,w)] += norm(S_p[i][index_spin(q,w)]);
                if (what & ARG)
                    S_BB[index_corr(q,w)] += norm(S_B[i][index_spin(q,w)]);
                if (what & IMAG)
                    S_bb[index_corr(q,w)] += norm(S_b[i][index_spin(q,w)]);
                if (what& SIN)
                    S_sinsin[index_corr(q,w)] += norm(S_sin[i][index_spin(q,w)]);
                if (what& COS)
                    S_coscos[index_corr(q,w)] += norm(S_cos[i][index_spin(q,w)]);
            }

    // Calculate noise correlator by summing all q per direction with kernel
    if (what & NOISE)
        for (size_t w = 0; w < n_freq; ++w) {
            // [100] direction
            for (size_t q = 1; q <= 2*n; ++q)
                noise[index_corr(0,w)] +=
                    noise_kernel(noise_vec[index_noise(q-1, w)],
                                 vec3_int(q,0,0), 4, n) * (q == 2*n ? 1:2);
            // [110] direction
            for (size_t q = 1; q <= 2*n; ++q)
                noise[index_corr(1,w)] +=
                    noise_kernel(noise_vec[index_noise(q+2*n-1, w)],
                                 vec3_int(q,q,0), 4*S2, n) * (q == 2*n ? 1:2);
            // [111] direction
            for (size_t q = 1; q <= n; ++q)
                noise[index_corr(2,w)] +=
                    noise_kernel(noise_vec[index_noise(q+4*n-1,w)],
                                 vec3_int(q,q,q), 2*S3, n) * (q == n ? 1:2);
        }

    ++n_corr;
}
#undef S2
#undef S3

// Filename extensions
const string x_sxx = ".sxx.dyn";
const string x_syy = ".syy.dyn";
const string x_szz = ".szz.dyn";
const string x_sin = ".sin.dyn";
const string x_cos = ".cos.dyn";
const string x_spm = ".spm.dyn";
const string x_arg  = ".arg.dyn";
const string x_imag = ".imag.dyn";
const string x_noise = ".noise.dyn";

void dynamic_correlator::save_corr(const char* root) const {
    // array in which to divide total correlators by n_corr * (#spins)
    double* div = new double[kspace_path.size()*n_freq];
    double rec = 1.0 / n_corr / n_time / (16*n*n*n);

    if (what & SXX) {
        FILE* f = fopen((root+x_sxx).c_str(), "wb");
        for (size_t i = 0; i < kspace_path.size()*n_freq; ++i)
            div[i] = S_xx[i] * rec;
        fwrite(div, sizeof(double), kspace_path.size()*n_freq, f);
        fclose(f);
    }

    if (what & SYY) {
        FILE* f = fopen((root+x_syy).c_str(), "wb");
        for (size_t i = 0; i < kspace_path.size()*n_freq; ++i)
            div[i] = S_yy[i] * rec;
        fwrite(div, sizeof(double), kspace_path.size()*n_freq, f);
        fclose(f);
    }

    if (what & SZZ) {
        FILE* f = fopen((root+x_szz).c_str(), "wb");
        for (size_t i = 0; i < kspace_path.size()*n_freq; ++i)
            div[i] = S_zz[i] * rec;
        fwrite(div, sizeof(double), kspace_path.size()*n_freq, f);
        fclose(f);
    }

    if (what & SPM) {
        FILE* f = fopen((root+x_spm).c_str(), "wb");
        for (size_t i = 0; i < kspace_path.size()*n_freq; ++i)
            div[i] = S_pm[i] * rec;
        fwrite(div, sizeof(double), kspace_path.size()*n_freq, f);
        fclose(f);
    }

    if (what & ARG) {
        FILE* f = fopen((root+x_arg).c_str(), "wb");
        for (size_t i = 0; i < kspace_path.size()*n_freq; ++i)
            div[i] = S_BB[i] * rec;
        fwrite(div, sizeof(double), kspace_path.size()*n_freq, f);
        fclose(f);
    }

    if (what & IMAG) {
        FILE* f = fopen((root+x_imag).c_str(), "wb");
        for (size_t i = 0; i < kspace_path.size()*n_freq; ++i)
            div[i] = S_bb[i] * rec;
        fwrite(div, sizeof(double), kspace_path.size()*n_freq, f);
        fclose(f);
    }

    if (what & NOISE) {
        FILE* f = fopen((root+x_noise).c_str(), "wb");
        for (size_t i = 0; i < 3*n_freq; ++i)
            div[i] = noise[i] * rec;
        fwrite(div, sizeof(double), 3*n_freq, f);
        fclose(f);
    }

    if (what & SIN) {
        FILE* f = fopen((root+x_sin).c_str(), "wb");
        for (size_t i = 0; i < 3*n_freq; ++i)
            div[i] = S_sinsin[i] * rec;
        fwrite(div, sizeof(double), kspace_path.size()*n_freq, f);
        fclose(f);
    }

    if (what & COS) {
        FILE* f = fopen((root+x_cos).c_str(), "wb");
        for (size_t i = 0; i < 3*n_freq; ++i)
            div[i] = S_coscos[i] * rec;
        fwrite(div, sizeof(double), kspace_path.size()*n_freq, f);
        fclose(f);
    }
    
    delete[] div;
}
