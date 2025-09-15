/* Implementation of class spin_manager
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

#include "spin_manager.hh"
#include <vec3.hh>
#include <direction.hh>
#include <misc.hh>
#include <cstdio>
#include <cstring>
#include <complex>
#include <cmath>

using namespace std;

#define TWOPI 6.283185307179586476925286766559005768394338798750211641949

spin_manager::spin_manager(int n):
    n(n),
    Q(TWOPI/n)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            S_zz[i][j] = new complex<double>[2ul*n*n];
            S_pm[i][j] = new complex<double>[2ul*n*n];
        }
}

spin_manager::~spin_manager() {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            delete[] S_zz[i][j];
            delete[] S_pm[i][j];
        }
}

void spin_manager::clear() {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            std::fill_n(S_zz[i][j], 2ul*n*n, std::complex<double>(0) );
            std::fill_n(S_pm[i][j], 2ul*n*n, std::complex<double>(0) );
            // memset(S_zz[i][j], 0, 2ul*n*n*sizeof(complex<double>));
            // memset(S_pm[i][j], 0, 2ul*n*n*sizeof(complex<double>));
        }
}

void spin_manager::load(FILE* f) {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) 
            fread(S_zz[i][j], sizeof(complex<double>), 2ul*n*n, f);
    
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) 
            fread(S_pm[i][j], sizeof(complex<double>), 2ul*n*n, f);
}

void spin_manager::add(FILE* f, double w) {
    complex<double>* temp = new complex<double>[2ul*n*n];
    
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            fread(temp, sizeof(complex<double>), 2ul*n*n, f);
            for (size_t q = 0; q < 2ul*n*n; ++q)
                S_zz[i][j][q] += temp[q]*w;
        }
    
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            fread(temp, sizeof(complex<double>), 2ul*n*n, f);
            for (size_t q = 0; q < 2ul*n*n; ++q)
                S_pm[i][j][q] += temp[q]*w;
        }
    
    delete[] temp;
}

size_t spin_manager::index(int x, int z) const {
    z = mod(z,2*n);
    if (z >= n) {
        x -= n;
        z -= n;
    }
    x = mod(x,2*n);
    return x*n + z;
}

void spin_manager::tensor(int x, int z, unsigned which) {
    qx = x;
    qz = z;
    qq = vec3(x,x,z)*Q;
    size_t q = index(x,z);
    for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
            m_tensor[m][n] = 0.0;
    // Miserere mei Deus secundum magnam misericordiam tuam
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            using namespace direction;
            complex<double> phase = polar(1.0, qq%(r[j]-r[i]));
            for (int m = 0; m < 3; ++m)
                for (int n = 0; n < 3; ++n) {
                    if (which & ISING)
                        m_tensor[m][n] += phase * S_zz[i][j][q] *
                            axis[i][2][m] * axis[j][2][n];
                    if (which & XY)
                        m_tensor[m][n] += 0.5 * phase * S_pm[i][j][q] *
                            (axis[i][0][m] * axis[j][0][n] +
                             axis[i][1][m] * axis[j][1][n]);
                }
        }
}

// Polarisation of neutrons
#define S2 1.414213562373095048801688724209698078569671875376948073176
const vec3 p(1/S2, -1/S2, 0);
#undef S2

double spin_manager::neutron_nsf() const {
    complex<double> retval = 0.0;
    for (int m = 0; m < 3; ++m)
        for (int n = 0; n < 3; ++n)
            retval += m_tensor[m][n] * p[m] * p[n];
    return retval.real();
}

double spin_manager::neutron_sf() const {
    complex<double> retval = 0.0;
    if (qx || qz) {
        // non-zero wave vector -> same as NSF with orthogonal vector
        vec3 o = p*qq;
        o.normalise();
        for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
                retval += m_tensor[m][n] * o[m] * o[n];
    } else {
        // zero wave vector -> delta^{mu,nu} - p^mu p^nu
        for (int m = 0; m < 3; ++m)
            for (int n = 0; n < 3; ++n)
                retval += m_tensor[m][n] * ((m==n?1:0) - p[m]*p[n]);
    }
    return retval.real();
}

double spin_manager::Szz(int x, int z) const {
    complex<double> retval = 0.0;
    size_t q = index(x,z);
    vec3 vq = vec3(x,x,z)*Q;
    
    using namespace direction;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            retval += polar(1.0, vq%(r[j]-r[i])) * S_zz[i][j][q];

    return retval.real();
}
