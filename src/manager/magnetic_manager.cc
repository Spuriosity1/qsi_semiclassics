/* Implementation of class mangetic_manager
 *
 * Created on 08/10/2018
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

#include "magnetic_manager.hh"
#include <vec3.hh>
#include <direction.hh>
#include <misc.hh>
#include <cstdio>
#include <cstring>
#include <complex>
#include <cmath>

using namespace std;

#define TWOPI 6.283185307179586476925286766559005768394338798750211641949

magnetic_manager::magnetic_manager(int n):
    n(n),
    Q(TWOPI/n)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) 
            S_BB[i][j] = new complex<double>[2ul*n*n];
}

magnetic_manager::~magnetic_manager() {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) 
            delete[] S_BB[i][j];
}

void magnetic_manager::clear() {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) 
            std::fill_n(S_BB[i][j], 2ul*n*n, std::complex<double>(0));
            // memset(S_BB[i][j], 0, 2ul*n*n*sizeof(complex<double>));
}

void magnetic_manager::load(FILE* f, unsigned block) {
    fseek(f, 32ul*n*n*sizeof(complex<double>)*block, SEEK_SET);
    
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) 
            fread(S_BB[i][j], sizeof(complex<double>), 2ul*n*n, f);
}

void magnetic_manager::add(FILE* f, double w, unsigned block) {
    fseek(f, 32ul*n*n*sizeof(complex<double>)*block, SEEK_SET);
    complex<double>* temp = new complex<double>[2ul*n*n];
    
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            fread(temp, sizeof(complex<double>), 2ul*n*n, f);
            for (size_t q = 0; q < 2ul*n*n; ++q)
                S_BB[i][j][q] += temp[q]*w;
        }
    
    delete[] temp;
}

size_t magnetic_manager::index(int x, int z) const {
    z = mod(z,2*n);
    if (z >= n) {
        x -= n;
        z -= n;
    }
    x = mod(x,2*n);
    return x*n + z;
}

/**
 * @brief Exaluates overall correlator S_BB. 
 * 
 * Convertes the raw SL-SL correlation blocks used as input into an overall correlator, 
 * re-adding the papropriate phases.
 * 
 * @param x Position h on the hhl plane
 * @param z Posision l on the hhl plane
 * @return double overall correlation
 */
double magnetic_manager::BB(int x, int z) const {    
    complex<double> retval = 0.0;
    size_t q = index(x,z);
    vec3 vq = vec3(x,x,z)*Q;
    
    using namespace direction;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            retval += polar(1.0, vq%(r[j]-r[i])) * S_BB[i][j][q];

    return retval.real();
}


double magnetic_manager::vison(int x, int z) const {
    complex<double> retval = 0.0;
    size_t q = index(x,z);
    vec3 vq = vec3(x,x,z)*Q;
    
    using namespace direction;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            retval += 4.0 * polar(1.0, vq%(r[j]-r[i])) *
                sin(vq%r[i]) * sin(vq%r[j]) * S_BB[i][j][q];

    return retval.real();
}
