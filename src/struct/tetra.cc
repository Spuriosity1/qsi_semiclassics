/* Implementation of functions in class tetra 
 *
 * Created on 11/06/2018
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

#include "tetra.hh"
#include "spin.hh"
#include <direction.hh>
#include <vec3.hh>
#include <complex>

void tetra::reg(spin* s) {neigh[s->sublat()] = s;}

vec3_int tetra::fcc() const {
    using namespace direction;
    switch(m_sl) {
    case sublattice::A: return m_pos - diamond[0];
    case sublattice::B: return m_pos - diamond[1];
    default: throw std::runtime_error("Invalid tetrahedron type");
    }
}

void tetra::rotate(double alpha) {
    using namespace std;
    complex<double> z = polar(1.0,alpha);

    for(int i = 0; i < 4; ++i) // TODO check eventual definition in spin (what does this mean????)
        neigh[i]->xy() *= z;
}
 