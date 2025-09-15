/* Implementation of functions in class ptetra 
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

#include "ptetra.hh"
#include "plaq.hh"
#include <direction.hh>
#include <vec3.hh>
#include <cmath>
#include <iostream>

void ptetra::reg(plaq* s) {neigh[s->sublat()] = s;}

vec3_int ptetra::fcc() const {
    using namespace direction;
    switch(m_sl) {
    case sublattice::A: return m_pos - diamond[0];
    case sublattice::B: return m_pos - diamond[1];
    default: throw std::runtime_error("Invalid tetrahedron type");
    }
}


inline double ptetra::vison_raw() const {
    double flux = 0.0;
    for (int i = 0; i < 4; ++i)
        flux += neigh[i]->B();
    flux *= (int) sublat_enum();
    return flux / (2*M_PI);
}

int ptetra::vison(bool memo) const {
    if (!memo) {
        m_vison = lround(vison_raw());
    }
    return m_vison;
}

double ptetra::energy() const {
    double retval = 0.0;
    for (int i = 0; i < 4; ++i)
        retval += neigh[i]->ring_energy();
    return retval;
}

double ptetra::energy(const vec3* a) const {
    double retval = 0.0;
    for (int i = 0; i < 4; ++i)
        retval += neigh[i]->ring_energy(a);
    return retval;
}

