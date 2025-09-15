/* Implementation of functions in class plaq 
 *
 * Created on 11/06/2018
 * Copyright (C) 2018, 2019 Attila Szabó <as2372@cam.ac.uk>
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

#include "spin.hh"
#include "ptetra.hh"
#include "plaq.hh"
#include <direction.hh>
#include <vec3.hh>
#include <misc.hh>
#include <complex>
#include <exception>

using namespace std;

void plaq::reg(ptetra* t) {
    switch(t->sublat_enum()) {
    case ptetra::sublattice::A: m_A = t;
        break;
    case ptetra::sublattice::B: m_B = t;
        break;
    default: throw std::runtime_error("Invalid tetrahedron type");
    }
}

ptetra* plaq::t_other(const ptetra* p) const { 
    if (m_A == p)
        return m_B;
    else if (m_B == p)
        return m_A;
    else
        throw std::runtime_error("Invalid tetrahedron neighbour");
}

void plaq::reg(spin* p, int index) {m_spin[mod(index,6)] = p;}

spin& plaq::operator[](int index) const {return *(m_spin[mod(index,6)] ); }

vec3_int plaq::fcc() const {
    using namespace direction;
    return m_pos - pyro[sublat()];
}



complex<double> plaq::ring() const {
    return
        m_spin[0]->xy() * conj(m_spin[1]->xy()) * 
        m_spin[2]->xy() * conj(m_spin[3]->xy()) * 
        m_spin[4]->xy() * conj(m_spin[5]->xy());
}

complex<double> plaq::ring(const vec3* a) const {
    return
        m_spin[0]->xy(a) * conj(m_spin[1]->xy(a)) * 
        m_spin[2]->xy(a) * conj(m_spin[3]->xy(a)) * 
        m_spin[4]->xy(a) * conj(m_spin[5]->xy(a));
}

/**
 * @brief Computes the energy $\sum_{<<ij>>} S^z_iS^z_j$ where one of $i$ or $j$ are on the ring.
 * Each spin has 12 2nn's (6 3nn's)
 * 
 * @return double 
 */
// double plaq::ising2nn_energy() const {
//    throw std::logic_error("Not Implemented");

// }

/**
 * @brief Computes the energy $\sum_{<<ij>>} S^z_iS^z_j$ where one of $i$ or $j$ are on the ring.
 * Each spin has 6 3nn's, all of which lie on the same sublattice.
 * 
 * @return double 
 */
//double plaq::ising3nn_energy() const {
//    double etot = 0;
//    for (int i=0; i<6; i++){
//        etot += m_spin[i]->zfield_ising_3nn()*m_spin[i]->ising();
//    }
//    return etot;
//}
//
//// array version
//double plaq::ising3nn_energy(const vec3* a) const {
//    double etot = 0;
//    for (int i=0; i<6; i++){
//        etot += m_spin[i]->zfield_ising_3nn(a)*m_spin[i]->ising(a);
//    }
//    return etot;
//}

/* Plaquettes: 
 * In the first approximation, sum "energies" of 6 spins making up the
 * plaquette. This results in overcounting:
 *     *this counted 6 times
 *     neighbours counted twice
 *     other affected spins counted once
 * Therefore, the "energies" of the two dual tetrahedra are subtracted, and
 * energy() 3 times more. 
 * 
 * Ising parts: m_spin[i] counts the Ising field, this is all we need to remember (no double counting)
 * */
double plaq::energy_MC() const {
    double retval = 0.0;
    for (int i = 0; i < 6; ++i)
        retval += m_spin[i]->energy();
    retval -= m_A->energy() + m_B->energy();
    retval -= 3*ring_energy();
    return retval;
}

double plaq::energy_MC(const vec3* a) const {
    double retval = 0.0;
    for (int i = 0; i < 6; ++i)
        retval += m_spin[i]->energy(a);
    retval -= m_A->energy(a) + m_B->energy(a);
    retval -= 3*ring_energy(a);
    return retval;
}

vec3 plaq::heat_current() const {
    using namespace direction;
    vec3 retval;
    for (int i = 0; i < 6; ++i)
        retval += 0.125 * plaqt[sublat()][i] *
            real( (*this)[i].diff().xy() * conj( (*this)[i+1].xy() ) *
                  (*this)[i+2].xy()      * conj( (*this)[i+3].xy() ) *
                  (*this)[i+4].xy()      * conj( (*this)[i+5].xy() ) );
    return retval;
}

vec3 plaq::heat_current(const vec3* a) const {
    using namespace direction;
    vec3 retval;
    for (int i = 0; i < 6; ++i)
        retval += 0.125 * plaqt[sublat()][i] *
            real( (*this)[i].diff(a).xy() * conj( (*this)[i+1].xy(a) ) *
                  (*this)[i+2].xy(a)      * conj( (*this)[i+3].xy(a) ) *
                  (*this)[i+4].xy(a)      * conj( (*this)[i+5].xy(a) ) );
    return retval;
}

#define PI4 0.785398163397448309615660845819875721049292349843776455243

void plaq::vison_hop(const ptetra* pos, double charge) {
    /* The hopping operation consists of rotating the XY angles of the 
     * surrounding spins by 45 degrees with alternating signs. 
     * The first sign should be negative if the positive vison is to be created
     * on the A tetrahedron. */
    double q=charge*PI4;

    complex<double> turn;
    if (pos == m_A)
	turn = polar(1.0, -q);
    else if (pos == m_B)
	turn = polar(1.0, q);
    else
	throw std::runtime_error("Invalid ptetra supplied");

    for (int i = 0; i < 6; ++i) {
	m_spin[i]->xy() *= turn;
	turn = conj(turn);
    }
}

double plaq::ring_energy(const vec3* array) const {
    std::complex O = ring(array);
    return -std::real(O*g_constant + 2*RK_potential*O*std::conj(O));
}

double plaq::ring_energy() const {
    std::complex O = ring();
    return -std::real(O*g_constant + 2*RK_potential*O*std::conj(O));
}

int plaq::spin_sign(const spin& s){
    for (int i=0; i<6; i++){
        if (this->m_spin[i] == &s){
            return i%2==0 ? 1 : -1;
        }
    }
    return 0;
}
