/* Implementation of functions in class spin 
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

#include "spin.hh"
#include "tetra.hh"
#include "plaq.hh"
#include <direction.hh>
#include <vec3.hh>
#include <misc.hh>
#include <random.hh>
#include <complex>
#include <random>




using namespace std;

void spin::register_plaquette(tetra* t) {
    switch(t->sublat_enum()) {
    case tetra::sublattice::A: m_A = t;
        break;
    case tetra::sublattice::B: m_B = t;
        break;
    default: throw std::runtime_error("Invalid tetrahedron type");
    }
}

tetra* spin::t_other(const tetra* p) const {
    if (m_A == p)
        return m_B;
    else if (m_B == p)
        return m_A;
    else
        throw std::runtime_error("Invalid tetrahedron neighbour");
}

void spin::register_plaquette(plaq* p, int index) {m_plaq[mod(index,6)] = p;}

plaq& spin::operator[](int index) const {return *(m_plaq[mod(index,6)] ); }

vec3_int spin::fcc() const {
    using namespace direction;
    return m_pos - pyro[sublat()];
}

// // Find and register second neighbours
// void spin::register_spin_neighbours(){
//     unsigned I2= 0;
//     unsigned I3 = 0;

//     for (char i=0; i<4; i++){
//         spin* s = this->tA()->neighbour(i);
//         if (s != this) {
//             for (char j = i+1; j!= i; j = (j+1)%4){
//                 spin* s2 = s->tB()->neighbour(j);
//                 if (s2->m_sl == this->m_sl){
//                     this->m_third_neighbour[I3] = s;
//                     I3++;
//                 } else {
//                     this->m_second_neighbour[I2] = s;
//                     I2++;
//                 }
//             }
//         }
//     }

//     for (char i=0; i<4; i++){
//         spin* s = this->tB()->neighbour(i);
//         if (s != this) {
//             for (char j = i+1; j!= i; j = (j+1)%4){
//                 spin* s2 = s->tA()->neighbour(j);
//                 if (s2->m_sl == this->m_sl){
//                     this->m_third_neighbour[I3] = s;
//                     I3++;
//                 } else {
//                     this->m_second_neighbour[I2] = s;
//                     I2++;
//                 }
//             }
//         }
//     }

// }

// // Get second-neighbour spin sites
// const std::vector<const spin*> spin::second_neighbours() const{

// }
//
//// 6 terms - not literally all geometric nearest neighbours (neglects cross-loop coupling)
//double spin::zfield_ising_3nn() const{
//    double retval = 0;
//
//    for (unsigned j=0; j<4; j++){
//        unsigned tmp = (j + (unsigned) this->m_sl)%4;
//        retval += this->tB()->neighbour(tmp)->tA()->neighbour((size_t) this->m_sl)->ising();
//        retval += this->tA()->neighbour(tmp)->tB()->neighbour((size_t) this->m_sl)->ising();
//    }
//    return Ising_3nn*retval;
//}
//
//double spin::zfield_ising_3nn(const vec3* array) const{
//    double retval = 0;
//    
//    for (unsigned j=0; j<4; j++){
//        unsigned tmp = (j + (unsigned) this->m_sl)%4;
//        retval += this->tB()->neighbour(tmp)->tA()->neighbour((size_t) this->m_sl)->ising(array);
//        retval += this->tA()->neighbour(tmp)->tB()->neighbour((size_t) this->m_sl)->ising(array);
//        ++j;
//    }
//    return Ising_3nn*retval;
//}


// Returns the local field arising from the six nearest neighbour hexagons
// Specifically, this is F = -dU/dS^+_j
// I.e. such that E = -spin. field() = -Re[w* . field_cplx()]
complex<double> spin::field_cplx() const {
    complex<double> retval = 0.0;

    // H = \sum_p g_constant * O + g_constant^* * O^dagger + \sum_{second_neighbour} S^zS^z
    // This function is only used in the angular XY step, which is the same regardles of any Ising interactions.

    // Iterate over all six neighbouring plaquettes
    for (int i = 0; i < 6; ++i) {
        // Determine the orientation
        // We are traversing the plaquette CCW around the dualA -> dualB vector
        // "this" site is "positive" if we go diamondA->diamondB, negative otherwise
        // (*this)[i][i] is this spin 
        // The convention defined when these pointers are registered (in qsi.cc)
        // is such that [0] *always* corresponds to an diamondA->diamondB pyrochlore site
        // using the plaquette's orientation
        // THEREFORE *this[i] is "positive" when i is even, "negative" when i is odd.
        // This distinction is only relevant for complex g: the O + O\dagger cancels it otherwise.

        

        complex<double> dOdS = 0; // represents d[gO + g* O\dagger]/d[S+_j]
        plaq& adj_plaq = (*this)[i];
        
        dOdS = adj_plaq[i+1].xy() * conj( adj_plaq[i+2].xy() ) *
                adj_plaq[i+3].xy() * conj( adj_plaq[i+4].xy() ) *
                adj_plaq[i+5].xy() *
                ( i%2==1 ? adj_plaq.g_constant : conj(adj_plaq.g_constant) );
        

        #ifdef DEBUG
        if (i%2 == 0){
                if (adj_plaq.spin_sign(*this) != 1){
                    throw std::runtime_error("spin sign is broken");
                }
        } else {
                if (adj_plaq.spin_sign(*this) != -1){
                    throw std::runtime_error("spin sign is broken");
                }
        }
        #endif

        retval += dOdS;
        

        // O.Od + Od.O = [O + Od]^2, since the double S+ and double S- vanish
        // Factor of 4: one factor of 2 from derivative, one from O + Od = 2 Re (O)
        retval += -2*adj_plaq.RK_potential*dOdS*adj_plaq.ring_energy();   
    }

    return retval;
}
                   
double spin::energy() const {
    double retval = 0.0;
    for (int i = 0; i < 6; ++i)
        retval += m_plaq[i]->ring_energy();

    // retval += this->ising()*this->zfield_ising_3nn();
    return retval;
}

// Returns the local field arising from the six nearest neighbour hexagons
complex<double> spin::field_cplx(const vec3* a) const {
    complex<double> retval = 0.0;
    for (int i = 0; i < 6; ++i) {
        plaq& adj_plaq = (*this)[i];
        complex<double> dOdS = adj_plaq[i+1].xy(a) * conj( adj_plaq[i+2].xy(a) ) *
                    adj_plaq[i+3].xy(a) * conj( adj_plaq[i+4].xy(a) ) *
                    adj_plaq[i+5].xy(a) *
                    ( i%2==1 ? adj_plaq.g_constant : conj(adj_plaq.g_constant) );

        retval += dOdS;
        retval += -2*adj_plaq.RK_potential*dOdS*adj_plaq.ring_energy(a);   
    }
    return retval;
}
                   
double spin::energy(const vec3* a) const {
    double retval = 0.0;
    for (int i = 0; i < 6; ++i){
        retval += m_plaq[i]->ring_energy(a);
    }
    // retval += this->ising(a)*this->zfield_ising_3nn(a);
    return retval;
}

// Rotates the spin about the Sz axis following a von Mises distribution guided by the size of h
void spin::MC_angle(double T, mt19937& g) {
    complex<double> h = field_cplx();
    complex<double> S = xy();
    if (T == 0.0) {
        // At zero temperature, need to align S^+ with h^+
        xy() = polar(abs(S), arg(h));
    } else {
        // The distribution sampled is ~exp(beta|h||S| cos(arg(S)-arg(h)))
        xy() = abs(S) * von_mises(g, arg(h), abs(h)*abs(S)/T);
    }
}

inline double sqr(double x) {return x*x;}

void spin::spinon_hop (double q, const tetra* pos) {
    if (pos == m_A)
        ising() += q;
    else if (pos == m_B)
        ising() -= q;
    else
        throw std::runtime_error("Invalid tetra supplied");
    // Normalise xy component
    xy() = std::polar(sqrt(1.0 - sqr(ising())),  std::arg(xy())  );
}

 vec3 spin::field() const{
    vec3 F(field_cplx());
    // F[2] += -zfield_ising_3nn();
    return F;
}
