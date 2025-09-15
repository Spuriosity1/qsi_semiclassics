/* Class spin
 * Describes a Heisenberg spin on a pyrochlore site.
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

#ifndef spin_hh
#define spin_hh

#include <vec3.hh>
#include <complex>
#include <random>

class tetra;
class plaq;

// static double plaq::Ising_3nn;

class spin {
public:
    // Possible sublattices of spins
    enum class sublattice {zero = 0, one = 1, two = 2, three = 3};

protected:
    // Heisenberg spin vectors (stored in contiguous arrays, hence pointers)
    // Actual spin
    vec3* const m_spin;
    // Index in any storage array (for complicated MC and dynamics)
    const size_t m_index;

    // Position (integer in unit-cell based coord.)
    const vec3_int m_pos;

    // // Second- and third-neighbour spins
    // spin* m_second_neighbour[12];
    // spin* m_third_neighbour[6];

    // Tetrahedron centres connected to this spin
    tetra* m_A;
    tetra* m_B;

    // Plaquettes the spin belongs to
    plaq* m_plaq[6];

    // Sublattice index
    const sublattice m_sl;

    // g-factor, set this to zero to exclude from DSF and SSF
    double g_factor=1.;
    
public:    

    //---------- CONSTRUCTION, FITTING INTO STRUCTURE -------------------------

    /* Can only be constructed with position, sublattice, array index,
     * and the vector array from which m_spin is taken */
    spin(const vec3_int& pos, const sublattice sl, const size_t index,
         vec3* array):
        m_spin(array+index),
        m_index(index),
        m_pos(pos),
        m_sl(sl)
    {
        // Initialise ordered along x-axis, a classical ground state for g > 0
        *m_spin = vec3(1.0, 0.0, 0.0);
    }

    
    
    // Virtual destructor to make sure descendants are destroyed
    virtual ~spin() {}

    // Register spin with tetrahedron centres
    void register_plaquette(tetra* t);
    // Registered tetrahedron centres
    tetra* tA() const {return m_A;}
    tetra* tB() const {return m_B;}
    tetra* t_other(const tetra* p) const;

    // Register with plaquette
    void register_plaquette(plaq* p, int index);
    // Registered plaquettes
    plaq& operator[](int index) const;

    // Get sublattice index
    int sublat() const {return (int)m_sl;}
    
    // Get position
    vec3_int pos() const {return m_pos;}

    // Get nearest fcc lattice site
    vec3_int fcc() const;

    // Find and register second neighbours
    // void register_spin_neighbours();

    // Get second-neighbour spin sites
    // const std::vector<const spin*> second_neighbours() const;


    //---------- ACCESSING & SETTING THE SPIN COMPONENTS ----------------------

    // The spin vector, constant & editable
    const vec3& heis() const {return *m_spin;}
    vec3& heis() {return *m_spin;}

    // The XY components of the spin as a complex number
    const std::complex<double>& xy() const {return m_spin->xy();}
    std::complex<double>& xy() {return m_spin->xy();}

    // The Ising component of the spin
    double ising() const {return (*m_spin)[2];}
    double& ising() {return (*m_spin)[2];}

    double get_g_factor() const {return this->g_factor;}
    void set_g_factor(double _g_factor) {this->g_factor = _g_factor;}

    //---------- PHYSICAL MEASUREMENTS ----------------------------------------

    // Effective field acting on the spin (always in the XY plane)
    std::complex<double> field_cplx() const;
    vec3 field() const;

    // Time derivative of spin direction
    vec3 diff() const {return (*m_spin) * field();}

//    double zfield_ising_3nn() const;
//    double zfield_ising_3nn(const vec3* array) const;

    // Effective energy of spin (total energy of plaquettes containing it)
    double energy() const;

    //---------- MONTE CARLO --------------------------------------------------

    // Find new XY angle using a von Mises distribution
    void MC_angle(double T, std::mt19937& g);

    /*---------- HANDLING OTHER ARRAYS ----------------------------------------
     * These functions do what their no-argument counterparts would if m_spin
     * came from the array argument. */
    
    const vec3& heis(const vec3* array) const {return array[m_index]; }
    vec3& heis(vec3* array) const {return array[m_index]; }
    
    const std::complex<double>& xy(const vec3* array) const {
        return array[m_index].xy();
    }
    std::complex<double>& xy(vec3* array) const {
        return array[m_index].xy();
    }

    double  ising(const vec3* array) const {return array[m_index][2];}
    double& ising(vec3* array) const {return array[m_index][2];}
    
    std::complex<double> field_cplx(const vec3* array) const;

    vec3 field(const vec3* array) const {
        vec3 f(field_cplx(array));
//        f[2] += -zfield_ising_3nn(array);
        return f; }

    vec3 diff(const vec3* array) const {return array[m_index] * field(array);}

    double energy(const vec3* array) const;

    //---------- ADDING SPINONS -----------------------------------------------
    // Add +q spinon charge to tetra "pos" and -q to the one on the other side
    void spinon_hop (double q, const tetra* pos);
};

#endif
