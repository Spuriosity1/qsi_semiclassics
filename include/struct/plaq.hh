/* Class plaq
 * Describes a plaquette (dual lattice site)
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

#ifndef plaq_hh
#define plaq_hh

#include <vec3.hh>
#include <complex>

class spin;
class ptetra;

class plaq {
public:
    // Possible sublattices of plaquettes
    enum class sublattice {zero = 0, one = 1, two = 2, three = 3};
    
protected:
    // Position (integer in unit-cell based coord.)
    const vec3_int m_pos;

    // Dual Tetrahedron centres connected to this plaquette
    ptetra* m_A;
    ptetra* m_B;

    // Constituent spin
    spin* m_spin[6];

    // Sublattice index
    const sublattice m_sl;

public:
    //---------- CONSTRUCTION, FITTING INTO STRUCTURE -------------------------
    
    // Can only be constructed with position and sublattice
    plaq(const vec3_int& pos, const sublattice sl): m_pos(pos), m_sl(sl) {}
    
    // Virtual destructor to make sure descendants are destroyed
    virtual ~plaq() {}

    // Register with tetrahedron centres
    void reg(ptetra* t);
    // Registered tetrahedron centres
    ptetra* tA() const {return m_A;}
    ptetra* tB() const {return m_B;}
    ptetra* t_other(const ptetra* p) const;

    // Register with spin
    void reg(spin* s, int index);
    // Registered spins
    spin& operator[](int index) const;

    // Get sublattice index
    int sublat() const {return (int)m_sl;}

    /**
     * @brief Returns signed orientation of this spin relative to the specified hexagon
     * Specifically: Constructing a vector from the centre of the A dual lattice to the B dual lattice,
     * traverse the hexagon 
     * 
     * @param s spin
     * @return +1 if positive, -1 if negative, 0 if spin is not a member of this plaquette
     */
    int spin_sign(const spin& s);
    
    // Get position
    vec3_int pos() const {return m_pos;}

    // Get nearest fcc lattice site
    vec3_int fcc() const;

    //------------------ COUPLING PARAMETERS-------------------
    std::complex<double> g_constant = 1.0;
    // Roskhar-Kivelson potential
    double RK_potential = 0;


    //---------- PHYSICAL MEASUREMENTS ----------------------------------------

    // Gauge field exp(iB), as seen from the A tetrahedron
    std::complex<double> ring() const;

    /* Gauge field B, as seen from the A tetrahedron
     * restricted between -pi..pi, as per c++ standard */
    double B() const {return std::arg(ring());}

    // Energy of the plaquette
    double ring_energy() const;

    /* Total energy of all plaquettes that have at least one spin in common
     * with this one.
     * Used to assess energetic cost of MC steps */
    double energy_MC() const;

    // double plaq::ising2nn_energy() const;
    
    double ising3nn_energy() const;
    double ising3nn_energy(const vec3* array) const;

    // Local heat current
    vec3 heat_current() const;

    //---------- HANDLING OTHER ARRAYS ----------------------------------------
    std::complex<double> ring(const vec3* array) const;

    double ring_energy(const vec3* array) const;

    double energy_MC(const vec3* array) const;

    vec3 heat_current(const vec3* array) const;

    //---------- VISON HOPPING ------------------------------------------------
    /* Introduces a positive vison on ptetra in parameter and a negative vison
     * on the other one. 
     * Throws an error message if given ptetra not next to this plaq. */
    void vison_hop(const ptetra* pos, double charge=1);
};

#endif
