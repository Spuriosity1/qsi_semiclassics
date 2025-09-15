/* Class ptetra
 * Describes a tetrahedron centre of the dual (plaquette) lattice
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

#ifndef ptetra_hh
#define ptetra_hh

#include <vec3.hh>

class plaq;

class ptetra {
public:
    // Possible sublattices of tetrahedron centres
    enum class sublattice {A = 1, B = -1};

private:
    // Neighbouring plaquettes
    plaq* neigh[4];

    // Sublattice tetrahedron belongs to
    const sublattice m_sl;
    
    // Position (integer in unit-cell based coord.)
    const vec3_int m_pos;

    // Memoised vison charge
    mutable int m_vison;

public:
    // Auxilliary pointer for predecessor in BFS, used in vison insertion
    plaq* prev;

    // Index in qsi simulation
    const size_t index;
    
    // Can only be constructed with position, sublattice, and index in qsi
    ptetra(const vec3_int& pos, sublattice sl, size_t index):
        m_sl(sl),
        m_pos(pos),
        index(index)
    {}

    // Register with a plaquette
    void reg(plaq* s);

    // Get sublattice index
    // int sublat() const {return (int)m_sl;}
    sublattice sublat_enum() const {return m_sl;}

    // Get position
    vec3_int pos() const {return m_pos;}

    // Get corresponding fcc lattice site
    vec3_int fcc() const;

    // Get a neighbour
    plaq* neighbour(size_t n) const {return neigh[n];}

    /* Get vison charge
     * If memo is true, returns memoised vison charge (use with care!)
     * If memo if false, evaluates vison charges and updates m_vison */
    int vison(bool memo = false) const;

    // Gets vison charge **without rounding**
    double vison_raw() const;

    // Energy of plaquettes making up tetrahedron
    double energy() const;
    double energy(const vec3* array) const;
};

#endif
