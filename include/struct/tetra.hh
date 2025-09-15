/* Class tetra
 * Describes a tetrahedron centre
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

#ifndef tetra_hh
#define tetra_hh

#include <vec3.hh>

class spin;

class tetra {
public:
    // Possible sublattices of tetrahedron centres
    enum class sublattice {A = 1, B = -1};

private:
    // Neighbouring spins
    spin* neigh[4];

    // Sublattice tetrahedron belongs to
    const sublattice m_sl;
    
    // Position (integer in unit-cell based coord.)
    const vec3_int m_pos;
    
public:
    // Auxilliary pointer for predecessor in BFS, used in vison insertion
    spin* prev;
    
    // Can only be constructed with position and sublattice
    tetra(const vec3_int& pos, sublattice sl): m_sl(sl), m_pos(pos) {}

    // Register with a spin
    void reg(spin* s);

    // Get sublattice index
    // int sublat() const {return (int)m_sl;}
    sublattice sublat_enum() const {return m_sl;}

    // Get position
    vec3_int pos() const {return m_pos;}

    // Get corresponding fcc lattice site
    vec3_int fcc() const;

    // Get a neighbour spin
    spin* neighbour(size_t n) const {return neigh[n];}
    
    // tetra* tet_neighbour(size_t n) const {return neigh[n]->t_other(this); }

    // Rotate XY angles of all neighbouring spins
    void rotate(double alpha);
};

#endif
