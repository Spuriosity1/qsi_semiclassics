/* Directions defining a tetrahedron and a plaquette in the pyrochlore lattice;
 * local axis directions; pyrochlore direction vectors.
 * Integer vectors in units of a_fcc/8.
 *
 * Created on 09/06/2018
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

#ifndef direction_hh
#define direction_hh

#include "vec3.hh"

namespace direction {
    // Positions of pyrochlore sites relative to "A" diamond sites
    extern const vec3_int pyro[4];

    // Positions of diamond sites relative to FCC sites
    extern const vec3_int diamond[2];

    // Positions of FCC lattice sites in the "rare earth" lattice
    extern const vec3_int fcc_Dy[4];

    // Positions of FCC lattice sites in the "titanium" lattice
    extern const vec3_int fcc_Ti[4];

    /* Directions of pyrochlore sites in a plaquette from the corresponding
     * dual pyrochlore site (plaquette centre)
     * First index: plaquette sublattice, second index: spin index
     * Works for neighbouring plaquettes of spins, too. */
    extern const vec3_int plaqt[4][6];

    /* Cubic lattice vectors. Trivial, but good practice to define these here.
     * First index: x,y,z.
     * Second index: also x,y, z. It's diagonal.
     * */
    extern const vec3_int cubic[3];

    // pyro scaled down by a factor of 8
    extern const vec3 r[4];

    /* Local axis vectors for each pyrochlore sublattice
     * First index: sublattice, second index: local axis direction */
    extern const vec3 axis[4][3];
}

#endif
