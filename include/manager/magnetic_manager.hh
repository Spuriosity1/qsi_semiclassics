/* Class magnetic_manager
 * Retrieves static B-field (any variant) correlators written to file by some
 * incarnation of static_correlator, and uses them to construct correlator,
 * vison density plots/cuts, etc.
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

#ifndef magnetic_manager_hh
#define magnetic_manager_hh

#include <cstdio>
#include <complex>

class magnetic_manager {
private:
    // Size of the original system, corresponding distance of q-points
    const int n;
    const double Q;

    /* The arrays of correlators <B_i(q) B_j(-q)>
     * Each pointer points to a dynamical complex array 2n^2 in size
     * Indexing: corr[i][j][alpha][beta][q] */
    std::complex<double>* S_BB[4][4];
    
    /* Generates the array index of an arbitrary q-point in the (hhk) plane
     * Unit of q_x,q_z is 2pi/n -> must be integer to be a valid q-point
     * in the periodic box
     * Takes fcc reciprocal lattice into account */
    size_t index(int x, int z) const;

public:
    magnetic_manager(int n);
    ~magnetic_manager();

    /* Load magnetic field correlators from the block of file specified by the
     * parameter <block>.
     * A "block" means 32n^2 complex numbers, the size of a full correlator
     * specification. This accommodates files containing everything. */
    void load(std::FILE* f, unsigned block = 0u);

    /* Add to field correlators <w> times the contents of block <block>
     * of file <f> */
    void add(std::FILE* f, double w, unsigned block = 0u);

    // Clear contents (fill with zeroes)
    void clear();

    // Evaluate overall BB correlator
    double BB(int qx, int qz) const;

    // Evaluate overall vison correlator (accurate only for "arg" B-field)
    double vison(int qx, int qz) const;
};

#endif
 
