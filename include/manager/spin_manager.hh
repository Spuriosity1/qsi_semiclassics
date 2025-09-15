/* Class spin_manager
 * Retrieves static spin correlators written to file by class static_correlator
 * and uses them to construct neutron scattering plots etc.
 *
 * Created on 01/08/2018
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

#ifndef spin_manager_hh
#define spin_manager_hh

#include <cstdio>
#include <complex>
#include <vec3.hh>

class spin_manager {
public:
    /* Flags for including Ising and transverse components
     * of the spin correlators in the tensor output */
    static const unsigned ISING = 1u, XY = 2u;
    
private:
    // Size of the original system, corresponding distance of q-points
    const int n;
    const double Q;
    
    // The tensor <S^mu(q) S^nu(-q)> in global indices
    std::complex<double> m_tensor[3][3];
    // Coordinates to which it corresponds
    int qx, qz;
    vec3 qq;

    /* The arrays of correlators <S^z_i(q) S^z_j(-q)> and <S^+_i(q) S^-_j(-q)>
     * Each pointer points to a dynamical complex array 2n^2 in size
     * Indexing: corr[i][j][alpha][beta][q] */
    std::complex<double>* S_zz[4][4];
    std::complex<double>* S_pm[4][4];
    
    /* Generates the array index of an arbitrary q-point in the (hhk) plane
     * Unit of q_x,q_z is 2pi/n -> must be integer to be a valid q-point
     * in the periodic box
     * Takes fcc reciprocal lattice into account */
    size_t index(int x, int z) const;

public:
    spin_manager(int n);
    ~spin_manager();

    // Load spin components from file
    void load(std::FILE* f);

    // Add to spin components w times of data from file
    void add(std::FILE* f, double w);

    // Clear contents (fill with zeroes)
    void clear();

    // Evaluate global spin correlators and store in tensor
    void tensor(int qx, int qz, unsigned which = ISING|XY);

    /* NSF and SF neutron scattering response for
     * most recently calculated tensor */
    double neutron_nsf() const;
    double neutron_sf() const;
    
    // Evaluate overall zz correlator
    double Szz(int qx, int qz) const;
};

#endif
