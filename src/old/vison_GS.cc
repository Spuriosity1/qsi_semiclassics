/* Executable to introduce visons and cool them to ground state in
 * semiclassical QSI simulations.
 *
 * Usage: vison_GS  <paramfile >outfile 2>garbage
 *
 * Contents of paramfile (layout doesn't matter, order does)
 ** N                           system size
 ** x1 y1 z1 x2 y2 z2           coordinates of visons
 ** T_hot sweep_hot             initial temperature, sweeps at that temperature
 ** N_step factor sweep_step    number of cooling steps, factor of temp. drop,
 **                             number of sweeps per step
 ** sweep_cold                  sweeps at lowest T that are output
 ** seed                        ID of seed in the seed file
 *
 * Outfile contains each measured energy point on a separate line
 * Garbage contains a lot of typeouts from various routines and some statistics
 * at the very end.
 *
 * Cooling protocol:
 **    sweep_hot sweeps at T_hot
 **    N_step steps
 ***       cooling by given factor
 ***       sweep_step sweeps at that temperature
 **    sweep_cold sweeps, energy of each listed in outfile
 *
 * Created on 23/08/2018
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

#include <qsi.hh>
#include <fourier.hh>
#include <vec3.hh>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <debug.hh>

using namespace std;

inline double sqr(double x) {return x*x;}

int main (int argc, char** argv) {
    // ---------- Reading parameters ----------
    unsigned n;
    vec3_int r1, r2;
    double T_hot, factor;
    unsigned sweep_hot,
        N_step, sweep_step,
        sweep_cold,
        seed;

    cin >> n >> r1[0] >> r1[1] >> r1[2] >> r2[0] >> r2[1] >> r2[2];
    cin >> T_hot >> sweep_hot;
    cin >> N_step >> factor >> sweep_step;
    cin >> sweep_cold >> seed;

    // --------- Setting up and performing simulation ----------
    
    fourier f(n); // FFS we need to fix this
    qsi simulate(n,&f,NULL,NULL,seed);

    // Create visons
    simulate.add_vison_pair(r1,r2);

    // Cooling protocol
    double T = T_hot;
    simulate.MC(T, sweep_hot, qsi::MC_SMALL, qsi::MC_SMALL);
    if ((simulate.vison(r1) != 1) || (simulate.vison(r2) != -1))
        throw "Lost visons!";
    db_printerr( "Energy: %.16f\n", simulate.energy()+16*n*n*n);

    for (unsigned i = 0; i < N_step; ++i) {
        T /= factor;
        db_printerr( "\nStep %d: T = %.5e\n", i, T); fflush(stderr);
        simulate.MC(T, sweep_step, qsi::MC_SMALL, qsi::MC_SMALL);
        db_printerr( "Energy: %.16f\n", simulate.energy()+16*n*n*n);
        fflush(stderr);
    }

    // Generate cold samples, print energies, take energy stats
    double
        ref = simulate.energy(),
        sum = 0.0,
        sqsum = 0.0;
    for (unsigned i = 0; i < sweep_cold; ++i) {
        simulate.MC(T, 1, qsi::MC_SMALL, qsi::MC_SMALL);
        double E = simulate.energy();
        sum += E-ref;
        sqsum += sqr(E-ref);
        printf("%.16f\n", 16*n*n*n + E);
    }
    db_printerr("\n");
    // sum -> average
    sum /= sweep_cold;
    // sum of squares -> variance
    sqsum = sqsum/sweep_cold - sqr(sum);
    
    // Print mean energy, stdev, exptal error, heat capacity (is this useful?)
    db_printerr("%4d %4d %4lu %.16f %.6e %.6e %.6e\n",
            simulate.vison(r1), simulate.vison(r2), simulate.num_visons(),
            ref + sum + 16*n*n*n, sqrt(sqsum), sqrt(sqsum/sweep_cold),
            sqsum/sqr(T));
    printf("%.16f %.6e\n", ref+sum+16*n*n*n, sqrt(sqsum));

    return 0;
}
