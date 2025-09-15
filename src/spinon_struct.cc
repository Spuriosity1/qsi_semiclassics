/* Introduces a pair of spinons, cools them to ground state in semiclassical
 * QSI, and compares them to quadratic approximation of ground state
 * field pattern
 * Assumes both spinons are on "up" tetrahedra
 *
 * Usage: spinon_struct <paramfile >outfile 2>garbage
 * Paramfile almost of the same structure as in vison_GS
 **  After N, specify spinon charge Q
 * Outfile contains a listing of all pyrochlore sites, one per line:
 ** x y z E_num E_quad
 ** E_num is the numerically found electric field of the hexagon
 * Garbage contains random printouts
 *
 * Same cooling protocol as in vison_GS
 *
 * Created on 03/09/2018
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


#include <debug.hh>
#include <qsi.hh>
#include <spin.hh>
#include <fourier.hh>
#include <direction.hh>
#include <vec3.hh>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_odeiv2.h>

#define TWOPI 6.283185307179586476925286766559005768394338798750211641949

using namespace std;

inline double sqr(double x) {return x*x;}

complex<double> gamma(const vec3& k) {
    return complex<double>(+cos(k[0]/4) * cos(k[1]/4) * cos(k[2]/4),
                           -sin(k[0]/4) * sin(k[1]/4) * sin(k[2]/4));
}

inline complex<double> cis(double x) {return polar(1.0, x);}

int main (int argc, char** argv) {
    // ---------- Reading parameters ----------
    unsigned n;
    double Q;
    vec3_int r1, r2;
    double T_hot, factor;
    unsigned sweep_hot,
        N_step, sweep_step,
        sweep_cold,
        seed;

    cin >> n >> Q >> r1[0] >> r1[1] >> r1[2] >> r2[0] >> r2[1] >> r2[2];
    cin >> T_hot >> sweep_hot;
    cin >> N_step >> factor >> sweep_step;
    cin >> sweep_cold >> seed;

    // --------- Setting up and performing simulation ----------

    fourier f(n,false); 
    qsi simulate(n,&f,NULL,NULL,seed);

    // Create spinons
    /*vec3_int r3(4,4,4);
    simulate.add_spinon_pair(r1,r3);
    simulate.add_spinon_pair(r3,r2);*/
    simulate.add_spinon_pair(r1,r2,Q);
    
    // Cooling protocol
    double T = T_hot;
    simulate.MC(T, sweep_hot, qsi::MC_SMALL, qsi::MC_NONE, qsi::MC_RANDOM);
    db_printerr( "Energy: %.16f\n", simulate.energy()+16*n*n*n);

    for (unsigned i = 0; i < 14; ++i) {
        T /= factor;
        db_printerr( "\nStep %d: T = %.5e\n", i, T); fflush(stderr);
        simulate.MC(T, sweep_step, qsi::MC_SMALL, qsi::MC_NONE, qsi::MC_RANDOM);
        db_printerr( "Energy: %.16f\n", simulate.energy()+16*n*n*n);
        fflush(stderr);
    }
    
    // Set up ODE solver, and move forward in time to spread photons
    /*db_printerr("\nTime evolution\n"); fflush(stderr);
    gsl_odeiv2_system ode = {qsi_evol, NULL, 3*simulate.n_spin(), &simulate};
    gsl_odeiv2_driver *drive =
        gsl_odeiv2_driver_alloc_y_new(&ode, gsl_odeiv2_step_rk8pd,
                                      5e-2, 1e-6, 0.0);
    double time = 0.0;
    gsl_odeiv2_driver_apply(drive, &time, 16.0*n,
                            (double*)simulate.state());
    gsl_odeiv2_driver_free(drive);
    db_printerr( "Energy: %.16f\n", simulate.energy()+16*n*n*n); */

    // Cooling protocol 2
    //T *= pow(factor,4.0);
    //simulate.MC(T, sweep_hot, qsi::MC_SMALL, qsi::MC_NONE);
    for (unsigned i = 14; i < N_step; ++i) {
        T /= factor;
        db_printerr( "\nStep %d: T = %.5e\n", i, T); fflush(stderr);
        simulate.MC(T, sweep_step, qsi::MC_SMALL, qsi::MC_NONE, qsi::MC_RANDOM);
        db_printerr( "Energy: %.16f\n", simulate.energy()+16*n*n*n);
        fflush(stderr);
    }
    
    simulate.MC(T, sweep_cold, qsi::MC_SMALL, qsi::MC_NONE, qsi::MC_RANDOM);

    // Generating the quadratic-level ground state and printing
    for (int sl = 0; sl < 4; ++sl) {
        // Populating the FT array
        for (unsigned x = 0; x < 2*n; ++x)
            for (unsigned y = 0; y < 2*n; ++y)
                for (unsigned z = 0; z < 2*n; ++z) {
                    vec3_int kZ(x,y,z);
                    vec3 k = TWOPI/n * kZ;
                    f(kZ) = 0.25 * Q *
                        (1.0 - cis(k%direction::r[sl] * 2) * conj(gamma(k))) *
                        (cis(k%r1 * -0.125) - cis(k%r2 * -0.125)) /
                        (1.0 - norm(gamma(k)));
                }
        // Correct for the Gamma points
        f(0,0,0) = 0.0;
        f(n,n,n) = 0.0;
        
        f.transform();

        size_t N = 8*n*n*n;

        // Print plaquette data and the appropriate quadratic estimates
        for (size_t i = sl; i < simulate.n_spin(); i+=4) {
            const spin* p = simulate.spin_no(i);
            printf("%4d %4d %4d %19.16f %19.16f\n",
                   p->pos()[0], p->pos()[1], p->pos()[2], p->ising(),
                   f( (p->pos() - direction::pyro[sl]) / 4).real()/N );
        }
    }
    
    return 0;
}
