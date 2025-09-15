/* Executable to evaluate heat current fluctuations for the QSI model.
 *
 * Usage: qsi_thermal_cond <N> <T> <file> <seed> <burn_MC> <burn_t> <step> <dt>
 *
 ** N:       system size
 ** T:       temperature
 ** file:    output filename
 ** seed:    index of seed in prestored random array (0..16ki)
 ** burn_MC: number of burn-in MC steps
 ** burn_t:  burn-in time [1/g]
 ** step:    number of tiem steps
 ** dt:      duration of time step [1/g]
 *
 * Created on 20/06/2019
 * Copyright (C) 2019 Attila Szabó <as2372@cam.ac.uk>
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
#include <plaq.hh>
#include <fourier.hh>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

using namespace std;

int main (int argc, char** argv) {
    // ---------- Reading parameters ----------
    if (argc < 9)
        throw std::runtime_error("Required parameters missing");
    
    int n = atoi(argv[1]);
    double T = atof(argv[2]);
    size_t seed = atoi(argv[4]);
    size_t burn_MC = atoi(argv[5]);
    double burn_t = atof(argv[6]);
    size_t step = atoi(argv[7]);
    double dt = atof(argv[8]);

    // ---------- Simulation setup ----------
    fourier f(n);
    qsi simulate(n, &f, NULL, NULL, seed);

    // Burn-in MC steps
    simulate.MC(T,burn_MC);
    db_printerr("MC done\n"); fflush(stderr);

    // Set up ODE solver
    gsl_odeiv2_system ode = {qsi_evol, NULL, 3*simulate.n_spin(), &simulate};
    gsl_odeiv2_driver *drive = 
        gsl_odeiv2_driver_alloc_y_new(&ode, gsl_odeiv2_step_rk8pd,
                                      5e-3, 1e-8, 0.0);
    
    // Burn-in time evolution
    double t = 0.0;
    if (burn_t > 0.0) {
        t = -burn_t;
        int status = gsl_odeiv2_driver_apply(drive, &t, 0.0,
                                             (double*)simulate.state());
        if (status != GSL_SUCCESS)
            throw status;
        db_printerr("Burn-in time done\n"); fflush(stderr);
    }

    // Simulate & record data
    FILE* fout = fopen(argv[3], "wb");
    //double B_sl[4]; // magnetic field summed over sublattices
    for (size_t i = 0; i < step; ++i) {
        int status = gsl_odeiv2_driver_apply(drive, &t, dt*i,
                                             (double*)simulate.state());
        if (status != GSL_SUCCESS)
            throw status;
        vec3 current = simulate.heat_current();
        /*B_sl[0] = B_sl[1] = B_sl[2] = B_sl[3] = 0;
        for (size_t pl = 0; pl < simulate.n_spin(); ++pl)
        B_sl[pl%4] += simulate.plaq_no(pl)->B();*/
        fwrite(&current, sizeof(vec3), 1, fout);
        //fwrite(B_sl, 4*sizeof(double), 1, fout);
        db_printerr( "%1lu", i%10); fflush(stderr);
    }
    db_printerr("\n"); fflush(stderr);

    // Clean up
    gsl_odeiv2_driver_free(drive);
    fclose(fout);
    return 0;
}
