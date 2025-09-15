/* Introduces 8 visons in a zinc blende motif and finds their ground state in the
 * semiclassical QSI simulation
 *
 * Usage: rocksalt paramfile outdir seed 2>garbage
 *
 * Contents of paramfile 
 ** system_size = N                           system size
 ** T_hot = T_hot                       initial temperature
 ** sweep_hot = sweep_hot               sweeps at that temperature
 ** N_step = N_step                     number of cooling steps,
 ** T_decay_factor = factor             factor of temp. drop,
 ** sweep_step = nstep                      number of sweeps per step
 ** sweep_cold = ncold                sweeps at lowest T that are output
 *
 * Outfile contains each measured energy point on a separate line,
 * and their average and stdev in the last line.
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
 * Created on 05/09/2018
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
#include <basic_parser.hh>
#include <qsi.hh>
#include <plaq.hh>
#include <spin.hh>
#include <ptetra.hh>
#include <fourier.hh>
#include <vec3.hh>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <direction.hh>
#include "parallel_io.hh"


typedef basic_parser<int, unsigned, double> param_t;

using namespace std;

inline double sqr(double x) {return x*x;}




int main (int argc, char** argv) {
    // ---------- Reading parameters ----------
    unsigned n = atoi(argv[1]);
    double gphase = atof(argv[2]);
    double RK = atof(argv[3]);

    if (n%2 != 0){
        fprintf(stderr, "System size must be even! %d \n",n);
        throw std::runtime_error("odd system size");
    }

    // --------- Setting up and performing simulation ----------
    
    fourier f(n, false); // TODO use it for quadratic estimates?
    qsi simulate(n,&f,NULL,NULL,0);
    
    simulate.set_uniform_g(polar<double>(1,M_PI*gphase));
    simulate.set_uniform_RK(RK);

    std::mt19937 random(6);
    
    std::cout<<simulate.n_spin();
    // Randomise the initial state
    for (size_t j=0; j<simulate.n_spin(); j++){
        std::normal_distribution normal = std::normal_distribution(0.,1.);
        vec3& v = simulate.state()[j];
        v[0] = normal(random);
        v[1] = normal(random);
        v[2] = normal(random);
        
        v.normalise();   
    }

    for (size_t j=0; j<simulate.n_spin(); j++){
        std::complex ee = simulate.spin_no(j)->field_cplx();
        printf("%+6.6f %+6.6f\n",ee, simulate.spin_no(j)->energy());
    }
    
    print_without_collision("test.txt",std::string(argv[2])+"\n");


    return 0;
}
