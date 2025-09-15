/* Introduces 8 visons in a zinc blende motif and finds their ground state in the
 * semiclassical QSI simulation
 *
 * Usage: gs_ZnS paramfile outdir seed 2>garbage
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
// #include <mutex>
#include <filesystem>


typedef basic_parser<int, unsigned, double> param_t;

using namespace std;

inline double sqr(double x) {return x*x;}

int main (int argc, const char** argv) {
    // ---------- Reading parameters ----------
    unsigned n;
    double T_hot, T_cold, factor, gphase, RK_potential, rhog;
    unsigned N_step, sweep_hot, sweep_step, sweep_cold, seed;
    std::string initialiser;
    bool save_visons, save_B, save_spins, silent;

    param_t par("gs_ZnS v1.7");
    par.declare("system_size", &n);
    par.declare("T_hot",&T_hot);
    par.declare("T_cold",&T_cold);
    par.declare("sweep_hot",&sweep_hot);
    par.declare("sweep_cold",&sweep_cold);
    par.declare("n_step",&N_step);
    par.declare("sweep_step",&sweep_step);
    // par.declare("prefix", &prefix);
    
    par.declare("seed", &seed);

    par.declare_optional("phi_g",&gphase, 0);
    par.declare_optional("rho_g",&rhog, 1);
    par.declare_optional("RK_potential",&RK_potential, 0);

    par.declare_optional("save_visons", &save_visons, true);
    par.declare_optional("save_B", &save_B, true);
    par.declare_optional("save_spins", &save_spins, false);
    par.declare_optional("silent", &silent, false);

    par.declare_optional("initialise", &initialiser, "");

    // read it into memory
    par.from_file(argv[1]);
    std::string prefix = std::filesystem::path(argv[1]).filename().replace_extension();
    prefix += par.from_argv(argc, argv, 3);
    par.assert_initialised();

    factor = exp( (log(T_hot) - log(T_cold))/N_step );

    if (argc<3) {
        fprintf(stderr, "Usage: gs_ZnS infile outdir [parameter overrides]\n");
        exit(1);
    }

    if (argv[2][0]=='-') {
        fprintf(stderr, "Usage: gs_ZnS infile outdir [parameter overrides]\n");
        exit(1);
    }

    std::filesystem::path outdir(argv[2]);
    
    

    std::string file; // declare for later

    // --------- Setting up and performing simulation ----------

    fourier f(n, false); // TODO use it for quadratic estimates?
    qsi simulate(n,&f,NULL,NULL,seed);

    //
    simulate.set_uniform_g(polar(1., gphase));
    simulate.set_uniform_RK(RK_potential);

    double T = T_hot;
    const mc_sampler sample_all(mc_sampler::sequential, simulate.n_spin());

    if (initialiser == ""){
        // default: randomise at high temperature
        if (!silent) fprintf(stderr, "Relaxing initial vison configuration at T=%.3f\n",T_hot);
        
        simulate.MC(T, sweep_hot, sample_all, mc_sz_mode::finite_t, mc_angle_mode::finite_t_gradient);
        if (!silent) fprintf(stderr, "Energy: %.16f\n", simulate.energy()+16*n*n*n);
    } else {
        simulate.load_spins(initialiser);
    }
    

    file = outdir / (prefix +  "_order.csv");
    FILE* out = fopen(file.c_str(), "w");

    fprintf(out, "# T order site-energy\n");


    // Anneal from high to low temperature
    for (unsigned i = 0; i < N_step; ++i) {
        T /= factor;
        simulate.MC(T, sweep_step, sample_all, mc_sz_mode::finite_t, mc_angle_mode::finite_t_gradient);
        if (!silent){
            fprintf(stderr, "\nStep %d: T = %.5e\n", i, T);
            fprintf(stderr, "Energy: %.16f\n", simulate.energy()+16*n*n*n);
        }

        // Print the vison order parameter
        fprintf(out, "%+9.6f %+9.6f %+9.6e\n",T, simulate.vison_aiao_order(),simulate.energy()/simulate.n_spin());

    }

    fclose(out);
    
    // Generate cold samples, print energies, take energy stats
    double
        ref = simulate.energy(),
        sum = 0.0,
        sqsum = 0.0;
    for (unsigned i = 0; i < sweep_cold; ++i) {
        simulate.MC(T, 1, sample_all, mc_sz_mode::finite_t, mc_angle_mode::finite_t_von_mises);
        double E = simulate.energy();
        sum += E-ref;
        sqsum += sqr(E-ref);
    }
    db_printerr("\n");
    // sum -> average
    sum /= sweep_cold;
    // sum of squares -> variance
    sqsum = sqsum/sweep_cold - sqr(sum);
    // Save the input file with the same stem
    par.into_file(outdir / (prefix+".inp"));

    if (save_B){
        simulate.save_B(outdir/(prefix+"_B.csv"));
    }

    if (save_visons){
        simulate.save_visons(outdir/(prefix +"_visons.csv"));
    }

    if (save_spins){
        simulate.save_spins(outdir/(prefix+"_spins.csv"));
    }

    return 0;
}
