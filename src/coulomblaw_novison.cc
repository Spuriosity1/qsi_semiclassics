/* Introduces 8 visons in a rocksalt motif and finds their ground state in the
 * semiclassical QSI simulation
 *
 * Usage: coulomblaw paramfile seed >outfile 2>garbage
 *
 * Contents of paramfile 
 ** system_size = N                           system size
 ** T_hot = T_hot                       initial temperature
 ** n_hot = n_hot               sweeps at that temperature
 ** n_anneal = n_anneal                     number of cooling steps,
 ** T_decay_factor = factor             factor of temp. drop,
 ** sweep_step = nstep                      number of sweeps per step
 ** n_cold = ncold                sweeps at lowest T that are output
 *
 * Outfile contains each measured energy point on a separate line,
 * and their average and stdev in the last line.
 * Garbage contains a lot of typeouts from various routines and some statistics
 * at the very end.
 *
 * Cooling protocol:
 **    n_hot sweeps at T_hot
 **    n_anneal steps
 ***       cooling by given factor
 ***       sweep_step sweeps at that temperature
 **    n_cold sweeps, energy of each listed in outfile
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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>


#include <debug.hh>
#include <basic_parser.hh>
#include <qsi.hh>
#include <plaq.hh>
#include <ptetra.hh>
#include <fourier.hh>
#include <vec3.hh>
#include <direction.hh>
#include "parallel_io.hh"
#include "spin.hh"

typedef basic_parser<int, unsigned, double> param_t;

using namespace std;

inline double sqr(double x) {return x*x;}


vec3 b_sum(const qsi& q){
    vec3 total_b(0,0,0);
    for (unsigned i=0; i<q.n_spin(); i++){
        const plaq* p = q.plaq_no(i);
        total_b += p->B()*direction::r[p->sublat()]*8/std::sqrt(3);
    }
    return total_b;
}


int main (int argc, const char** argv) {
    // ---------- Reading parameters ----------
    unsigned nx, ny, nz, L;
    // double T_temper;
    double T_hot, T_cold, factor, gphase, rhog, RK_potential;
    unsigned n_anneal, n_hot, sweep_step, n_cold;
    unsigned seed;

    bool save_B, save_spins;

    std::string ofile;

    param_t par("Coulomb law novison v1.0");
    par.declare_optional("Nx", &nx, 0);
    par.declare_optional("Ny", &ny, 0);
    par.declare_optional("Nz", &nz, 0);
    par.declare_optional("L", &L, 0);
    // par.declare("T_temper",&T_temper);
    // par.declare("n_temper",&n_temper);
    par.declare("T_hot",&T_hot);
    par.declare("n_hot",&n_hot);
    par.declare("n_anneal",&n_anneal);
    par.declare("n_cold",&n_cold);
    par.declare_optional("T_cold",&T_cold, 1e-4);
    par.declare("sweep_anneal",&sweep_step);
    par.declare("seed", &seed);

    par.declare("summary_file", &ofile);

    par.declare_optional("phi_g",&gphase, 0);
    par.declare_optional("rho_g",&rhog, 1);
    par.declare_optional("RK_potential",&RK_potential, 0);

    par.declare_optional("save_B", &save_B, true);
    par.declare_optional("save_spins", &save_spins, false);

    
    

    // read it into memory
    // The stem for where to save files
    
    if (argc < 3){
        fprintf(stderr, "Usage: coulomblaw_novison [infile.toml] [output directory] (options...)");
    }

    std::filesystem::path in(argv[1]);

    par.from_file(in.c_str());
    std::filesystem::path output_path(argv[2]);

    // the filename (different quantitties will be appended to this)
    std::string prefix(in.filename().replace_extension());
    prefix += par.from_argv(argc, argv, 3);
    par.assert_initialised();

    // make sure that exactly one of L or Nx,Ny,Nz are declared
    bool valid_L = (L!=0);
    bool valid_Lxyz = (nx !=0 && ny != 0 && nz != 0);
    if (valid_L ^ valid_Lxyz ){
        if (valid_L) {
            nx = L;
            ny = L;
            nz = L;
        }
    } else {
        std::cerr<<"Must specify exactly one of L or (Nx, Ny, Nz)\n";
        throw std::runtime_error("Bad system size specification");
    }

    // sanitise prefix
    std::replace( prefix.begin(), prefix.end(), '/', '_');
    std::replace( prefix.begin(), prefix.end(), '\\', '_');
    std::replace( prefix.begin(), prefix.end(), ' ', '_');

    // Slash operator will handle DOS paths correctly
    std::string outstem(output_path / prefix);
        

        
    // --------- Setting up and performing simulation ----------
    unsigned nmax = std::max({nx, ny, nz});
    // fourier f(nmax, false); // TODO use it for quadratic estimates?
    qsi simulate(nx, ny, nz, NULL,NULL,NULL,seed);

    simulate.set_uniform_g(polar(rhog, gphase));
    simulate.set_uniform_RK(RK_potential);
    
    
    

    // Add in a nontrivial loop in the z direction
    vec3_int x0 = simulate.nearest_ptetra(vec3_int(nx*4,ny*4,0), 0);
    vec3_int x1 = simulate.nearest_ptetra(vec3_int(nx*4,ny*4,nz*8/3), 0);
    vec3_int x2 = simulate.nearest_ptetra(vec3_int(nx*4,ny*4,nz*8*2/3), 0);


    
    simulate.add_vison_pair(x0, x1,1.0,false);
    simulate.add_vison_pair(x1, x2,1.0,false);
    simulate.add_vison_pair(x2, x0,1.0,false);
    
    
    mc_sampler all_samples(mc_sampler::sequential,simulate.n_spin());

    std::cerr<<"Introducing visons.\n";
    std::cerr<<"Hot step annealing...\n";
    
    simulate.MC(T_hot, n_hot, all_samples, mc_sz_mode::finite_t, mc_angle_mode::finite_t_gradient);

    // make sure it worked
    if (simulate.num_visons() != 0){
        std::cerr<<"Simulation added visons\n";
        return 1;
    }
    if (b_sum(simulate).len2() < 1e-3 ){
        std::cerr<<"Simulation failed to polarise the system\n";
        return 2;
    }

    fprintf(stderr, "Energy: %.16f\n", simulate.energy()+16*nx*ny*nz);
    double T = T_hot;

    std::cerr<<"Total magnetic fields are "<< b_sum(simulate) <<"\n";
    if (b_sum(simulate).len() > 1e-10){
        std::cerr<<"Visons introduced with finite polarisation, which will hamper any attempts to anneal.\n" <<
            outstem+"_B.error.csv" <<"\n" << outstem+"_visons.error.csv\n";
        simulate.save_B(outstem+"_B.error.csv");
        simulate.save_visons(outstem+"_visons.error.csv");
        return 2;
    }

      

    

    // T_cold = T_hot * (factor)^-n_anneal
    factor = std::exp(-std::log(T_hot/T_cold) / static_cast<double>(n_anneal));
    // Anneal from high to low temperature
    for (unsigned i = 0; i < n_anneal; ++i) {
        T *= factor;
        fprintf(stderr, "\nStep %d: T = %.5e\n", i, T); fflush(stderr);
        simulate.MC(T, sweep_step, all_samples, mc_sz_mode::finite_t, mc_angle_mode::finite_t_gradient);
        fprintf(stderr, "Energy: %.16f\n", simulate.energy()+16*nx*ny*nz);
        fflush(stderr);
    }

    // make sure it worked
    if (simulate.num_visons() != 0){
        std::cerr<<"Simulation added visons\n";
        return 3;
    }
    if (b_sum(simulate).len2() < 1e-3 ){
        std::cerr<<"Annealing removed polarisation from the system\n";
        return 4;
    }


    // Generate cold samples, print energies, take energy stats
    double
        ref = simulate.energy(),
        sum = 0.0,
        sqsum = 0.0;
    for (unsigned i = 0; i < n_cold; ++i) {
        simulate.MC(T, 1, all_samples, mc_sz_mode::finite_t, mc_angle_mode::finite_t_von_mises);
        double E = simulate.energy();
        sum += E-ref;
        sqsum += sqr(E-ref);
        // printf("%.16f\n", 16*n*n*n + E);
    }
    fprintf(stderr, "\n");
    // sum -> average
    sum /= n_cold;
    // sum of squares -> variance
    sqsum = sqsum/n_cold - sqr(sum);


    
    std::filesystem::path summary_out = output_path / ofile;

    // File Format:                     seed, rho_g, phi_g, Nx, Ny, Nz, mean energy, var
    std::string summary = string_format("%8u %.16f %.16f %4u %4u %4u %.16f %.16e\n", 
                                        seed, rhog, gphase, nx, ny, nz, ref + sum + 16*nx*ny*nz, sqsum);
    std::cout << summary;
    print_without_collision(summary_out, summary);

    par.into_file(outstem+".inp");



    // auxiliary info
    if (save_B){
        simulate.save_B(outstem+"_B.csv");
    }
    if (save_spins){
        simulate.save_spins(outstem+"_spins.csv");
    }
    return 0;
}
