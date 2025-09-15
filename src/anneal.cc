/* Introduces 8 visons in a zinc blende motif and finds their ground state in the
 * semiclassical QSI simulation
 *
 * Usage: anneal paramfile outdir seed 2>garbage
 * 
 *
 * Outfile contains each measured energy point on a separate line,
 * and their average and stdev in the last line.
 * Garbage contains a lot of typeouts from various routines and some statistics
 * at the very end.
 *
 * Created on 05/09/2018
 * Copyright (C) 2018 Attila Szabó <as2372@cam.ac.uk> 2023 Alaric Sanders <als217@cam.ac.uk>
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
#include "basic_parser.hh"
#include <qsi.hh>
#include <plaq.hh>
#include <spin.hh>
#include <ptetra.hh>
#include <fourier.hh>
#include <vec3.hh>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <direction.hh>

typedef basic_parser<int, unsigned, double> param_t;

using namespace std;


inline double sqr(double x) {return x*x;}


int main (int argc, const char** argv) {
    // ---------- Reading parameters ----------
    unsigned n;
    double gphase, rhog, RK_potential;
    unsigned num_sample, sweep_step;
    unsigned seed;
    std::string initialiser, tempfile, couplingfile;
    bool save_visons, save_B, save_spins, silent;


    if (argc<3 || argv[2][0]=='-' ) {
        fprintf(stderr, "Usage: %s infile outdir [parameter overrides]", argv[0]);
        exit(1);
    }

    param_t par("anneal v1.3");
    par.declare("system_size", &n);
    par.declare("seed", &seed);
    par.declare("num_anneal",&sweep_step); // annealing steps per temperature
    par.declare("num_sample",&num_sample); // sampling steps per temperature
    
    par.declare_optional("phi_g",&gphase, 0);
    par.declare_optional("rho_g", &rhog, 1);
    par.declare_optional("RK_potential",&RK_potential, 0);

    par.declare("temperature_file", &tempfile);
    

    par.declare_optional("save_visons", &save_visons, true);
    par.declare_optional("save_B", &save_B, true);
    par.declare_optional("save_spins", &save_spins, false);
    par.declare_optional("silent", &silent, false);

    par.declare_optional("initial_spins", &initialiser, "");


    std::filesystem::path output_path(argv[2]);

    // Extract the file root
    std::filesystem::path in(argv[1]);
    std::string prefix(in.filename().replace_extension());
    // read it into memory
    par.from_file(argv[1]);
    prefix += par.from_argv(argc, argv, /*start idx of argv*/3);
    par.assert_initialised();


    // --------- Setting up and performing simulation ----------

    fourier f(n, false); // TODO use it for quadratic estimates?
    qsi simulate(n,&f,NULL,NULL,seed);
    mc_sampler sampling(mc_sampler::sequential, simulate.n_spin());

    // Load the parameters for uniform simulation
    simulate.set_uniform_g(polar(rhog, gphase));
    simulate.set_uniform_RK(RK_potential);

    if (initialiser != ""){
        simulate.load_spins(initialiser);   
    }
    
    if (couplingfile != ""){
        simulate.load_couplings(couplingfile);
    }
    

    // Set up temperature grid
    std::vector<std::pair<bool, double> > temp_grid;
    
    // load temperatures from file
    std::filesystem::path tfpath = in.replace_filename(tempfile);
    ifstream ifs(tfpath);
    if (!ifs.good()){
        std::cerr << "Cannot open temperature file " << tfpath << std::endl;
        throw std::runtime_error("Could not open temperature file");
    }

    while (ifs){
        double T;
        char c;
        ifs >> c >> T;
        if (c == 'S') {
            temp_grid.push_back(std::pair(true, T));
        } else if (c == 'A') {
            temp_grid.push_back(std::pair(false, T));
        } else {
            // assume we wanted to save it
            std::cerr<< "[ WARN ] All tempfile lines should start with S or A";
            std::cerr<< "for (S)ample or (A)nneal. Only S lines are reported.\n";
            temp_grid.push_back(std::pair(true, T));
        }

    }
    

    std::string outstem (output_path / prefix);

    

	FILE* out = std::fopen((outstem + "_order.csv").c_str(), "w");


	fprintf(out, "#T sumE sumE2 sumE3 sumE4 sumV sumV2 sumS sumS2 sumCosb sumSinb sumCos2b sumSin2b\n");

    for (auto& [save_step, T] : temp_grid) {
        
        // Equlibrate to the new temperatue
        // TODO is there a way to be sure that this worked?
        simulate.MC(T, sweep_step, sampling);

        if (!silent){
            fprintf(stderr, "\nT = %.5e\n Energy: %.16f\n", T, simulate.energy());
        }
        
        if (!save_step) continue;
        
        double sum_energy[4] = {0,0,0,0}; 
        double sum_vison[2] = {0,0}; 
        double sum_orderparam[2] = {0,0};
        double sum_spinon[2] = {0,0};
        double sum_sinb[2] = {0,0};
        double sum_cosb[2] = {0,0};

        double tmp;
        std::complex<double> tmp2;
        for (unsigned i=0; i<num_sample; i++){
            tmp = simulate.energy();
            sum_energy[0]  += tmp;
            sum_energy[1] += tmp*tmp;
            sum_energy[2] += tmp*tmp*tmp;
            sum_energy[3] += tmp*tmp*tmp*tmp;

            tmp = simulate.num_visons();
            sum_vison[0]  += tmp;
            sum_vison[1] += tmp*tmp;

            tmp = simulate.vison_aiao_order();
            sum_orderparam[0]  += tmp;
            sum_orderparam[1] += tmp*tmp;

            tmp = simulate.spin_aiao_order();
            sum_spinon[0]  += tmp;
            sum_spinon[1] += tmp*tmp;

            // average value of the Wilson loop
            tmp2 = simulate.avg_cisb();
            sum_cosb[0] += std::real(tmp2);
            sum_cosb[1] += std::real(tmp2)*std::real(tmp2);
            sum_sinb[0] += std::imag(tmp2);
            sum_sinb[1] += std::imag(tmp2)* std::imag(tmp2);

            simulate.MC(T, 1, sampling);
        }
 
		fprintf(out, "%+9.6e ", T);
		for (int i=0; i<4; i++){
			fprintf(out, "%+9.6e ", sum_energy[i] / simulate.n_spin());
		}

		for (int i=0; i<2; i++){
			fprintf(out, "%+9.6e ", sum_vison[i] / simulate.n_tetra());
		}

		for (int i=0; i<2; i++){
			fprintf(out, "%+9.6e ", sum_spinon[i] / simulate.n_tetra());
		}

		for (int i=0; i<2; i++){
			fprintf(out, "%+9.6e ", sum_cosb[i]);
			fprintf(out, "%+9.6e ", sum_sinb[i]);
		}

		fprintf(out, "\n");

	}

	fclose(out); 
    
    
    if (save_B){
        simulate.save_B(outstem+"_B.csv");
    }

    if (save_visons){
        simulate.save_visons(outstem+"_visons.csv");
    }

    if (save_spins){
        simulate.save_spins(outstem+"_spins.csv");
    }

    return 0;
}
