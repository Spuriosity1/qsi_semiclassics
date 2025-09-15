/* Executable to evaluate static & dynamic structure factors for the QSI model.
 *
 * Usage: qsi_dynamic INFILE OUTPATH SEED
 *
 * SEED: RNG seed
 * 
 * INFILE: simulator settings, with structure
 ** label:     name for the parameter set. Choose it wisely!
 ** N:         system size, required
 ** T:         temperature, required
 ** n_sample:  number of MC samples, default: 2048
 ** sweep:     number of MC sweeps between samples, default: 4
 ** n_burnin:  number of MC sweeps for burning in, default: 128
 ** n_time:    number of time steps on each MC sample, default: 512
 ** n_freq:    number of energy/frequency points to output per point on the hardcoded dispersion, default: 128
 ** dt:        duration of time step [1/g], default: 0.0625
 ** theta:     complex 'twist' to add to g
 *
 * OUTPATH: output folder (Output files are named [label]+[SEED].szz.dyn, [label]+[SEED].noise.dyn )
 **      
 * 
 * Created on 06/08/2018
 * Copyright (C) 2018,2019 Attila Szabó <as2372@cam.ac.uk>
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

#include <basic_parser.hh>
#include <debug.hh>
#include <spin.hh>
#include <plaq.hh>
#include <qsi.hh>
#include <fourier.hh>
#include <dynamic_correlator.hh>
#include <corrtypes.hh>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <time.h>
#include <filesystem>

const std::string VERSION = "1.9";

using namespace std;


typedef basic_parser<int, unsigned, double> param_t;


int main (int argc, const char** argv) {
    // ---------- Reading parameters ----------
    // Required parameters
    if (argc < 3){
        fprintf(stderr, "Usage: qsi_dynamic INFILE OUTPATH [overrides....]\n");
        throw std::runtime_error("Required parameters missing");
    }

    int n;
    double T;

    // parameter variables
    unsigned sample, sweep, burnin, n_time, n_freq, seed;
    double delta;
    
    std::string initialiser;
    

    param_t p("Dynamical Structure Factors (qsi_dynamic) v" + VERSION + currentDateTime() );
    p.declare("system_size",&n);
    p.declare("temperature",&T);
    p.declare("n_sample", &sample);
    p.declare("n_sweep", &sweep);
    p.declare("n_burnin", &burnin);
    p.declare("n_timesteps", &n_time);
    p.declare("n_Egrid", &n_freq);
    p.declare("dt", &delta);
    p.declare("seed", &seed);

    double gphase, rhog, RK_potential;
    p.declare_optional("phi_g", &gphase, 0);
    p.declare_optional("rho_g", &rhog, 1);
    // p.declare_optional("Ising_3nn", &Ising_3nn, 0.);
    p.declare_optional("RK_potential", &RK_potential, 0);

    p.declare_optional("initialise", &initialiser, "");

    // Spin component correlators
    bool save_sxx, save_syy, save_szz, save_spm;
    // B correlators
    // these two are redundant
    bool save_sinb, save_cosb;
    // saves Arg(ringflip)
    bool save_bb;
    bool save_noise;
    
    p.declare_optional("save_sxx", &save_sxx, false);
    p.declare_optional("save_syy", &save_syy, false);
    p.declare_optional("save_szz", &save_szz, true);
    p.declare_optional("save_spm", &save_spm, false);

    p.declare_optional("save_sinb", &save_sinb, false);
    p.declare_optional("save_cosb", &save_cosb, false);
    p.declare_optional("save_noise", &save_noise, false);

    p.declare_optional("save_bb", &save_bb, true);
    
    

    //read the params
    p.from_file(argv[1]);

    std::filesystem::path output_path(argv[2]);
    

    // Extract the file root
    std::filesystem::path in(argv[1]);
    std::string prefix(in.filename().replace_extension());
    prefix += p.from_argv(argc, argv, 3);
    p.assert_initialised();

    // parameter parsing done

    // Decide which correlators to store
    unsigned which_correlators = 0;
    if (save_sxx)
        which_correlators |= corrtypes::SXX;
    if (save_syy)
        which_correlators |= corrtypes::SYY;
    if (save_szz)
        which_correlators |= corrtypes::SZZ;
    if (save_spm)
        which_correlators |= corrtypes::SPM;
    if (save_sinb)
        which_correlators |= corrtypes::SIN;
    if (save_cosb)
        which_correlators |= corrtypes::COS;
    // if (save_imag)
    //     which_correlators |= corrtypes::IMAG;
    if (save_bb)
        which_correlators |= corrtypes::ARG;
    if (save_noise)
        which_correlators |= corrtypes::NOISE;

    fourier f(n);
    dynamic_correlator dyn(n,n_time,n_freq,which_correlators,&f);
    qsi simulate(n,&f,NULL,&dyn,seed);

    mc_sampler sample_all(mc_sampler::sequential, simulate.n_spin());

    simulate.set_uniform_g(std::polar(rhog, gphase));
    simulate.set_uniform_RK(RK_potential);
    // simulate.set_uniform_Ising_3nn(Ising_3nn);

    cerr << "QSI dynamic version " << VERSION << "\n";

    // Load the initialiser file
    if (initialiser != ""){
        simulate.load_spins(initialiser);
    }

    // Burn-in
    fprintf(stderr, "[qsi_dynamic] Burning in...\n");
    simulate.MC(T, burnin, sample_all);
    fprintf(stderr,"\n"); fflush(stderr);

    // Set up ODE solver
    gsl_odeiv2_system ode = {qsi_evol, NULL, 3*simulate.n_spin(), &simulate};
    gsl_odeiv2_driver *drive = 
        gsl_odeiv2_driver_alloc_y_new(&ode, gsl_odeiv2_step_rk8pd, 5e-3, 1e-8, 0.0);

    std::vector<double> GSE(sample);

    // Simulate & record data
    for (unsigned i = 0; i < sample; ++i) {
        fprintf(stderr, "[qsi_dynamic] Simulating: Sample %u/%u\n", i+1,sample);fflush(stderr);
        // Refresh magnetic noise input
        dyn.clear_input();
        // Monte Carlo
        simulate.MC(T,sweep, sample_all);
        simulate.save_dynamic(0);
        double t = 0.0;

        // Evolve ODE (using state vector as evolution vector - should be fine)
        for (unsigned zt = 1; zt < n_time; ++zt) {
            int status = gsl_odeiv2_driver_apply(drive, &t, delta*zt,
                                                 (double*)simulate.state());
            if (status != GSL_SUCCESS)
                throw status;
            simulate.save_dynamic(zt);
            db_printerr("%1d",zt%10); fflush(stderr);
        }
        dyn.add_correlator();
        // save the ground state energy of the sytem
        GSE[i] = 1. + (simulate.energy()/16/n/n/n);

    }
    
    std::string out = output_path / prefix;
    std::cout<<"Saving binary data to"<<out<<std::endl;
    dyn.save_corr(out.c_str());

    // save the GSE vector
    ofstream ofs(out + ".gse.csv");
    for (unsigned i = 0; i< sample; i++ ){
        ofs << GSE[i] << std::endl;
    }
    ofs.close();

    

    // Clean up
    gsl_odeiv2_driver_free(drive);
    return 0;
}
