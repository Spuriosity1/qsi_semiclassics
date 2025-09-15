/* Executable to evaluate static & dynamic structure factors for the QSI model.
 *
 * Usage: qsi_dynamic INFILE OUTPATH
 *
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
#include <misc.hh>
#include <filesystem>

const std::string VERSION = "1.0";

using namespace std;


typedef basic_parser<int, unsigned, double> param_t;

int main (int argc, const char** argv) {
    // ---------- Reading parameters ----------

    double T;

    // Optional parameter variables
    unsigned
        sample = 2048,
        sweep  = 4,
        burnin = 128,        
        n_time = 512, 
        n_freq = 128,
        seed;
    double
        delta = 0.0625;

    bool
	save_couplings = false;
    
    std::string initialiser = "";
    

    param_t p("Dynamical Structure Factors (qsi_dynamic) v" + VERSION + currentDateTime() );
    // p.declare("system_size",&n); This is deduced from the defectfile
    p.declare("temperature",&T);
    p.declare("n_sample", &sample);
    p.declare("n_sweep", &sweep);
    p.declare("n_burnin", &burnin);
    p.declare("n_timesteps", &n_time);
    p.declare("n_Egrid", &n_freq);
    p.declare("dt", &delta);
    p.declare("seed",&seed);

	p.declare_optional("save_couplings", &save_couplings, false);
    
    // the defectfile
    
    // csv file containing a saved spin state for initialisation
    // p.declare_optional("initialise", &initialiser, "");
    // csv file containing a QSI configuration [x][y][z] [phig][RK]
    // p.declare("defectfile", &defectfile);
    
    if (argc < 4){
        std::cout << "Usage: " <<argv[0] << " INFILE DEFECTFILE OUTDIR [INITIALISER] [options...]\n";
        return 1;
    }


    //read the params
    p.from_file(argv[1]);
    std::filesystem::path defectfile(argv[2]);
    std::filesystem::path output_path(argv[3]);

    // Extract the file root
    std::filesystem::path in(argv[1]);

    // Load an initialiser if it is passed
    if (argc >= 5 && argv[4][0]!='-'){
        initialiser = argv[4];
        p.from_argv(argc, argv, 5);
    } else {
        p.from_argv(argc, argv, 4);
    }
    
    std::string prefix = p.get_outfile_prefix(false, false);
    prefix += "%df=";
    prefix += defectfile.filename().replace_extension();
    
    p.assert_initialised();

    //Extract the system size from the defectfilename
    vec3_int L = get_system_size(defectfile.filename().replace_extension());
    my_assert(L.x() == L.y() && L.x() == L.z(), "Only cubes are supported");
    int n = L.x();
    // parameter parsing done
    fourier f(n);
    unsigned what = corrtypes::SZZ | corrtypes::ARG;
    dynamic_correlator dyn(n,n_time,n_freq,what,&f);
    qsi simulate(n,&f,NULL,&dyn,seed);
    
    simulate.set_uniform_g(0);
    simulate.set_uniform_RK(0);
    simulate.load_couplings(defectfile, false);

    mc_sampler sample_all(mc_sampler::sequential, simulate.n_spin());
    

    cerr << "QSI defect dynamic version " << VERSION << "\n";

    // Load the initialiser file
    if (initialiser != ""){
        cerr << "Loading initialiser spins from " << initialiser << "\n";

        simulate.load_spins(initialiser);
        prefix += "%if=";
        prefix += std::filesystem::path(initialiser).filename().replace_extension();
        std::cout <<  std::filesystem::path(initialiser).filename().replace_extension();
    }


    // Burn-in
    fprintf(stderr, "[qsi_dynamic] Burning in...\n");
    simulate.MC(T,burnin, sample_all, mc_sz_mode::finite_t, mc_angle_mode::finite_t_von_mises);
    fprintf(stderr,"\n"); fflush(stderr);

    // Set up ODE solver
    gsl_odeiv2_system ode = {qsi_evol, NULL, 3*simulate.n_spin(), &simulate};
    gsl_odeiv2_driver *drive = 
        gsl_odeiv2_driver_alloc_y_new(&ode, gsl_odeiv2_step_rk8pd,
                                      5e-3, 1e-8, 0.0);




    std::vector<double> GSE(sample);

    // Simulate & record data
    for (unsigned i = 0; i < sample; ++i) {
        fprintf(stderr, "[qsi_dynamic] Simulating: Sample %u/%u\n", i+1,sample);
        // Refresh magnetic noise input
        dyn.clear_input();
        // Monte Carlo
        simulate.MC(T,sweep, sample_all, mc_sz_mode::finite_t, mc_angle_mode::finite_t_von_mises);
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

    
    if (save_couplings){
		std::string out = output_path / prefix;
		std::cout<<"Saving coupling data to"<<out<<std::endl;
		simulate.save_couplings(out + ".couplings.csv");
	}

    // Clean up
    gsl_odeiv2_driver_free(drive);
    return 0;
}
