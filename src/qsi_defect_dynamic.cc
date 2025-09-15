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
#include <defect_strategy.hh>
#include <spin.hh>
#include <DO_QSI.hh>
#include <direction.hh>
#include <corrtypes.hh>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <time.h>
#include <misc.hh>
#include <filesystem>

const std::string VERSION = "1.2";

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
        anneal = 128,
        L,
        seed;

    double delta;
    double Jpm;
    double B[3];


    double g[4];
    double g_defect[4];

    double disorder_concentration;
    
    std::string initialiser = "";
    std::string defect_stratname;
    

    param_t p("Dynamical Structure Factors (qsi_dynamic) v" + VERSION + currentDateTime() );
    // p.declare("system_size",&n); This is deduced from the defectfile
    // physical params
    p.declare("disorder_concentration", &disorder_concentration);
    p.declare_optional("defect_strategy", &defect_stratname, "random");

    p.declare("L", &L);
    p.declare("Jpm", &Jpm);

    p.declare_optional("Bx", B  , 0.);
    p.declare_optional("By", B+1, 0.);
    p.declare_optional("Bz", B+2, 0.);

    p.declare("temperature",&T);
    // QSI numerical
    p.declare("n_sample", &sample);
    p.declare("n_anneal", &anneal);
    p.declare("n_sweep", &sweep);
    p.declare("n_burnin", &burnin);
    p.declare("n_timesteps", &n_time);
    p.declare("n_Egrid", &n_freq);
    p.declare("dt", &delta);
    p.declare("seed",&seed);
    
    
    if (argc < 3){
        std::cout << "Usage: " <<argv[0] << " INFILE OUTDIR [options...]\n";
        return 1;
    }


    //read the params
    p.from_file(argv[1]);
    std::filesystem::path output_path(argv[2]);
    p.from_argv(argc, argv, 3);

    // Extract the file root
    std::filesystem::path in(argv[1]);
 
    std::string prefix = p.get_outfile_prefix(false, false);

    std::string out = output_path / prefix;
    
    p.assert_initialised();

    //Extract the system size from the defectfilename
    // parameter parsing done
    fourier f(L);
    unsigned what = corrtypes::SZZ | corrtypes::ARG;
    dynamic_correlator dyn(L,n_time,n_freq,what,&f);
    qsi simulate(L,&f,NULL,&dyn,seed);


    mc_sampler sampling(mc_sampler::sequential, simulate.n_spin());


    cerr << "QSI defect dynamic version " << VERSION << "\n";

    // calculate the g values
    calc_g_values(g, Jpm, B);
    calc_g_defect_values(g_defect, Jpm, B);

    // Gnerate the defects
    std::srand(seed); // randomness less important here

    // set up the uniform g values
    for (unsigned j=0; j<simulate.n_spin(); j++){
        // evil casting away const
        plaq* pl = (plaq *) simulate.plaq_no(j);
        pl->g_constant = std::complex<double>(g[pl->sublat()], 0); 
    }

    // add random nonmagnetic sites (implement this by setting g=0 on their sites)
    std::vector<size_t> defect_spin_idxes;
    std::random_device rd1;
    auto dstrat = defect_strategy(simulate, defect_stratname,
            disorder_concentration, rd1);
    
    
    fprintf(stderr, "[defect] generating defect sites\n");
    std::string defect_list_name = out + ".defects.csv";
    fprintf(stderr, "[defect] saving to %s\n", defect_list_name.c_str()); 
    std::ofstream d_ofs(defect_list_name);
    d_ofs << "spin_idx, X, Y, Z" << std::endl;
    for (auto& site : dstrat.choose_spin_indices()) {
        // evil casting away const
        spin* s = (spin *) simulate.spin_no(site);
        d_ofs << site << ", " << s->pos()[0] << ", " << s->pos()[1] <<", " << s->pos()[2] << std::endl;
        s->set_g_factor(0); // erase the spin from any correlators
        sampling.ban_spin(site); // stop the sampling of the nonexistent spin in MC
                                 // remove contribution of this spin to CMC
        for (unsigned j=0; j<6; j++){
            const plaq* pl = simulate.c_plaq_at(s->pos() + direction::plaqt[s->sublat()][j]);
            plaq* p = (plaq *) pl; // evil casting away const
            p->g_constant = std::complex<double>(g_defect[p->sublat()], 0); 
        }
    }
    d_ofs.close();

    // Run annealing procedure
    //
    // hardcoded max_temp
    double T_hot = (g[0] + g[1] + g[2] + g[3]) /2.;
    double curr_temp = T_hot;
    // T = T_hot * factor^n_anneal => log(factor) = (log(T)-log(T_hot)) /n_anneal 

    fprintf(stderr, "[anneal] Annealing...\n");
    double factor = exp( (log(T)-log(T_hot) ) / anneal );
    for (unsigned x =0; x<anneal; x++) {
        curr_temp *= factor;
        simulate.MC(curr_temp, sweep, sampling,
                mc_sz_mode::finite_t,
                mc_angle_mode::finite_t_von_mises);
    }

    // done 
    auto spinfile_name = out + ".spins";
    auto couplingfile_name = out + ".couplings";
    fprintf(stderr, "[anneal] Anneal complete, state saved to\n");
    fprintf(stderr, "[anneal] %s\n", spinfile_name.c_str());
    fprintf(stderr, "[anneal] %s\n", couplingfile_name.c_str());
    simulate.save_spins(spinfile_name);
    simulate.save_couplings(couplingfile_name);

    // Burn-in
    fprintf(stderr, "[qsi_dynamic] Burning in...\n");
    simulate.MC(T,burnin, sampling, mc_sz_mode::finite_t, mc_angle_mode::finite_t_von_mises);
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
        simulate.MC(T,sweep, sampling, mc_sz_mode::finite_t, mc_angle_mode::finite_t_von_mises);
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
        GSE[i] = 1. + (simulate.energy()/16/L/L/L);

    }
    
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
