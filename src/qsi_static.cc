/* Executable to evaluate static structure factors for the QSI model.
 *
 * Usage: qsi_static N T output_root seed [samples] [sweep] [burn_in] [g_phase] > outfile
 ** N:           system size, required
 ** T:           temperature, required
 ** output_root: root of all output filenames, required
 ** seed:        index of seed in the prestored array, required
 ** samples:     number of MC samples, default: 2048
 ** sweep:       number of MC sweeps between samples, default: 4
 ** burn_in:     number of MC sweeps for burning in, default: 128
 ** g_phase:     complex twist to add to g constant, default: 0
 *
 * output_file.* will contain the static spin and B-field correlators
 **      in binary format, as described in static_correlator.hh
 * outfile is a 3-column table, each line gives, for one sample:
 **      #visons,  #nn_dipoles,  energy
 *
 * Created on 01/08/2018
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
#include <corrtypes.hh>
#include <fourier.hh>
#include <static_correlator.hh>
#include <cstdio>
#include <cstdlib>
#include <spin.hh>
#include <string>
#include <basic_parser.hh>

using namespace std;

const std::string VERSION = "1.3";
typedef basic_parser<int, unsigned, double> param_t;

int main (int argc, const char** argv) {
    // ---------- Reading parameters -----------
    // Required parameters
    if (argc < 3){
        fprintf(stderr, "Usage: qsi_static INFILE OUTPATH [overrides....]\n");
        throw std::runtime_error("Required parameters missing");
    }

    int n;
    double T_cold, T_hot = 2;

    // these values are never used in the code, but are good to hold in your mind
    unsigned
        sample = 2048,
        sweep  = 4,
        burnin = 128,        
        n_anneal, seed; 

    
    // std::string initialiser;
    

    param_t p("Static Structure Factors (qsi_static) v" + VERSION + currentDateTime() );
    p.declare("system_size",&n);
    p.declare("temperature",&T_cold);
    p.declare("n_sample", &sample);
    p.declare("n_sweep", &sweep);
    p.declare("n_burnin", &burnin);
    p.declare("seed", &seed);


    double gphase=0, rhog=1;
    p.declare("phi_g", &gphase);
    p.declare_optional("rho_g", &rhog, 1.);

    p.declare("T_hot", &T_hot);
    p.declare("n_anneal", &n_anneal);

    
    
    // Spin component correlators
    bool save_sxx, save_syy, save_szz, save_spm;
    // B correlators
    bool save_sinb, save_imag, save_cosb;
    // saves Arg(ringflip)
    bool save_bb;
    bool save_noise;
    bool save_state;
    bool mc_random_sample;
    
    p.declare_optional("save_sxx", &save_sxx, false);
    p.declare_optional("save_syy", &save_syy, false);
    p.declare_optional("save_szz", &save_szz, true);
    p.declare_optional("save_spm", &save_spm, false);

    p.declare_optional("save_sinb", &save_sinb, false);
    p.declare_optional("save_cosb", &save_cosb, false);
    p.declare_optional("save_imag", &save_imag, false);
    p.declare_optional("save_noise", &save_noise, false);

    p.declare_optional("save_b_state", &save_state, true);

    p.declare_optional("save_bb", &save_bb, true);
    
    p.declare_optional("random_MC_sample", &mc_random_sample, false);

    //read the params
    p.from_file(argv[1]);

    std::filesystem::path output_path(argv[2]);

    // Extract the file root
    std::filesystem::path in(argv[1]);
    std::string prefix(in.filename().replace_extension());
    prefix += p.from_argv(argc, argv, /*start idx of argv*/3);
    p.assert_initialised();

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
    if (save_imag)
        which_correlators |= corrtypes::IMAG;
    if (save_bb)
        which_correlators |= corrtypes::ARG;
    if (save_noise)
        which_correlators |= corrtypes::NOISE;

    fourier f(n);
    static_correlator corr(n,which_correlators,&f);
    qsi simulate(n,&f,&corr,NULL,seed);

    simulate.set_uniform_g(std::polar(rhog, gphase));

    // sample everywhere
    mc_sampler all_samples( mc_random_sample ? mc_sampler::random : mc_sampler::sequential, simulate.n_spin());

    // Burn-in
    simulate.MC(T_hot,burnin, all_samples);

    double T = T_hot;
    double factor = 1;

    if (n_anneal > 0) factor = exp((log(T_hot)-log(T_cold)) / n_anneal);
    fprintf(stderr, "#T     V_t\\eta_t    <E>/N\n");
    // Anneal from high to low temperature
    for (unsigned i = 0; i < n_anneal; ++i) {
        T /= factor;
        simulate.MC(T, sweep, all_samples, mc_sz_mode::finite_t, mc_angle_mode::finite_t_gradient);

        // Print the vison order parameter (redirect this to a file in a bash script if needed)
        fprintf(stderr, "%+9.6f %+9.6f %+9.6e\n",T, simulate.vison_aiao_order(),simulate.energy()/simulate.n_spin());

    }
    
    // Simulate & record data


    std::string metastem;
    metastem = (output_path/prefix).generic_string()+"_energy.csv";
    FILE* metafile = fopen(metastem.c_str(),"w");

    fprintf(metafile, "#n_vison n_dipole E/N\n");
    for (unsigned i = 0; i < sample; ++i) {
        fprintf(stderr, "[ sim ] %u/%u\n",i+1, sample);
        simulate.MC(T,sweep, all_samples, mc_sz_mode::finite_t, mc_angle_mode::finite_t_von_mises);
        simulate.save_static();
        size_t vison, pair;
        simulate.vison_pair(vison,pair);
        fprintf(metafile, "%8lu %8lu %.12e\n", vison, pair, simulate.energy()/simulate.n_spin());
        fflush(stdout);
    }

    // Save data
    std::string outstem(output_path/prefix);
    std::cout<< "Saving to "<< outstem << "\n";
    corr.save_corr(outstem.c_str());

    if (save_state){
        std::string bfile = outstem + ".final_B.csv";
        std::cout<< "Saving B file to "<< bfile << "\n";
        simulate.save_B(bfile);
    }


    return 0;
}


