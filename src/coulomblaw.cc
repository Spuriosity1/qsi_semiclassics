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

typedef std::vector<vec3_int> dirac_string;

bool lost_vison(const qsi& simulate, const std::vector<dirac_string>& vison_locs){
    bool die = false;
    if (simulate.num_visons() != vison_locs.size()*2){
        die = true;
    }
    std::cerr<< "There are " <<simulate.num_visons()<<" visons in the cube, expected "<<vison_locs.size()*2<<"\n";
    for (const auto& vp : vison_locs) {
        std::cerr << vp.front() << " \tsl " << (int)simulate.c_ptetra_at(vp.front())->sublat_enum()
        << "\tcharge[+] " << simulate.vison(vp.front()) << "\n"
        << vp.back() << " \tsl " << (int)simulate.c_ptetra_at(vp.back())->sublat_enum() 
        << "\tcharge[-] " << simulate.vison(vp.back()) << "\n";
        
        die = die || (simulate.vison(vp.front()) != 1);
        die = die || (simulate.vison(vp.back()) != -1);
    }
    return die;
}


std::vector<dirac_string> visons_from_file(const std::string& vison_file, const qsi& simulate, bool strict=true){
    // std::filesystem::path vfile(vison_file);
    // expect format
    //    +        -
    // [0 0 0]  [4 4 4]
    // any amount of whitespace is ignored
    std::vector<dirac_string> vison_list;

    std::ifstream vfin(vison_file, std::ios_base::in);
    if (!vfin) {
        throw std::runtime_error(std::string("Could not open visons file ") + vison_file);
    }
    std::string line;
    unsigned i=0;
    while ( std::getline(vfin,line) ){
        i++;
        // skip comments
        std::istringstream iss(line);
        if ((iss >> std::ws).peek() == '#') continue;
        dirac_string p;
        while((iss>>std::ws).peek() == '['){
            p.push_back(vec3_int::read_from(iss));
        }
        if (p.size() < 2){
            std::cerr<< "Bad line :"<<i<<", not enough visons\n";
            throw std::runtime_error("Vison file parse error: insufficient valid sites");
        }

        for (auto& site : p){
            try {
                ptetra::sublattice sl = simulate.c_ptetra_at(site)->sublat_enum();
                
                std::cerr << "Loaded " << ((sl == ptetra::sublattice::A) ? "(A)" : "(B)");
                if (site == p.front())
                    std::cerr <<  " + vison";
                else if (site == p.back())
                    std::cerr <<  " - vison";
                else
                    std::cerr<<" waypoint";
                    
                std::cerr<<" at " << site <<std::endl;
                
            } catch (const std::range_error& e){
                std::cerr << "Vison location on line "<<i<<" is not on a dual diamond site: "<< site <<"\n" ;
                std::cerr << "Nearest +ve vison site is at " << simulate.nearest_ptetra(site, 0) <<std::endl;
                std::cerr << "Nearest -ve vison site is at " << simulate.nearest_ptetra(site, 1) <<std::endl;
                throw e;
            }
        }

        vison_list.push_back(p);
        
    }

    // Snap to grid
    // for (auto& [v1, v2] : vison_list){
    //     std::cout << "Snapping... "
    //     v1 = simulate.nearest_ptetra(v1, 0);
    //     v2 = simulate.nearest_ptetra(v2, 1);
    // }

    return vison_list;
}


void exclude_from_sampling(mc_sampler& restricted_samples, const qsi& ice, const std::vector<dirac_string> vison_list){
    for ( const auto& vp : vison_list ) {
        
        
        // int sl1 = (simulate.ptetra_sl(v1) == 0) ? 1 : -1;
        // int sl2 = (simulate.ptetra_sl(v2) == 0) ? 1 : -1;

        const ptetra* p1 = ice.c_ptetra_at(vp.front());
        const ptetra* p2 = ice.c_ptetra_at(vp.back());

        for (int mu=0; mu<4; mu++){
            const plaq* pp1 = p1->neighbour(mu);
            const plaq* pp2 = p2->neighbour(mu);

            size_t plaq1_idx = ice.plaq_idx_at( pp1->pos() );
            size_t plaq2_idx = ice.plaq_idx_at( pp2->pos() );
            restricted_samples.ban_plaquette(plaq1_idx);
            restricted_samples.ban_plaquette(plaq2_idx);

            for (int a=0; a<6; a++){
                size_t s1a = ice.spin_idx_at((*pp1)[a].pos());
                size_t s2a = ice.spin_idx_at((*pp2)[a].pos());
                std::cout<<ice.spin_no(s1a)->pos() <<"\t" << ice.spin_no(s2a)->pos()<<"\n";
                restricted_samples.ban_spin(s1a);
                restricted_samples.ban_spin(s2a);
            }
            std::cout<<"\n";
        }
    }
}

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
    double T_hot, T_cold, factor, gphase, rhog, RK_potential, confinement_potential;
    unsigned n_anneal, n_hot, sweep_step, n_cold;
    unsigned n_relax;
    unsigned seed;

    bool save_B, save_visons, save_spins;

    std::string vison_dir, ofile;

    param_t par("Coulomb law v2.3");
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
    par.declare("n_relax",&n_relax);
    par.declare("seed", &seed);

    par.declare("summary_file", &ofile);
    par.declare("vison_lattice_dir", &vison_dir);

    par.declare_optional("phi_g",&gphase, 0);
    par.declare_optional("rho_g",&rhog, 1);
    par.declare_optional("RK_potential",&RK_potential, 0);
    par.declare_optional("confinement_potential",&confinement_potential, 10);

    par.declare_optional("save_B", &save_B, true);
    par.declare_optional("save_visons", &save_visons, false);
    par.declare_optional("save_spins", &save_spins, false);

    
    

    // read it into memory
    // The stem for where to save files
    
    if (argc < 3){
        fprintf(stderr, "Usage: coulomblaw [infile.toml] [output directory] (options...)");
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

    
    // Vison List
    std::filesystem::path vison_file = in.parent_path() / vison_dir;
    vison_file /= string_format("%03d_%03d_%03d.vison", nx, ny, nz);
    
    std::vector<dirac_string> vison_list = visons_from_file(vison_file, simulate, false);
    
    // Ensure that confinement potential has the same sign as g'
    // If g' > 0, positive-sublattice must be swapped due to internal sign convention
    if (sin(gphase) > 0){
        for (auto& v : vison_list){
            std::reverse(v.begin(), v.end());
        }
        confinement_potential =   abs(confinement_potential);
    } else  {
        confinement_potential =  -abs(confinement_potential);
    }

    // Add the visons
    for (const auto& vp : vison_list){
        for (unsigned i=1; i<vp.size(); i++){
            simulate.add_vison_pair(vp[i-1], vp[i],1.0,false);
        }
    }

    // Add a strong onsite potential
    for (const auto& vp : vison_list){
        for (int mu=0; mu<4; mu++){
            simulate.c_ptetra_at(vp.front())->neighbour(mu)->g_constant = std::complex<double>(0,confinement_potential);
            simulate.c_ptetra_at(vp.back())->neighbour(mu)->g_constant = std::complex<double>(0,confinement_potential);
        }
    }

    // Check that visons were successfully inserted
    if (lost_vison(simulate, vison_list)){
        std::cerr <<"Vison introduction failed\n" << outstem+"_B.error.csv" <<"\n" << outstem+"_visons.error.csv\n";
        simulate.save_B(outstem+"_B.error.csv");
        simulate.save_visons(outstem+"_visons.error.csv");
        return 1;
    }
    
    mc_sampler all_samples(mc_sampler::sequential,simulate.n_spin());

    std::cerr<<"Introducing visons.\n";
    std::cerr<<"Hot step annealing...\n";
    
    simulate.MC(T_hot, n_hot, all_samples, mc_sz_mode::finite_t, mc_angle_mode::finite_t_gradient);

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

      // check if the annealing step broke everything... :/
    if (lost_vison(simulate, vison_list)){
        std::cerr << "Hot step has moved visons from their initial positions.\n"<<
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
    
    // Check for lost visons
    if (lost_vison(simulate, vison_list)){
        std::cerr << "Simulation has moved visons from their initial positions.\n"<<
            outstem+"_B.error.csv" <<"\n" << outstem+"_visons.error.csv\n";
        simulate.save_B(outstem+"_B.error.csv");
        simulate.save_visons(outstem+"_visons.error.csv");
        return 3;
    }

    // Remove the onsite potential
    for (const auto& vp : vison_list){
        for (int mu=0; mu<4; mu++){
            simulate.c_ptetra_at(vp.front())->neighbour(mu)->g_constant = std::polar(rhog, gphase);
            simulate.c_ptetra_at(vp.back())->neighbour(mu)->g_constant = std::polar(rhog, gphase);
        }
    }

    
    // Relax
    simulate.MC(T, n_relax, all_samples, mc_sz_mode::finite_t, mc_angle_mode::finite_t_gradient);

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

    // check if the annealing step broke everything... :/
    if (lost_vison(simulate, vison_list)){
        std::cerr << "Cold sample generation has moved visons from their initial positions.\n"<<
            outstem+"_B.error.csv" <<"\n" << outstem+"_visons.error.csv\n";
        simulate.save_B(outstem+"_B.error.csv");
        simulate.save_visons(outstem+"_visons.error.csv");
        return 4;
    }


    
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
    if (save_visons){
        simulate.save_visons(outstem+"_visons.csv");
    }
    if (save_spins){
        simulate.save_spins(outstem+"_spins.csv");
    }
    return 0;
}
