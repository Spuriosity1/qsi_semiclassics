#include <basic_parser.hh>
#include <cstdlib>
#include <qsi.hh>
#include <ptetra.hh>
#include <plaq.hh>
#include <spin.hh>
#include <stdexcept>
#include <string>
#include <iostream>
#include <direction.hh>
#include <defect_strategy.hh>
#include <DO_QSI.hh>

typedef basic_parser<int, unsigned, double> parser_t;

/*
// calculates g0,1,2,3 from given B and J_\pm
// Uses the expression given in Equation (4) of
// arXiv:2312.11641v1
// Gives the thing in units of 1.5 Jpm^3
void calc_g_values(double* g, double Jpm, const double B[]){
    // safe code is for cowards
    for (uint8_t mu=0; mu<4; mu++){
        double z_dot_B = -1* (
            direction::axis[mu][2][0]*B[0]+
            direction::axis[mu][2][1]*B[1]+
            direction::axis[mu][2][2]*B[2]
        );
        
        g[mu] = 1.5*Jpm*Jpm*Jpm + 1.25*Jpm*Jpm*z_dot_B*z_dot_B;
        g[mu] /= fabs(1.5*Jpm*Jpm*Jpm);
    }
}

void calc_g_defect_values(double* g_defect, double Jpm, const double B[]){
    for (uint8_t mu=0; mu<4; mu++){
        double z_dot_B = -1*(
            direction::axis[mu][2][0]*B[0]+
            direction::axis[mu][2][1]*B[1]+
            direction::axis[mu][2][2]*B[2]
        );
        
        g_defect[mu] = -8.*Jpm*Jpm*z_dot_B;
        g_defect[mu] /= fabs(1.5*Jpm*Jpm*Jpm);
    }
}
*/


/**
 * @brief Generates a DO-QSI with uniform plaquette values.
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */
int main(int argc, const char* argv[]){
    parser_t params;
    // unsigned Lx, Ly, Lz;
    unsigned L;

    double B[3]; // units of J_yy
    double Jpm;  // units of J_yy

    double g[4];
    double g_defect[4];

    unsigned int num_defects, seed;
    std::string defect_strategy;
    double defect_concentration;
    std::string ofile = "defect_doqsi";

    // params.declare("Lx", &Lx);
    // params.declare("Ly", &Ly);
    // params.declare("Lz", &Lz);

    params.declare("L", &L);

    params.declare("num_defects", &num_defects);
    params.declare_optional("defect_strategy", &defect_strategy, "random");
    params.declare("defect_concentration", &defect_concentration);

    params.declare("Bx", B  );
    params.declare("By", B+1);
    params.declare("Bz", B+2);

    params.declare("Jpm", &Jpm);
    params.declare("seed", &seed);




    if (argc < 2){
        std::cerr << "Usage: generate_couplings OUTDIR [options...]" << std::endl;
        return 1;
    }
    ofile += params.from_argv(argc, argv, 2);
    params.assert_initialised();

    if(Jpm == 0){
        std::cerr << "Cannot have zero Jpm"<<std::endl;
        return 1;
    }

    // calculate the uniform g values
    calc_g_values(g, Jpm, B);
    calc_g_defect_values(g_defect, Jpm, B);

    std::filesystem::path outdir(argv[1]);

    qsi simulate(L, L, L, NULL, NULL, NULL);
    std::srand(seed);

    // set up the uniform g values
    for (unsigned j=0; j<simulate.n_spin(); j++){
        // evil casting away const
        plaq* p = (plaq *) simulate.plaq_no(j);
        p->g_constant = std::complex<double>(g[p->sublat()], 0); 
    }

    // add random nonmagnetic sites (implement this by setting g=0 on their sites)
    std::vector<size_t> defect_spin_idxes;
    std::random_device rd2;
    auto dstrat = defect_strategy(simulate, "random", defect_concentration, rd2);
    
    
    for (auto& site : dstrat.choose_spin_indices()) {
        // evil casting away const
        spin* s = (spin *) simulate.spin_no(site);
        s->set_g_factor(0); // erase the spin from any correlators
        for (unsigned j=0; j<6; j++){
            const plaq* pl = simulate.c_plaq_at(s->pos() + direction::plaqt[s->sublat()][j]);
            plaq* p = (plaq *) pl; // evil casting away const
            p->g_constant = std::complex<double>(g_defect[p->sublat()], 0); 
        }
    }

    // std::string pref = couplingfile_name+params.get_outfile_prefix(false, false);] 
    
    std::cerr << "Saving couplings...\n";
    const auto& fpath = outdir / (ofile + ".couplings");
    simulate.save_couplings(fpath);
    std::cout << std::string(fpath) << "\n";

    const auto& spath = outdir / (ofile + ".spins");
    std::cerr << "Saving spin initialiser to " << spath <<std::endl; 
    simulate.save_spins(spath);

    return 0;
}
