
#include <basic_parser.hh>
#include <qsi.hh>
#include <ptetra.hh>
#include <plaq.hh>
#include <spin.hh>
#include <string>
#include <iostream>


typedef basic_parser<int, unsigned, double> parser_t;

int main(int argc, const char* argv[]){
    parser_t params;
    unsigned Lx, Ly, Lz;
    // unsigned L;

    double g[4];
    std::string ofile;
    // float gp0, gp1, gp2, gp3 RK;


    params.declare("Lx", &Lx);
    params.declare("Ly", &Ly);
    params.declare("Lz", &Lz);

    // params.declare("L", &L);
    
    params.declare_optional("g0", g,1.0);
    params.declare_optional("g1", g+1,1.0);
    params.declare_optional("g2", g+2,1.0);
    params.declare_optional("g3", g+3,1.0);

    params.declare("name", &ofile);


    if (argc < 3){
        std::cerr << "Usage: generate_couplings OUTDIR [options...]" << std::endl;
        return 1;
    }
    params.from_argv(argc, argv, 2);
    params.assert_initialised();

    std::filesystem::path outdir(argv[1]);

    qsi simulate(Lx, Ly, Lz, NULL, NULL, NULL);

    for (unsigned j=0; j<simulate.n_spin(); j++){
        // evil casting away const
        plaq* p = (plaq *) simulate.plaq_no(j);
        p->g_constant = std::complex<double>(g[p->sublat()], 0); 
    }

    // std::string pref = couplingfile_name+params.get_outfile_prefix(false, false);]
    
    
    std::cerr << "Saving couplings...\n";
    const auto& fpath = outdir / (ofile + ".couplings");
    simulate.save_couplings(fpath);
    std::cout << std::string(fpath) << "\n";

    return 0;
}
