
#include <basic_parser.hh>
#include <qsi.hh>
#include <ptetra.hh>
#include <plaq.hh>
#include <spin.hh>
#include <tetra.hh>
#include <string>
#include <iostream>

typedef basic_parser<int, unsigned, double> parser_t;


/**
 * @brief Generates a ground state by hand based on four fluxes
 * 
 * @param argc 
 * @param argv 
 * @return int 
 */


int main(int argc, const char* argv[]){
    parser_t p;
    unsigned Lx, Ly, Lz;
    // unsigned L;

    double B[4]; // magnetic field (units of pi)
    std::string ofile;
    double T_hot, T_cold;
    unsigned burnin, n_anneal, seed, sweep;
    // float gp0, gp1, gp2, gp3 RK;


    p.declare("Lx", &Lx);
    p.declare("Ly", &Ly);
    p.declare("Lz", &Lz);

    // p.declare("L", &L);
    
    p.declare_optional("b0", B,0.);
    p.declare_optional("b1", B+1,0.);
    p.declare_optional("b2", B+2,0.);
    p.declare_optional("b3", B+3,0.);

    p.declare("name", &ofile);

    p.declare("temperature",&T_cold);
    p.declare("n_sweep", &sweep);
    p.declare("n_burnin", &burnin);
    p.declare("seed", &seed);
    p.declare("T_hot", &T_hot);
    p.declare("n_anneal", &n_anneal);




    if (argc < 3){
        std::cerr << "Usage: generate_initialiser INDIR OUTDIR [options...]" << std::endl;
        return 1;
    }
    p.from_file(argv[1]);
    p.from_argv(argc, argv, 3);
    p.assert_initialised();

    // Check that the fluxes are achievable
    if (std::abs(std::remainder( (B[0] + B[1] + B[2] + B[3]), 2 )) > 1e-6  ) {
        std::cerr << "FATAL: specified fluxes\n";
        std::cerr << B[0] <<"pi "<<B[1]<<"pi "<<B[2]<<"pi "<<B[3]<<"pi \n";
        std::cerr << "Do not sum to a multiple of 2pi\n";
        return 1;
    }

    std::filesystem::path outdir(argv[2]);

    qsi simulate(Lx, Ly, Lz, NULL, NULL, NULL,seed);

    for (int j=0; j<simulate.n_spin(); j++){
        // evil casting away const
        // delightfully devilish
        plaq* p = (plaq*) simulate.plaq_no(j);
        p->g_constant = std::polar(1., -B[p->sublat()]*M_PI);
    }
    
     // sample everywhere
    mc_sampler all_samples(mc_sampler::sequential, simulate.n_spin());

    // Burn-in
    simulate.MC(T_hot, burnin, all_samples);

    double T = T_hot;
    double factor = 1;

    if (n_anneal > 0) factor = exp((log(T_hot)-log(T_cold)) / n_anneal);
    fprintf(stderr, "#T    ");
    double errs[4] = {0,0,0,0};
    // Anneal from high to low temperature
    for (unsigned i = 0; i < n_anneal; ++i) {
        T /= factor;
        simulate.MC(T, sweep, all_samples, mc_sz_mode::finite_t, mc_angle_mode::finite_t_gradient);
        // check for success
        for (int mu=0; mu<4; mu++) {errs[mu] = 0;}

        for (unsigned J=0; J<simulate.n_spin(); J++){
            auto p = simulate.plaq_no(J);
            int nu = p->sublat();
            errs[nu] += abs(p->ring() - std::polar(1., -B[nu] * M_PI) );
        }

        // Print the vison order parameter (redirect this to a file in a bash script if needed)
        fprintf(stderr, "%+9.6f  %+.6f %+.6f %+.6f %+.6f %+9.6e\n", T, errs[0], errs[1], errs[2], errs[3], simulate.energy()/simulate.n_spin());
    }

    // check that the process succeeded
    // check deviation from commanded fluxes
    for (int mu=0; mu<4; mu++) {
        errs[mu] *= 4. / simulate.n_spin();
    }

    fprintf(stderr, "Average flux errors:\nb0            b1             b2            b3\n");
    fprintf(stderr, "%+6.6e %+6.6e %+6.6e %+6.6e\n", errs[0], errs[1], errs[2], errs[3]);

    // look for outliers and other statistics
    double largest_error = 0;
    const plaq* erroneous_plaq = NULL;
    double sum_e = 0;
    double sum_e2 = 0;
    for (unsigned i=0; i<simulate.n_spin(); i++){
        const plaq* p = simulate.plaq_no(i);
        double e = abs(p->ring() - std::polar(1., -B[p->sublat()] * M_PI) );
        if (e > largest_error){
            largest_error = e;
            erroneous_plaq = p;
        }
        sum_e += e;
        sum_e2 += e*e;
    }
    fprintf(stderr, "Worst plaquette: sl %1d, error %f, actual flux %f pi\n", 
        erroneous_plaq->sublat(), largest_error, erroneous_plaq->B()/M_PI);

    sum_e /= simulate.n_spin();
    sum_e2 /= simulate.n_spin();
    
    fprintf(stderr, "Mean error %f, stdev %f\n", 
        sum_e, sqrt(sum_e2 - sum_e*sum_e));
    
    
    std::cerr << "Saving fluxes...\n";
    
    const auto& bpath = outdir / (ofile + "_B.csv");
    const auto& spath = outdir / (ofile + ".spins");
    const auto& gpath = outdir / (ofile + ".gauge");
    const auto& gpath2 = outdir / (ofile + ".gauge2");

    std::cout << std::string(bpath) << "\n";
    std::cout << std::string(spath) << "\n";
    std::cout << std::string(gpath) << "\n";
    std::cout << std::string(gpath2) << "\n";

    simulate.save_B(bpath);
    simulate.save_spins(spath);
    simulate.save_A(gpath);
    
    FILE* fp = fopen(gpath2.c_str(), "w");
    fprintf(fp, "up_tetra_idx pyro_sl A\n");
    for (int i=0; i< simulate.n_spin()/4; i++){
        const tetra* t = simulate.tetra_no(2*i);
        for (int mu=0; mu<4; mu++){
            fprintf(fp, "%3d %1d %9.6f\n", i, mu,std::arg(t->neighbour(mu)->xy()));
        }
    }

    return 0;
}
