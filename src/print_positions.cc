#include <fourier.hh>
#include <qsi.hh>
#include <string>
#include <spin.hh>
#include <plaq.hh>
#include <iostream>
#include <random>
#include <complex>
// USAGE: print_positions nx ny nz (spins | plaquettes)
#define EPS 1e-8




vec3 rand_unit(std::default_random_engine g){
    std::normal_distribution<double> n(0,1);
    vec3 retval;
    for (int i=0; i<3; i++){
        retval[i] = n(g);
    }
    
    return retval.normalise();
}


int main(int argc, const char** argv) {
    

    int n;
    if (argc < 3){
        std::cerr << "USAGE: print_positions n OUTFILE [randomise]" << std::endl;
        return 1;
    }
    n = atoi(argv[1]);
    fourier f(n);
    qsi qq(n, &f);
    
    std::string spinfile(argv[2]);
    spinfile += ".spins";
    std::string plaqfile(argv[2]);
    plaqfile += ".coupling";


    if (argc >= 4){
        std::default_random_engine g(atoi(argv[2]));
    

        std::normal_distribution<double> nd(0,1);

        for (int i=0; i<qq.n_spin();i++){
            // HORRENDOUS hack to access protected informaiton (and fill it with garbage)
            spin* s = (spin*)((void*)(qq.spin_no(i)));
            plaq* p = (plaq*)((void*)(qq.plaq_no(i)));
            s->heis() = rand_unit(g);
            
            p->g_constant = std::complex<double>(nd(g), nd(g));
            p->RK_potential = nd(g);

        }

    } else {
        qq.set_uniform_g(1.0);
        qq.set_uniform_RK(0);
    }

    qq.save_spins(spinfile);
    qq.save_couplings(plaqfile);

    qsi qq2(n, &f);

    qq2.load_spins(spinfile, true);
    qq2.load_couplings(plaqfile, true);

    for (int i=0; i<qq.n_spin();i++){
        vec3 diff = qq.spin_no(i)->heis() - qq2.spin_no(i)->heis();
        std::complex z = qq.plaq_no(i)->g_constant - qq2.plaq_no(i)->g_constant;
        double zz = qq.plaq_no(i)->RK_potential - qq2.plaq_no(i)->RK_potential;
        if ( diff.len() > EPS || std::abs(z) > EPS || std::abs(zz) > EPS){
            std::cout<<diff[0]<<" "<<diff[1]<<" "<<diff[2]<<" "<<z<<" "<<zz <<std::endl;
            // throw std::runtime_error("Broken IO code");
        }
    }
}