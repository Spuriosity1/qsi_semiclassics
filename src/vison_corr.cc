/* Executable to generate full real-space vison correlators
 *
 * Usage: vison_corr N T output seed [samples [sweep [burn_in]]] > outfile
 ** N:           system size, required
 ** T:           temperature, required
 ** output:      output filename, required
 ** seed:        index of seed in the prestored array, required
 ** samples:     number of MC samples, default: 16384
 ** sweep:       number of MC sweeps between samples, default: 4
 ** burn_in:     number of MC sweeps for burning in, default: 1024
 *
 * output will contain the real space vison correlators as a 3d binary array
 **      of double precision floats.
 * outfile is a 3-column table, each line gives, for one sample:
 **      #visons,  #nn_dipoles,  energy
 *
 * Created on 08/10/2018
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
#include <fourier.hh>
#include <ptetra.hh>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <basic_parser.hh>

using namespace std;

const static string VERSION = "1.1";

int main (int argc, char** argv) {
    // ---------- Reading parameters ----------
    if (argc < 5)
        throw std::runtime_error("Required parameters missing");
    
    int n ; // system size
    double T ; //temperature
    unsigned seed ; // seed index
    
    // Optional parameters
    unsigned
        sample = 16384,
        sweep  = 4,
        burnin = 1024;

    basic_parser<int, unsigned, double> p("Static Structure Factors (qsi_static) v" + VERSION + currentDateTime());
    p.declare("system_size",&n);
    p.declare("temperature",&T);
    p.declare("seed", &seed);
    p.declare("sample", &sample);
    p.declare("n_sweep", &sweep);
    p.declare("n_burnin", &burnin);


    double gphase=0;
    p.declare("phi_g", &gphase);

    p.declare("T_hot", &T);
    // p.declare("n_anneal", &n_anneal);

    //read the params
    p.from_file(argv[1]);

    std::filesystem::path output_path(argv[2]);
    std::filesystem::path in(argv[1]);
    std::string prefix(in.filename().replace_extension());
    prefix += p.from_argv(argc, argv, 3);
    p.assert_initialised();


    fourier f(2*n); // since we need a point every a_0/4
    qsi simulate(n,&f,NULL,NULL,seed);



    simulate.set_uniform_g(polar(1.,gphase));

    size_t N = 64ul*n*n*n;
    double* corr = new double[N];
    memset(corr, 0, N*sizeof(double));
    complex<double>* fft = &f(0,0,0);
    
    simulate.MC(T,burnin);
    
    // Simulate & record data
    for (unsigned i = 0; i < sample; ++i) {
        simulate.MC(T,sweep);
        f.clear_input();
        for (size_t v = 0; v < simulate.n_spin()/2; ++v) 
            f( (simulate.ptetra_no(v)->pos()) /2 ) =
                simulate.ptetra_no(v)->vison() * 1.0;
        f.transform();
        for (size_t v = 0; v < N; ++v)
            corr[v] += norm(fft[v]);
    }

    // Load recorded correlators into fourier and transform them "back"
    for (size_t v = 0; v < N; ++v)
        fft[v] = corr[v];
    f.transform();
    for (size_t v = 0; v < N; ++v)
        corr[v] = fft[v].real();

    // Save data
    FILE* out = fopen(argv[3],"wb");
    fwrite(corr, sizeof(double), N, out);

    return 0;
}
