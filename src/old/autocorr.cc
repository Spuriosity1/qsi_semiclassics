/* Executable to obtain vison autocorrelators of MC protocol
 * This program does pure classical MC based on the ringflip Hamiltonian.
 *
 * Usage: autocorr N T output_file [--seed seed] [--sample n_samples] [--burnin burn_in] [--gphase g_phase]
 ** N:           system size, required
 ** T:           temperature, required
 ** output_file: name of output file, required
 ** seed:        index of seed in the prestored array, default: 0
 ** n_samples:     number of MC samples, default: 16384
 ** burn_in:     number of MC sweeps for burning in, default: 1024
 ** g_phase:     Complex 'twist' to add to g. default: 0
 *
 * output_file will contain the values of S^z for all non-burn-in samples
 ** in a binary file
 *
 * Created on 16/09/2018
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
#include <spin.hh>
#include <plaq.hh>
#include <ptetra.hh>
#include <fourier.hh>
#include <cstdio>
#include <cstdlib>
// #include <iostream>

using namespace std;

int main (int argc, char** argv) {
    // ---------- Reading parameters ----------
    if (argc < 4)
        throw std::runtime_error("Required parameters missing");

    int n = atoi(argv[1]); // system size
    double T = atof(argv[2]); //temperature


    // Optional Command line parameters
    
    const string ps_seed = "--seed";
    const string ps_sample = "--sample";
    const string ps_burnin = "--burnin";
    const string ps_gphase = "--gphase";

    int seed = 0;
    unsigned  sample = 16384;
    unsigned  burnin = 1024;
    plaq::g_constant = 1.;

    int argi = 4;
    
    while (argi < argc) {
        if (argv[argi] == ps_seed) {
            seed = atoi(argv[argi+1]);
            argi += 2;
        } else if(argv[argi] == ps_sample) {
            sample = atoi(argv[argi+1]);
            argi += 2;
        } else if(argv[argi] == ps_burnin) {
            burnin = atoi(argv[argi+1]);
            argi += 2;
        } else if(argv[argi] == ps_gphase){
            plaq::g_constant = std::polar(1., atof(argv[argi+1]));
            argi += 2;
        } else throw std::runtime_error("Invalid optional parameters");
    }




    fprintf(stderr, "Simulating with parameters\n");
    fprintf(stderr, "System Size   %d\n", n);
    fprintf(stderr, "Temperature   %e g\n", T);
    fprintf(stderr, "Seed          %d\n", seed);
    fprintf(stderr, "Sample        %d\n", sample);
    fprintf(stderr, "Burn in       %d\n", burnin);
    fprintf(stderr, "g phase       %f e^i%f\n", abs(plaq::g_constant), arg(plaq::g_constant));

    fourier f(n);
    qsi simulate(n,&f,NULL,NULL,seed);

    fprintf(stderr, "SIM SIZE %lu \n", simulate.n_spin());


    char* Sz = new char[simulate.n_spin()/2];

    // Burn-in
    fprintf(stderr, "burning in...\n");
    simulate.MC(T,burnin);

    FILE* out = fopen(argv[3],"wb");

    // Simulate & record data
    for (unsigned i = 0; i < sample; ++i) {
        fprintf(stderr, "[ sim ] %d / %d \n", i+1, sample);
        simulate.MC(T,1);
        for (size_t j = 0; j < simulate.n_spin()/2; ++j){
            Sz[j] = simulate.ptetra_no(j)->vison();
        }

        fwrite(Sz, sizeof(char), simulate.n_spin()/2, out);
    }

    fclose(out);
    delete[] Sz;

    return 0;
}
