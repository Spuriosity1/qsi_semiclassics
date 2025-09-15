/* Filters free visons from vison autocorrelation data
 *
 * Usage: autocorr_free N infile outfile
 ** N: system size
 *
 * infile & outfile of the same form as the vison incarnation of autocorr
 *
 * Created on 01/10/2018
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
#include <plaq.hh>
#include <ptetra.hh>
#include <fourier.hh>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

int main (int argc, char** argv) {
    try{
        // ---------- Reading parameters ----------
        int n = atoi(argv[1]); // system size
        size_t N = 8ul*n*n*n; // number of dual diamond sites

        fourier f(n);
        qsi simulate(n,&f,NULL,NULL,0);

        FILE* in = fopen(argv[2], "rb");
        FILE* out = fopen(argv[3], "wb");

        char* inc = new char[N];
        char* outc = new char[N];
        
        while (!feof(in)) {
            fread(inc, N*sizeof(char), 1, in);
            memcpy(outc, inc, N*sizeof(char));
            for (size_t i = 0; i < 2*N; ++i) {
                const plaq* p = simulate.plaq_no(i);
                outc[p->tA()->index] += inc[p->tB()->index];
                outc[p->tB()->index] += inc[p->tA()->index];
            }
            fwrite(outc, N*sizeof(char), 1, out);
        }
        
        fclose(in);
        fclose(out);
        
        delete[] inc;
        delete[] outc;
        
        return 0;
    } catch (const char* s) {
        fprintf(stderr,"%s",s);
    }
}
