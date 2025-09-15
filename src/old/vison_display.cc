/* Writes out a vison configuration from an output of autocorr
 *
 * Usage: vison_display infile N index
 ** infile: binary file of vison numbers stored as signed char
 ** N:      system size
 ** index:  of the stored vison config to list
 *
 * output is a list of sites with a vison on it, one per line, of format
 ** x y z charge
 *
 * Created on 16/12/2018
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

#include <vec3.hh>
#include <direction.hh>
#include <cstdio>
#include <cstdlib>

using namespace std;

// System size, needed for indexing routine
unsigned int N;

/* Position of dual diamond site #i
 * Order of sites in listing: x,y,z,fcc,diamond_sl */
vec3_int pos(size_t i) {
    unsigned dia = i%2; i/=2;
    unsigned fcc = i%4; i/=4;
    unsigned z = i%N;   i/=N;
    unsigned y = i%N;   i/=N;
    unsigned x = i%N;

    using namespace direction;
    return 8*vec3_int(x,y,z) + fcc_Ti[fcc] + diamond[dia];
}

int main(int argc, char** argv) {
    N  = atoi(argv[2]);
    size_t index = atoi(argv[3]);
    size_t sites = 8ul*N*N*N;

    FILE* f = fopen(argv[1], "rb");
    fseek(f, index*sites, SEEK_SET);
    char* raw = new char[sites];
    fread(raw, 1, sites, f);
    fclose(f);

    for (size_t i = 0; i < sites; ++i) 
        if (raw[i] != 0) {
            vec3_int r = pos(i);
            printf("%4d %4d %4d %2d\n",r[0],r[1],r[2],raw[i]);
        }
    
    delete[] raw;
    return 0;
}
