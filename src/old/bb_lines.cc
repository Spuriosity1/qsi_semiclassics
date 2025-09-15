/* Executable to list bb correlator along pinch point lines
 *
 * Usage: bb_lines N infile outfile [block]
 **   block specifies which block of infile is to be used (default: 0)
 *
 * Structure of outfile: each cut a separate block, each line is
 **  qx qy qz dq bb qq
 **  q*: wave vector in units of 2pi/Na
 **  dq: distance to pinch point in units of 1/a
 **  bb: the correlator
 **  qq: approx. vison charge correlator
 * Lines listed:
 **  [00k] k =  0..4
 **  [hh2] h = -1..1
 **  [hhh] h =  0..2
 *
 * Created on 11/08/2018
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

#include <magnetic_manager.hh>
#include <cstdio>
#include <cstdlib>

using namespace std;

#define TWOPI 6.283185307179586476925286766559005768394338798750211641949
#define SQRT2 1.414213562373095048801688724209698078569671875376948073176
#define SQRT3 1.732050807568877293527446341505872366942805253810380628055

void print_line(FILE* f, const magnetic_manager& c, int h, int k, double d) {
    fprintf(f,"%4d %4d %4d %10.6f %.6e %.6e\n",h,h,k,d,c.BB(h,k),c.vison(h,k));
}

int main(int argc, char** argv) {
    unsigned block = 0u;
    if (argc > 4)
        block = atoi(argv[4]);
    
    int n = atoi(argv[1]);
    magnetic_manager corr(n);

    FILE* in = fopen(argv[2], "rb");
    corr.load(in, block);
    fclose(in);

    FILE* out = fopen(argv[3], "w");
    for (int i = 0; i <= 4*n; ++i)
        print_line(out,corr,0,i,TWOPI*(i-2*n)/n);
    fprintf(out,"\n");
    for (int i = -n; i <= n; ++i)
        print_line(out,corr,i,2*n,TWOPI*SQRT2*i/n);
    fprintf(out,"\n");
    for (int i = 0; i <= 2*n; ++i)
        print_line(out,corr,i,i,TWOPI*SQRT3*(i-n)/n);
    fclose(out);

    return 0;
}
