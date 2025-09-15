/* Executable to average over several files containing real numbers
 *
 * Assumes all files are binary containing a stream of double precision floats
 * and of the same length.
 * These are loaded, each position is averaged over, and then printed into
 * the output file.
 *
 * Usage: average infile1 [infile2 [infile3 [...]]] outfile
 *
 * Created on 10/08/2018
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
#include <stdexcept>

using namespace std;

int main (int argc, const char** argv) {
    if (argc < 3)
	throw std::runtime_error("No input files");

    // First file
    FILE* f = fopen(argv[1],"rb");
    fseek(f, 0L, SEEK_END);
    size_t N = ftell(f) / sizeof(double);
    rewind(f);

    double* sum  = new double[N]; // array for summing up the final result
    double* read = NULL; // array for reading each extra file
    if (argc > 3)
	read = new double[N]; // only needed if there *are* extra files

    fread(sum, sizeof(double), N, f);
    fclose(f);

    // Further files
    for (int i = 2; i < argc-1; ++i) {
	f = fopen(argv[i],"rb");
	fread(read, sizeof(double), N, f);
	fclose(f);

	for (size_t k = 0; k < N; ++k)
	    sum[k] += read[k];
    }

    // turn sum into an average
    for (size_t k = 0; k < N; ++k)
	sum[k] /= (argc-2);

    f = fopen(argv[argc-1],"wb");
    fwrite(sum, sizeof(double), N, f);
    fclose(f);

    return 0;
}
