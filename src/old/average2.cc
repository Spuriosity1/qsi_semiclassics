/* Executable to average over several files containing real numbers
 * 
 * Assumes all files are binary containing a stream of double precision floats.
 * These are loaded, each position is averaged over, and then printed into 
 * the output file. 
 * 
 * If any file has a different length to the original 
 * 
 * Visualise the following:
 * 
 * infile1 = [ x11 x12 x13 ... ] +
 * infile2 = [ x21 x22 x23 ... ] +
 * infile3 = [ x31 x32 x33 ... ]
 *            =
 * average = [ x#1 x#2 x#3 ... ]
 *
 * 
 * Usage: average infile1 [infile2 [infile3 [...]]] > outfile
 * Pipe this into `od -F` to get the float representation.
 *
 * Created on 07/04/2022
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
#include <cstdlib>
#include <sys/stat.h>

using namespace std;

long GetFileSize(const char* filename)
{
    struct stat stat_buf;
    int rc = stat(filename, &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}

int main (int argc, const char** argv) {
    if (argc < 2)
	throw std::runtime_error("No input files");

    // Get the size of the first argument
    long Nl = GetFileSize(argv[1]);
    if (Nl == -1) {
        fprintf(stderr, "File not found:%s",argv[1]);
        throw std::runtime_error("File Not Found");
    }
    // Explicit cast to avoid signed compaorison warnings
    size_t N = (size_t) Nl;

    // printf("%ld =?= %lu", Nl, N);

    double* sum = new double[N];
    double* buffer = new double[N];

    // Read the first file into memory
    FILE* f = fopen(argv[1],"rb");
    fread(sum, sizeof(double), N, f);
    fclose(f);


    for (size_t i=2; i<argc; i++){
        // check that the new file has the correct size
        long Ml = GetFileSize(argv[i]);
        if (Ml == -1) {
            fprintf(stderr, "File not found:%s\n",argv[1]);
            throw std::runtime_error("File Not Found");
        } else if (Ml != Nl) {
            fprintf(stderr, "File %s has length (%ld) inconsistent with file %s (%ld)\n",argv[i],Ml,argv[1],Nl);
            throw std::runtime_error("File Not Found");
        }
        FILE* f = fopen(argv[i], "rb");
        fread(buffer, sizeof(double), N, f);
        for (size_t j=0; j<N; j++){
            sum[j] += buffer[j];
        }

        fclose(f);
        
    }

    
    delete[] buffer;

    // Normalise
    for (size_t j=0; j<N; j++){
        sum[j] /= (double) (argc-1);
    }

    // output
    fwrite(sum, sizeof(double), N, stdout);
    
    delete[] sum;
    return 0;
}
