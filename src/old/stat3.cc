/* Executable to estimate several thermodynamic observables and their errors
 * for the classical QSI simulations by Gaussian and unbiased estimation.
 *
 * Usage: stat3 index size N Ns [L] <infile >>outfile
 **  index: of probed temperature, T = 1/(2.5 - i/36)
 **  size:  of the system (linear)
 **  Ns:    number of lines in file
 **  N:     number of unbiased samples
 **  L:     size of unbiased samples (defaults to half of available data)
 *
 * Infile assumed to have 3 columns of data: n=#vison, n2=#dipole, e=energy
 * Define n1 = n-2*n2, an estimate of #free vison
 * Prints a single line containing estimates and errors of
 **  index, temperature;
 **  n, muE=de/dn, muArrh, muArrh/muE = var(n)/n;
 **  same for n2; same for n1;
 **  e; C=de/dT
 *
 * Created on 10/10/2018
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
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cmath>

using namespace std;

double T, E0;

inline double sqr(double x) {return x*x;}

// Indexing double[4] of raw data and double[11] of means and (co)variances
enum :char {Nt = 0, N2 = 1, N1 = 2, E = 3};
inline char mean(char i) {return i;}
inline char var(char i) {return i+4;}
inline char cov(char i) {return i+8;}

// Reads a line of infile into a double[4]; contents will be n, n2, n1, e
void get_raw (double* store) {
    scanf("%lf %lf %lf", store, store+1, store+3);
    store[2] = store[0]-2*store[1];
    store[3] += E0;
}

// Swaps two double[4]'s
void swap (double* a, double* b) {
    double s[4];
    memcpy(s, a, 4*sizeof(double));
    memcpy(a, b, 4*sizeof(double));
    memcpy(b, s, 4*sizeof(double));
}

// Generates means and (co)variances of input raw data into double[11]
void sampling (const double* raw, size_t num, double* stat) {
    // Zeroing
    for (int i = 0; i < 11; ++i) stat[i]=0.0;
    // Means
    for (size_t i = 0; i < num; ++i)
        for (int j = 0; j < 4; ++j) 
            stat[mean(j)] += raw[4*i+j];
    for (int i = 0; i < 4; ++i) stat[i] /= num;
    // (Co)variances
    for (size_t i = 0; i < num; ++i) {
        for (int j = 0; j < 4; ++j)
            stat[var(j)] += sqr(raw[4*i+j]-stat[j]);
        for (int j = 0; j < 3; ++j)
            stat[cov(j)] += (raw[4*i+j]-stat[j]) * (raw[4*i+3]-stat[3]);
    }
    for (int i = 4; i < 11; ++i) stat[i] /= num;
}

/* Mean and stdev of some observable deducible from means & (co)variances
 * for all samples
 * f should take a double[11] and an optional integer parameter and
 **   give an estimate for the observable */
void mean_std (const double* samples, size_t num,
               double(*f)(const double*, int),
               int param, double* mean, double* std) {
    *mean = 0;
    for (size_t i = 0; i < num; ++i) 
        *mean += f(samples + 11*i, param);
    *mean /= num;

    *std = 0;
    for (size_t i = 0; i < num; ++i)
        *std += sqr(f(samples + 11*i, param) - *mean);
    *std = sqrt(*std / num);
}

// Returns the mean
inline double avg(const double* sample, int which) {return sample[which];}

// Calculates mu_E for a given type of N
inline double muE(const double* sample, int which) {
    return sample[cov(which)]/sample[var(which)];
}

// Calculates mu_Arrh for a given type of N
inline double muA(const double* sample, int which){
    return sample[cov(which)]/sample[mean(which)];
}

// Calculates the ratio muA/muE for a given type of N
inline double muR(const double* sample, int which){
    return sample[var(which)]/sample[mean(which)];
}

// Calculates the heat capacity
inline double C(const double* sample, int dum) {
    return sample[var(E)]/sqr(T);
}



int main(int argc, char** argv) {
    srand(time(NULL));

    int index;
    switch (argv[1][0]) {
    default:
        index = atoi(argv[1]);
        T = 1.0/(2.5 - index/36.0);
        break;
    case 'w':
        index = atoi(argv[1]+1);
        T = 1.0+0.05*index;
        break;
    case 'h':
        index = atoi(argv[1]+1);
        T = 2.0+0.1*index;
    }
    int size = atoi(argv[2]);
    E0 = 16.0*size*size*size;
    size_t N_raw = atoi(argv[3]);
    size_t N_sample = atoi(argv[4]);
    size_t L = N_raw/2;
    if (argc > 5)
        L = atoi(argv[5]);

    double* raw = new double[N_raw*4];
    for (size_t i = 0; i < N_raw; ++i)
        get_raw(raw + 4*i);

    double* sample = new double[N_sample*11];
    for (size_t i = 0; i < N_sample; ++i) {
        /* Select L random samples, scramble the dataset, record averages */
        for (size_t j = N_raw - 1; j >= N_raw - L; --j) {
            size_t chosen = rand() % j;
            swap(raw + 4*chosen, raw + 4*j);
        }
        sampling(raw + 4*(N_raw - L), L, sample + 11*i);
    }

    // PRINTING
    printf("%2d %10.8f ",index,T);
    
    double mn, std;
    // Various observables for N's
    for (int i = 0; i < 3; ++i) {
        mean_std(sample, N_sample, avg, i, &mn, &std);
        printf("%.8e %.8e ", mn, std);
        mean_std(sample, N_sample, muE, i, &mn, &std);
        printf("%11.8f %11.8f ", mn, std);
        mean_std(sample, N_sample, muA, i, &mn, &std);
        printf("%11.8f %11.8f ", mn, std);
        mean_std(sample, N_sample, muR, i, &mn, &std);
        printf("%11.8f %11.8f ", mn, std);
    }
    mean_std(sample, N_sample, avg, E, &mn, &std);
    printf("%.8e %.8e ", mn, std);
    mean_std(sample, N_sample, C, -1, &mn, &std);
    printf("%.8e %.8e\n", mn, std);

    delete[] raw;
    delete[] sample;

    return 0;
}
