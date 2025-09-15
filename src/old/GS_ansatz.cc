/* Curl-free ansatz ground state energy of the semiclassical QSI model with two
 * visons in it
 *
 * Usage: GS_ansatz x y z N --gphase g_phase
 ** x,y,z: coordinates of one vison (has to be a valid diamond lattice site,
 **                                  the other one is in the origin)
 ** N:     system size in which to evaluate
 *
 * Outputs three numbers:
 **  the ground state energy,
 **  the same in quadratic approximation,
 **  the estimate 2mu-alpha/r, where mu is obtained from fitting
 *
 * Created on 11/09/2018
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

#include <fourier.hh>
#include <direction.hh>
#include <misc.hh>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cmath>

using namespace std;

inline complex<double> gamma(const vec3& k) {
    return complex<double>(+cos(k[0]/4) * cos(k[1]/4) * cos(k[2]/4),
                           -sin(k[0]/4) * sin(k[1]/4) * sin(k[2]/4));
}
inline complex<double> cis(double x) {return polar(1.0, x);}

inline complex<double> cis(const vec3& k, unsigned i) {
    return polar(1.0, k % direction::pyro[i] * 0.25);
}

inline double W1(double B) {return 1-cos(B);}
inline double W2(double B) {return 0.5*B*B;}

#define TWOPI 6.283185307179586476925286766559005768394338798750211641949
#define PI2 1.570796326794896619231321691639751442098584699687552910487
#define MU 7.872367608

int main(int argc, char** argv) {
    if ( argc < 5 ){
        fprintf(stderr, "Insuffieintly many arguments supplied\n");
    }

    const string ps_gphase= "--gphase";

    complex<double> g_constant = 1.;

    if ( argc >= 6 ){
        // Ignore extra arguments silently
        if ( argc == 6 || argv[5] != ps_gphase ){
            fprintf(stderr, "Unrecognised optional argument specifier %s \n", argv[5]);
            throw std::runtime_error("Unrecognised argument");
        }

        g_constant = std::polar(1., atof(argv[6]));
        // g_constant = atof(argv[6]);
        fprintf(stderr, "Using g phase  = %f e^i%f\n", abs(g_constant), arg(g_constant));
    }

    int x = atoi(argv[1]);
    int y = atoi(argv[2]);
    int z = atoi(argv[3]);
    vec3_int r(x,y,z);
    size_t N = atoi(argv[4]);
    size_t NN = 8ul*N*N*N;

    if ( (mod(x,2)!=mod(y,2)) || (mod(x,2)!=mod(z,2)) ){
        fprintf(stderr, "Invalid diamond site: x,y,z must be congruent mod 2\n");
        throw std::runtime_error("Invalid diamond site");
    }

    bool sublat;
    switch(mod(x+y+z, 4)) {
    case 0:
        sublat = true;
        break;
    case 3:
        sublat = false;
        break;
    default:
        fprintf(stderr, "Invalid diamond site: x+y+z must be 0 mod 4\n");
        throw std::runtime_error("Invalid diamond site");
    }

    fourier f(N);
    double E1 = 0.0, E2 = 0.0;

    for (unsigned i = 0; i < 4; ++i) {
        for (unsigned x = 0; x < 2*N; ++x)
            for (unsigned y = 0; y < 2*N; ++y)
                for (unsigned z = 0; z < 2*N; ++z) {
                    vec3_int kZ(x,y,z);
                    vec3 k = TWOPI/N * kZ;
                    if (sublat)
                        f(kZ) = (1.0 - cis(k,i)*conj(gamma(k))) *
                            (1.0 - cis(k%r * -0.25));
                    else
                        f(kZ) = (1.0 - cis(k,i)*conj(gamma(k))) +
                            cis(k%r * -0.25) * (cis(k,i)-gamma(k));
                    f(kZ) *= PI2 / (1.0 - norm(gamma(k)));
                }
        f(0,0,0) = 0.0;
        f(N,N,N) = 0.0;

        f.transform();

        for (unsigned x = 0; x < 2*N; ++x)
            for (unsigned y = 0; y < 2*N; ++y)
                for (unsigned z = 0; z < 2*N; ++z) {
                    E1 += W1( (g_constant*f(x,y,z)).real() / NN);
                    E2 += W2( (g_constant*f(x,y,z)).real() / NN);
                }
    }

    printf("%.16f %.16f %.16f\n", E1, E2, 2*MU-4*M_PI/sqrt(x*x+y*y+z*z));

    return 0;
}
