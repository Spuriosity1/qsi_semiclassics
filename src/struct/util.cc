/* Implementation of functions in util.hh
 *
 * Created on 07/08/2018
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
#include "util.hh"
#include <cmath>

using namespace std;

#define PI4   0.785398163397448309615660845819875721049292349843776455243
#define TWOPI 6.283185307179586476925286766559005768394338798750211641949

inline double sqr (double x) {return x*x;}

/* The formula given in PRB 95, 134439, after unwrapping the layers of
 * unwieldy formalism, is
 **    omega = sqrt(12) g * sqrt[sum_ij sin^2(k.Delta_ij)]
 * where Delta_ij is the vector from a plaquette centre of sublattice i
 * to a spin of sublattice j.
 * Delta_ij are listed in the array plaqt in this code in units of a/8. */
double freq_analytic (const vec3& k) {
    // The sum under the square root
    double xi = 0.0;
    for (int i = 0; i < 4; ++i)
        /* We don't worry much which sublattice j actually belongs to
         * plaqt[i][j] (j<3) should traverse all non-i sublattices */
        for (int j = 0; j < 3; ++j)
            // [k] = 2pi/a, [plaqt]=a/8, hence the pi/4
            xi += sqr(sin(PI4 * k % direction::plaqt[i][j]));

    return sqrt(12*xi);
}

/* General purpose interpolator
 * Looks up x in xp and returns a linear intrepolation based on yp
 * The length of both xp and yp are given by n
 * xp is assumed to be sorted ascending */
double interp(double x, const double* xp, const double* yp, unsigned n) {
    if (x < xp[0])
        throw std::runtime_error("x is smaller than first element of xp");
    // x is in the right range, check segments one by one
    for (unsigned i = 0; i < n-1; ++i)
        // if it so far, x > xp[i]
        if (x <= xp[i+1])
            // x belongs to this segment
            return yp[i] + (x-xp[i]) * (yp[i+1]-yp[i]) / (xp[i+1]-xp[i]);
    // if it got so far, x > xp[n-1], which is no good
    throw std::runtime_error("x is greater than last element of xp");
}

// Path of dynamic_correlator::path4 as an array (units: 2pi/a)
const double p4[4][8] = {
    {    0,    4,    6,    7,   10,   12,   13,   14},
    {  0.0,  1.0,  1.0, 0.75,  0.0,  0.5, 1.00,  1.0},
    {  0.0,  0.0,  0.5, 0.75,  0.0,  0.5, 0.25,  0.0},
    {  0.0,  0.0,  0.0, 0.00,  0.0,  0.5, 0.25,  0.0}};

// Path of dynamic_correlator::path2 as an array (units: 2pi/a)
const double p2[4][6] = {
    {   0,   2,   3,   4,   5,   7},
    { 0.0, 1.0, 1.0, 0.5, 0.0, 1.0},
    { 0.0, 0.0, 0.5, 0.5, 0.0, 1.0},
    { 0.0, 0.0, 0.0, 0.5, 0.0, 0.0}};

vec3 path4 (double l) {
    return vec3(interp(l,p4[0],p4[1],8),
                interp(l,p4[0],p4[2],8),
                interp(l,p4[0],p4[3],8));
}

vec3 path2 (double l) {
    return vec3(interp(l,p2[0],p2[1],6),
                interp(l,p2[0],p2[2],6),
                interp(l,p2[0],p2[3],6));
}

/**
 * @brief A coarse approximation to the Coulomb connection integrated along the link connecting the (real) tetrahedron t0 to t1.
 * 
 * @param plaq_loc Origin of the plaquette
 * @param n a normal vector connecting the (-) monopole to the (+) monopole
 * @param t0 tetrahedron 0 location
 * @param t1 tetrahedron 1 location
 * @return double 
 */
double Dirac_dipole(const vec3& plaq_loc, const vec3& normal_vector, const vec3& t0, const vec3& t1){
    vec3 n = normal_vector;
    n.normalise();
    vec3 r = 0.5*(t0 + t1) - plaq_loc;

    const double a = normal_vector.len() / 2;
    
    double z = r % n;
    double p = sqrt(r.len2() - z*z);
    if(p < 0.1 || std::isnan(p)){
        throw std::logic_error("Impossible plaquette distance");
    }
    // cross product
    vec3 e_phi = (n * r).normalise();

    double retval = e_phi % (t1 - t0) * 0.5* ( 
          (z - a)/sqrt((z - a)*(z - a) + p*p) 
        - (z + a)/sqrt((z + a)*(z + a) + p*p)
    ) / p ;
    
    if (std::isnan(retval)){
        throw std::logic_error("nan detected");
    }
    return retval;
}



double t3_metric(const vec3_int& x, const vec3_int& y, const vec3_int& L){
    vec3_int D(x);
    vec3_int Y(y);
    D.wrap(L);
    D -= Y.wrap(L);

    for (int i=0; i<3; i++)
        D[i] = abs(D[i]) < L[i]/2 ? D[i] : L[i] - abs(D[i]);

    return D.len();
}

int t3_metric2(const vec3_int& x, const vec3_int& y, const vec3_int& L){
    vec3_int D(x);
    vec3_int Y(y);
    D.wrap(L);
    D -= Y.wrap(L);

    for (int i=0; i<3; i++)
        D[i] = abs(D[i]) < L[i]/2 ? D[i] : L[i] - abs(D[i]);

    return D.len2();
}