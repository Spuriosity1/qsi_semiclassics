/* Miscellaneous mathematical functions.
 *
 * Created on 12/12/2017
 * Copyright (C) 2017, 2018 Attila Szabó <as2372@cam.ac.uk>
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

#ifndef misc_hh
#define misc_hh

#include "vec3.hh"

// Checks if all components of an integer vector is divisble by some number
inline bool alldiv(const vec3_int v, int n) {
    return (v[0]%n==0) && (v[1]%n==0) && (v[2]%n==0); 
}

// Proper modular division for potentially negative inputs
inline int mod(int a, int b) {
    if (b<0) b=-b;
    int m = a%b;
    if (m<0) m+=b;
    return m;
}

// Modular division in the range [-b/2, b/2)
inline int rem(int a, int b) {
    if (b<0) b=-b;
    int r = mod(a,b);
    if (r < (b+1)/2)
	return r;
    else
	return r-b;
}

// TODO move to basic_parser.hh
inline void my_assert(bool statement, const std::string& msg){
    if (!statement){
        throw std::runtime_error(msg);
    }
}

inline vec3_int get_system_size(const std::string& dfname){
    // this is garbage code
    std::string s = dfname;
    std::string delimiter = "%";
    std::string token;
    unsigned pos;
    // add a delimiter at the end for fun
    s += delimiter;

    int L  = -1;
    int Lx = -1;
    int Ly = -1;
    int Lz = -1;

    while ((pos = s.find(delimiter)) != std::string::npos && s.length() > 0) {
        token = s.substr(0, pos);
        size_t pos_equals = token.find("=");
        if (pos_equals != std::string::npos){
            if (token.substr(0,pos_equals)=="L"){
                L = std::stoi(token.substr(pos_equals+1));
            } else if (token.substr(0,pos_equals)=="Lx") {
                Lx = std::stoi(token.substr(pos_equals+1));
            } else if (token.substr(0,pos_equals)=="Ly") {
                Ly = std::stoi(token.substr(pos_equals+1));
            } else if (token.substr(0,pos_equals)=="Lz") {
                Lz = std::stoi(token.substr(pos_equals+1));
            }
        }
        s.erase(0, pos + delimiter.length());
    }

    if (L == -1) {
        my_assert( Lx > 0 && Ly > 0 && Lz > 0, "If L is unspecified, all three of Lx, Ly, Lz must be specified." );
    } else {
        my_assert( Lx ==-1 && Ly ==-1 && Lz ==-1, "If L is specified, none of Lx, Ly, Lz may be specified." );
        Lx = L;
        Ly = L;
        Lz = L;
    }

    return vec3_int(Lx,Ly,Lz);
}

#endif
