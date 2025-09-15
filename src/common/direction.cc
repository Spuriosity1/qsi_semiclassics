/* Implementation of directions.hh
 *
 * Created on 09/06/2018
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

#include "direction.hh"

const vec3_int direction::pyro[4] = {
    vec3_int( 1, 1, 1),
    vec3_int( 1,-1,-1),
    vec3_int(-1, 1,-1),
    vec3_int(-1,-1, 1)};

const vec3 direction::r[4] = {pyro[0]*0.125, pyro[1]*0.125,
                              pyro[2]*0.125, pyro[3]*0.125};

const vec3_int direction::diamond[2] = {vec3_int(), vec3_int(2,2,2)};

const vec3_int direction::fcc_Dy[4] = {
    vec3_int(0,0,0),
    vec3_int(0,4,4),
    vec3_int(4,0,4),
    vec3_int(4,4,0)};

const vec3_int direction::fcc_Ti[4] = {
    vec3_int(4,4,4),
    vec3_int(4,0,0),
    vec3_int(0,4,0),
    vec3_int(0,0,4)};

// Vectors from a plaquette centre to its six vertices
// first index: plaquette sublattice
// second index: enumeration
const vec3_int direction::plaqt[4][6] = {
    {
	vec3_int( 0,-2, 2),
	vec3_int( 2,-2, 0),
	vec3_int( 2, 0,-2),
	vec3_int( 0, 2,-2),
	vec3_int(-2, 2, 0),
	vec3_int(-2, 0, 2)},
    {
	vec3_int( 0, 2,-2),
	vec3_int( 2, 2, 0),
	vec3_int( 2, 0, 2),
	vec3_int( 0,-2, 2),
	vec3_int(-2,-2, 0),
	vec3_int(-2, 0,-2)},
    {
	vec3_int( 0,-2,-2),
	vec3_int(-2,-2, 0),
	vec3_int(-2, 0, 2),
	vec3_int( 0, 2, 2),
	vec3_int( 2, 2, 0),
	vec3_int( 2, 0,-2)},
    {
	vec3_int( 0, 2, 2),
	vec3_int(-2, 2, 0),
	vec3_int(-2, 0,-2),
	vec3_int( 0,-2,-2),
	vec3_int( 2,-2, 0),
	vec3_int( 2, 0, 2)}
};

// Square roots of 2, 3, and 6 for normalisation
#define S2 1.414213562373095048801688724209698078569671875376948073176
#define S3 1.732050807568877293527446341505872366942805253810380628055
#define S6 2.449489742783178098197284074705891391965947480656670128432

const vec3 direction::axis[4][3] = {
    {vec3( 1, 1,-2)/S6, vec3(-1, 1, 0)/S2, vec3( 1, 1, 1)/S3},
    {vec3( 1,-1, 2)/S6, vec3(-1,-1, 0)/S2, vec3( 1,-1,-1)/S3},
    {vec3(-1, 1, 2)/S6, vec3( 1, 1, 0)/S2, vec3(-1, 1,-1)/S3},
    {vec3(-1,-1,-2)/S6, vec3( 1,-1, 0)/S2, vec3(-1,-1, 1)/S3}
};

const vec3_int direction::cubic[3] = {
	vec3_int(8,0,0), vec3_int(0,8,0), vec3_int(0,0,8)
};
