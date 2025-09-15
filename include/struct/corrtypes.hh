/* Bit flags for the different things whose correlators can be evaluated
 *
 * Created on 11/04/2019
 * Copyright (C) 2019 Attila Szabó <as2372@cam.ac.uk>
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

#ifndef corrtypes_hh
#define corrtypes_hh

namespace corrtypes {
    // <S^z_i(q) S^z_j(-q)>
    static const unsigned SZZ = 0x01u;
    static const unsigned SXX = 0x02u;
    static const unsigned SYY = 0x04u;
    // <S^+_i(q) S^-_j(-q)>
    static const unsigned SPM = 0x08u;
    // < O_ring > = <B(q)> = arg(S^+S^-...)
    static const unsigned ARG = 0x10u;
    // = Im(S^+S^-...)
    static const unsigned IMAG = 0x20u;
    // = sin(S^+S^-...)
    static const unsigned SIN = 0x40u;
    static const unsigned NOISE = 0x80u;
    static const unsigned COS = 0x100u;

    static const unsigned ALL = SZZ | SXX | SYY | SPM | ARG | IMAG | SIN | NOISE | COS;
}

#endif
