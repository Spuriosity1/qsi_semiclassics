/* Generates an analytic photon dispersion curve along the path of
 * dynamic_correlator::path4
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

#include <util.hh>
#include <cstdio>

int main(void) {
    for (int i = 0; i <= 700; ++i) {
        vec3 k = path4(0.02*i);
        printf("%5.2f %5.3f %5.3f %5.3f %10.8f\n",
               0.02*i,
               k[0], k[1], k[2],
               freq_analytic(k));
    }

    return 0;
}
