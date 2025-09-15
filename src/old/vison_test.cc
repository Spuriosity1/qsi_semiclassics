/* Executable to test vison insertion routines
 *
 * Created on 22/08/2018
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

#include <qsi.hh>
#include <fourier.hh>
#include <cstdio>

using namespace std;

int main() {
    try {
        fourier f(4);
        qsi simulate(4,&f);

        simulate.add_vison_pair(vec3_int(4,4,4),
                           vec3_int(24,24,20));
        printf("%d %d %lu\n",
               simulate.vison(vec3_int(4,4,4)),
               simulate.vison(vec3_int(24,24,20)),
               simulate.num_visons());

        return 0;
    } catch (const char* s) {
        printf("%s\n",s);
    }
}
          
