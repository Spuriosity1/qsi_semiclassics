#pragma once
#include <direction.hh>



// calculates g0,1,2,3 from given B and J_\pm
// Uses the expression given in Equation (4) of
// arXiv:2312.11641v1
// Gives the thing in units of 1.5 Jpm^3
void calc_g_values(double* g, double Jpm, const double B[]){
    // safe code is for cowards
    for (uint8_t mu=0; mu<4; mu++){
        double z_dot_B = -1* (
            direction::axis[mu][2][0]*B[0]+
            direction::axis[mu][2][1]*B[1]+
            direction::axis[mu][2][2]*B[2]
        );
        
        g[mu] = 1.5*Jpm*Jpm*Jpm + 1.25*Jpm*Jpm*z_dot_B*z_dot_B;
        g[mu] /= fabs(1.5*Jpm*Jpm*Jpm);
    }
}

void calc_g_defect_values(double* g_defect, double Jpm, const double B[]){
    for (uint8_t mu=0; mu<4; mu++){
        double z_dot_B = -1*(
            direction::axis[mu][2][0]*B[0]+
            direction::axis[mu][2][1]*B[1]+
            direction::axis[mu][2][2]*B[2]
        );
        
        g_defect[mu] = -8.*Jpm*Jpm*z_dot_B;
        g_defect[mu] /= fabs(1.5*Jpm*Jpm*Jpm);
    }
}

