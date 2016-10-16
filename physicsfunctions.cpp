// Header for this file
#include "physicsfunctions.h"
// C System-Headers
//
// C++ System headers
#include <cmath>       //sqrt, abs, M_PI
// Boost Headers
//
// Miscellaneous Headers
//
//Project Specific Headers
//


#define ALPHA 7.2973525376e-3 // strong coupling constant = 1/137
#define H 4.13566766225e-15 //h in eV s
#define G_KSVZ 0.97 // eV
#define KB 1.3806488e-23 //Watts / Hz / K


double axion_width ( double frequency ) {
    return frequency*10.0e-6/2.0;
}

double KSVZ_axion_coupling( double frequency ) {
    //compute mass in eV
    double mass_ev = frequency*H*1e6;
    //compute coupling in GeV^-1
    return 1e-7*(mass_ev/0.62)*(ALPHA*G_KSVZ/M_PI);
}
double lorentzian (double f0, double omega, double Q ) {

    double gamma = omega/(2.0*Q);
    return pow(gamma,2.0)/( pow( (omega-f0), 2.0 )+pow(gamma,2.0) );

}

double max_ksvz_power( double effective_volume, double b_field, double frequency, double Q) {
    return 2.278e-33*pow(b_field, 2.0)*effective_volume*frequency*Q;
}

double power_per_bin( double noise_temperature, double bin_width ) {
    return KB*noise_temperature*bin_width*1e6;
}

double dbm_to_watts ( double power_dbm ) {
    return pow10(power_dbm / 10.0)/1000;
}

#undef ALPHA
#undef H
#undef G_KSVZ
#undef KB
