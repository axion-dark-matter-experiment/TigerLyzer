// Header for this file
#include "physicsfunctions.h"
// C System-Headers
//
// C++ System headers
#include <cmath>       //sqrt, abs, M_PI
// Boost Headers
#include <boost/algorithm/string.hpp>  //split() and is_any_of for parsing .csv files
#include <boost/lexical_cast.hpp>  //lexical cast (unsurprisingly)
#include <dirent.h>
// Miscellaneous Headers
#include <omp.h>  //OpenMP pragmas
#include "/home/bephillips2/gnuplot-iostream/gnuplot-iostream.h"
//Project Specific Headers
#include "physicsfunctions.h"


#define ALPHA 7.2973525376e-3 // strong cooupling constant = 1/137
#define H 4.13566766225e-15 //h in eV s
#define G_KSVZ 0.97 // eV
#define KB 1.3806488e-23 //Watts / Hz / K


double axion_width ( double frequency ) {
    return frequency*10.0e-6/2.0;
}

//See Ed Daw's theses, page 23, eq.2.25
double KSVZ_axion_coupling( double frequency ) {
    //compute mass in GeV^-1
    double mass_ev = frequency*H*1e6;
    //compute coupling in GeV^-1
    return 1e-7*(mass_ev/0.62)*(ALPHA*G_KSVZ/M_PI);
}

//Aply to uncertainties and powers
//double freq_mhz = bin_mid_frequency();

double estimate_G_2 ( double freq_mhz ) {
    return pow( KSVZ_axion_coupling( freq_mhz ), 2.0);
}

double lorentzian (double f0, double omega, double Q ) {

    double gamma = omega/(2.0*Q);
    return pow(gamma,2.0)/( pow( (omega-f0), 2.0 )+pow(gamma,2.0) );

}

double max_ksvz_power( double effective_volume, double b_field, double frequency, double Q) {
    return 2.278e-33*pow(b_field, 2.0)*effective_volume*frequency*Q;
}

//double noise_power=kB*runs[onrun].noise_temperature *(thespectrum.getBinWidth()*1e6);

//Noise power per bin, see Daw These pg. 72 eq. 4.10
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
