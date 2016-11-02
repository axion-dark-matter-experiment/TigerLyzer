// Header for this file
#include "spectrum.h"
// C System-Headers
#include <termios.h> /* POSIX terminal control definitions */
#include <sys/ioctl.h>
#include <fcntl.h>   //fopen(),fclose()
#include <unistd.h>  //read(), write()
#include <stdio.h>

// C++ System headers
#include <vector>      //vector
#include <string>      //string
#include <fstream>     //iss* ofstream
#include <chrono>      // timing functions
#include <cmath>       //sqrt, abs
#include <iostream>    //cout
#include <typeinfo>    //typeid
#include <algorithm>   // transform, find, count, erase
#include <functional>  // plus/minus/multiplies
#include <utility>     //std::make_pair
#include <map>         //std::map
#include <mutex> //protect against concurrent access when using (unordered) parallel for loops

// Boost Headers
#include <boost/algorithm/string.hpp>  //split() and is_any_of for parsing .csv files
#include <boost/lexical_cast.hpp>  //lexical cast (unsurprisingly)
#include <dirent.h>

// Miscellaneous Headers
#include <omp.h>  //OpenMP pragmas
#include "/home/bephillips2/gnuplot-iostream/gnuplot-iostream.h"

//Project Specific Headers
#include "singlespectrum.h"
#include "physicsfunctions.h"


Spectrum::Spectrum() {}
Spectrum::~Spectrum() {}

uint Spectrum::size() {
    return spectra.size();
}


SingleSpectrum Spectrum::BlankGrandSpectrum () {

    if( size() == 0 ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nCannot build a Grand Spectrum- no spectra loaded.";
        throw std::invalid_argument(err_mesg);
    }

    uint total_bins = 0;
    std::vector<double> min_freqs;
    std::vector<double> max_freqs;

    std::mutex guard;

    #pragma omp parallel for
    for (uint i = 0; i < spectra.size() ; i ++ ) {
        guard.lock();
        total_bins += spectra.at(i).size();
        min_freqs.push_back( spectra.at(i).min_freq() );
        max_freqs.push_back( spectra.at(i).max_freq() );
        guard.unlock();
    }

    auto min = std::min_element(min_freqs.begin(),min_freqs.end());
    double min_frequency = static_cast<double>(*min);

    auto max = std::max_element(max_freqs.begin(),max_freqs.end());
    double max_frequency = static_cast<double>(*max);

    return SingleSpectrum (total_bins, min_frequency, max_frequency);
}

inline double overlap_power_weight( double power_a, double power_b, double delta_a, double delta_b ) {
    double tau_a = 1.0/pow( delta_a, 2.0 );
    double tau_b = 1.0/pow( delta_b, 2.0 );
    return (tau_a*power_a+tau_b*power_b)/(tau_a + tau_b);
}

inline double overlap_uncertainity_weight( double delta_a, double delta_b ) {
    double tau_a = 1.0/pow( delta_a, 2.0 );
    double tau_b = 1.0/pow( delta_b, 2.0 );
    return sqrt((1.0)/(tau_a + tau_b));
}

inline bool check_frequency( double to_check , SingleSpectrum& check_against) {
    return !( to_check > check_against.max_freq() || to_check < check_against.min_freq());
}


SingleSpectrum Spectrum::GrandSpectrum() {

    auto grand_spectrum = BlankGrandSpectrum();
    uint g_size = grand_spectrum.size();

    #pragma omp parallel for
    for(uint i=0; i< g_size; i++) {

        double g_frequency_at_i = grand_spectrum.bin_mid_freq(i);

        for(uint k = 0; k < size() ; k++ ) {

            if ( check_frequency(g_frequency_at_i, spectra.at(k)) ) {

                uint bin_number = spectra.at(k).bin_at_frequency( g_frequency_at_i );

                if ( grand_spectrum.sa_power_list.at(i) != 0.0 ) {

                    double current_power = grand_spectrum.sa_power_list.at(i);
                    double overlap_power = spectra.at(k).sa_power_list.at(bin_number);

                    double current_uncertainty = grand_spectrum.uncertainties.at(i);
                    double overlap_uncertainity = spectra.at(k).uncertainties.at(bin_number);

                    grand_spectrum.sa_power_list.at(i) = overlap_power_weight( current_power,\
                                                         overlap_power,\
                                                         current_uncertainty,\
                                                         overlap_uncertainity );

                    grand_spectrum.uncertainties.at(i) = overlap_uncertainity_weight(\
                                                         current_uncertainty,\
                                                         overlap_uncertainity );

                } else {

                    grand_spectrum.sa_power_list.at(i) = spectra.at(k).sa_power_list.at(bin_number);
                    grand_spectrum.uncertainties.at(i) = spectra.at(k).uncertainties.at(bin_number);

                }

            } else {
                continue;
            }

        }
    }

    grand_spectrum.current_units = Units::AxionPower;
    return grand_spectrum;
}

template<typename T> inline T positive_part ( T& x ) {
    return ( x < 0 )? 0 : x;
}

SingleSpectrum Spectrum::Limits() {

    auto g_spectrum = GrandSpectrum();

    for (uint i = 0; i < g_spectrum.size() ; i++ ) {

        double gc_power = positive_part( g_spectrum.sa_power_list[i]);
        double excl_90_watts = gc_power + 1.282*g_spectrum.uncertainties[i];

        g_spectrum.sa_power_list[i] = KSVZ_axion_coupling( g_spectrum.bin_mid_freq(i) )*sqrt( excl_90_watts );
        g_spectrum.uncertainties[i] = KSVZ_axion_coupling( g_spectrum.bin_mid_freq(i) );
    }

    g_spectrum.rebin( 600 );
    g_spectrum.current_units = Units::ExclLimit90;

    return g_spectrum;
}

inline double axion_coupling_power( double g_spec_power, double g_spec_mid_freq ) {
    return  g_spec_power*pow(KSVZ_axion_coupling(g_spec_mid_freq),2.0);
}

inline double axion_coupling_uncertainity ( double g_spec_uncertainity, double g_spec_mid_freq ) {
    return g_spec_uncertainity*pow(KSVZ_axion_coupling(g_spec_mid_freq),2.0);
}

//Operation does not seem to benefit from parallelism
void Spectrum::dBmToWatts() {

    for( auto& spec : spectra ) {
        spec.dBmToWatts();
    }
}

//Operation does not seem to benefit from parallelism
void Spectrum::WattsToExcessPower() {

    for( auto& spec : spectra ) {
        spec.WattsToExcessPower();
    }
}

//Operation does not seem to benefit from parallelism
void Spectrum::KSVZWeight() {

    for( auto& spec : spectra ) {
        spec.KSVZWeight();
    }

}

void Spectrum::LorentzianWeight() {

    for( auto& spec : spectra ) {
        spec.LorentzianWeight();
    }

}

//Operation does not seem to benefit from parallelism
SingleSpectrum Spectrum::at( uint idx ) {
    return spectra.at( idx );

}

Spectrum &Spectrum::operator+=(SingleSpectrum& spec) {

    spectra.push_back(spec);

    return *this;
}

Spectrum &Spectrum::operator-=(SingleSpectrum& spec) {

    for( uint i = 0; i < spectra.size() ; i++ ) {

        if( spectra[i] == spec ) {
            spectra.erase( spectra.begin() + i );
        }
    }

    return *this;
}
