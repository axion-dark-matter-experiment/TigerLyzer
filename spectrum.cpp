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

class recurse_mean {

  public:
    recurse_mean() {}
    recurse_mean( double initial_sample ) {
        held_average = initial_sample;
        counter ++;
    }

    double operator()(double current_sample) {
        if ( counter == 0 ) {
            counter ++;
            held_average = current_sample;
        } else {
            counter ++;
            held_average = (counter - 1)/counter*held_average + current_sample/counter;
        }

        return held_average;
    }

    double uncertainity (double current_sample) {
        return sqrt(pow((counter - 1)/counter*held_average,2.0) + pow((current_sample/counter),2.0));
    }

  private:
    double held_average = 0.0;
    uint counter = 0;
};

//grand_g_prediction.power[i]*(pow(get_axion_KSVZ_coupling(grand_g_prediction.getBinMidFreq(i)),2.0))
//grand_g_prediction.uncertainty[i]*(pow(get_axion_KSVZ_coupling(grand_g_prediction.getBinMidFreq(i)),2.0))

//grand_spectrum.sa_power_list.at(i) *= pow( KSVZ_axion_coupling( grand_spectrum.bin_mid_freq(i). , 2.0);
//grand_spectrum.uncertainties.at(i) *= pow( KSVZ_axion_coupling( grand_spectrum.bin_mid_freq(i). , 2.0);

SingleSpectrum Spectrum::GrandSpectrum() {

    auto grand_spectrum = BlankGrandSpectrum();
    uint g_size = grand_spectrum.size();


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

    return grand_spectrum;
}

template<typename T> inline T positive_part ( T& x ) {
    return ( x < 0 )? 0 : x;
}


//    thelimits.power[i]=x+2.0*grand_g_prediction.uncertainty[i];

//    thelimits.power[i] = sqrt(thelimits.power[i]) *get_axion_KSVZ_coupling(thelimits.getBinMidFreq(i));
//    thelimits.uncertainty[i] = get_axion_KSVZ_coupling(thelimits.getBinMidFreq(i));

SingleSpectrum Spectrum::Limits() {

    auto g_spectrum = GrandSpectrum();

    for (uint i = 0; i < g_spectrum.size() ; i++ ) {

        //double x=grand_g_prediction.power[i];
        //if(x<0) x=0;
        //thelimits.power[i]=x+2.0*grand_g_prediction.uncertainty[i];

        //thelimits.power[i] = sqrt(thelimits.power[i]) *get_axion_KSVZ_coupling(thelimits.getBinMidFreq(i));
        //thelimits.uncertainty[i] =get_axion_KSVZ_coupling(thelimits.getBinMidFreq(i));


        double g_power = positive_part( g_spectrum.sa_power_list[i]);

        g_spectrum.sa_power_list[i] = g_power + 2.0*g_spectrum.uncertainties[i];
        g_spectrum.sa_power_list[i] = sqrt( g_spectrum.sa_power_list[i] )* KSVZ_axion_coupling( g_spectrum.bin_mid_freq(i) );

        g_spectrum.uncertainties[i] = KSVZ_axion_coupling(g_spectrum.bin_mid_freq(i));
    }

    g_spectrum.rebin( 600 );

    return g_spectrum;
}

inline double axion_coupling_power( double g_spec_power, double g_spec_mid_freq ) {
    return  g_spec_power*pow(KSVZ_axion_coupling(g_spec_mid_freq),2.0);
}

inline double axion_coupling_uncertainity ( double g_spec_uncertainity, double g_spec_mid_freq ) {
    return g_spec_uncertainity*pow(KSVZ_axion_coupling(g_spec_mid_freq),2.0);
}

SingleSpectrum Spectrum::GSquaredPrediction() {

    auto g_spec = GrandSpectrum();

    for( uint i = 0 ; i < g_spec.size() ; i++ ) {

        double power = g_spec.sa_power_list.at(i);
        double uncertainty = g_spec.uncertainties.at(i);
        double mid_freq = g_spec.bin_mid_freq(i);

        double ax_power = axion_coupling_power ( power, mid_freq );
        g_spec.sa_power_list.at(i) = ax_power;
        double ax_uncertainity = axion_coupling_power ( uncertainty, mid_freq );
        g_spec.uncertainties.at(i) = ax_uncertainity;
    }

    return g_spec;
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
