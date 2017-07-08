// Header for this file
#include "singlespectrum.h"
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
#include <mutex>       //std::mutex
// Boost Headers
#include <boost/algorithm/string.hpp>  //split() and is_any_of for parsing .csv files
#include <boost/lexical_cast.hpp>  //lexical cast (unsurprisingly)
#include <dirent.h>
// Miscellaneous Headers
#include <omp.h>  //OpenMP pragmas
#include "/home/bephillips2/gnuplot-iostream/gnuplot-iostream.h"
//Project Specific Headers
#include "physicsfunctions.h"


SingleSpectrum::SingleSpectrum(std::string raw_data) {
    //Load relevent parameters from string
    ParseRawData(raw_data);

    //convert from natives units of dBm to absolute power in watts
    dBmToWatts();
}

SingleSpectrum::SingleSpectrum(uint size) {
    sa_power_list = std::vector<double> (size, 0.0);
    uncertainties = std::vector<double> (size, 0.0);
}

SingleSpectrum::SingleSpectrum(uint size, double min_freq, double max_freq) {

    center_frequency = (max_freq - min_freq)/2.0 + min_freq;
    frequency_span = (max_freq - min_freq);

    sa_power_list = std::vector<double> (size, 0.0);
    uncertainties = std::vector<double> (size, 0.0);
}

SingleSpectrum::~SingleSpectrum() {
    sa_power_list.clear();
}

uint SingleSpectrum::size() {
    return sa_power_list.size();
}

SingleSpectrum &SingleSpectrum::operator*=(double scalar) {

    for(unsigned int i=0; i<sa_power_list.size(); i++) {
        sa_power_list[i]*=scalar;
        uncertainties[i]*=scalar;
    }

    return *this;
}

SingleSpectrum &SingleSpectrum::operator+=(double scalar) {

    for(unsigned int i=0; i<sa_power_list.size(); i++) {
        sa_power_list[i]+=scalar;
    }
    return *this;
}


bool operator== (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b) {

    return ( spectra_a.sa_power_list == spectra_b.sa_power_list);
}

bool operator!= (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b) {

    return ( spectra_a.sa_power_list != spectra_b.sa_power_list);
}

std::ostream& operator << (std::ostream& stream, SingleSpectrum& spectrum) {

    for(unsigned int i=0; i < spectrum.size(); i++) {
        stream << spectrum.sa_power_list[i]\
               << ","\
               << spectrum.uncertainties[i]\
               << std::endl;
    }

    return stream;
}

std::ofstream& operator << (std::ofstream& stream, SingleSpectrum& spectrum) {

    double delta_f = spectrum.max_freq() - spectrum.min_freq();
    double min_f = spectrum.min_freq();
    double n = static_cast<double>( spectrum.size() );

    for(unsigned int i=0; i<spectrum.size(); i++) {

        double freq = delta_f*( i / n ) + min_f;
        stream << freq << "," << spectrum.sa_power_list[i] << std::endl;
    }

    return stream;
}

SingleSpectrum operator+ (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b) {

    if( spectra_a.size() != spectra_b.size() ) {
        std::string err_mesg = "Spectra are not the same size ";
        err_mesg += boost::lexical_cast<std::string> (spectra_a.size());
        err_mesg += " vs ";
        err_mesg += boost::lexical_cast<std::string> (spectra_b.size());
        throw std::length_error(err_mesg);
    }

    auto spectra_c = spectra_a;

    std::transform(spectra_c.sa_power_list.begin(),\
                   spectra_c.sa_power_list.end(),\
                   spectra_b.sa_power_list.begin(),\
                   spectra_c.sa_power_list.begin(),\
                   std::plus<double>());

    return spectra_c;
}

SingleSpectrum operator+ (SingleSpectrum& spectra_a, std::vector<double>& spectra_b) {

    if( spectra_a.size() != spectra_b.size() ) {
        std::string err_mesg = "Spectra and vector are not the same size ";
        err_mesg += boost::lexical_cast<std::string> (spectra_a.size());
        err_mesg += " vs ";
        err_mesg += boost::lexical_cast<std::string> (spectra_b.size());
        throw std::length_error(err_mesg);
    }

    auto spectra_c = spectra_a;

    spectra_c.sa_power_list.resize(spectra_b.size());

    std::transform(spectra_c.sa_power_list.begin(),\
                   spectra_c.sa_power_list.end(),\
                   spectra_b.begin(),\
                   spectra_c.sa_power_list.begin(),\
                   std::plus<double>());

    return spectra_c;
}

SingleSpectrum operator+ (SingleSpectrum& spectra_a, double scalar) {

    auto spectra_b = spectra_a;

    for(unsigned int i=0; i<spectra_b.sa_power_list.size(); i++) {
        spectra_b.sa_power_list[i] += scalar;
        spectra_b.uncertainties[i] += scalar;
    }

    return spectra_b;
}

SingleSpectrum operator- (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b) {

    if( spectra_a.size() != spectra_b.size() ) {
        std::string err_mesg = "Spectra are not the same size ";
        err_mesg += boost::lexical_cast<std::string> (spectra_a.size());
        err_mesg += " vs ";
        err_mesg += boost::lexical_cast<std::string> (spectra_b.size());
        throw std::length_error(err_mesg);
    }

    auto spectra_c = spectra_a;

    std::transform(spectra_c.sa_power_list.begin(),\
                   spectra_c.sa_power_list.end(),\
                   spectra_b.sa_power_list.begin(),\
                   spectra_c.sa_power_list.begin(),\
                   std::minus<double>());

    return spectra_c;
}

SingleSpectrum operator- (SingleSpectrum& spectra_a, std::vector<double>& spectra_b) {

    if( spectra_a.size() != spectra_b.size() ) {
        std::string err_mesg = "Spectra and vector are not the same size ";
        err_mesg += boost::lexical_cast<std::string> (spectra_a.size());
        err_mesg += " vs ";
        err_mesg += boost::lexical_cast<std::string> (spectra_b.size());
        throw std::length_error(err_mesg);
    }

    auto spectra_c = spectra_a;

    std::transform(spectra_c.sa_power_list.begin(),\
                   spectra_c.sa_power_list.end(),\
                   spectra_b.begin(),\
                   spectra_c.sa_power_list.begin(),\
                   std::minus<double>());

    return spectra_c;
}

SingleSpectrum operator- (SingleSpectrum& spectra_a, double scalar) {

    auto spectra_b = spectra_a;

    for(unsigned int i=0; i<spectra_b.sa_power_list.size(); i++) {
        spectra_b.sa_power_list[i] -= scalar;
        spectra_b.uncertainties[i] -= scalar;
    }

    return spectra_b;
}

SingleSpectrum operator* (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b) {

    if( spectra_a.size() != spectra_b.size() ) {
        std::string err_mesg = "Spectra are not the same size ";
        err_mesg += boost::lexical_cast<std::string> (spectra_a.size());
        err_mesg += " vs ";
        err_mesg += boost::lexical_cast<std::string> (spectra_b.size());
        throw std::length_error(err_mesg);
    }

    auto spectra_c = spectra_a;

    for(unsigned int i=0; i<spectra_a.size(); i++) {
        spectra_c.sa_power_list[i] =\
                                    spectra_a.sa_power_list[i] *\
                                    spectra_b.sa_power_list[i];
    }

    return spectra_c;
}

SingleSpectrum operator* (SingleSpectrum& spectra_a, std::vector<double>& spectra_b) {

    if( spectra_a.size() != spectra_b.size() ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nSpectra and vector are not the same size ";
        err_mesg += boost::lexical_cast<std::string> (spectra_a.size());
        err_mesg += " vs ";
        err_mesg += boost::lexical_cast<std::string> (spectra_b.size());
        throw std::length_error(err_mesg);
    }

    auto spectra_c = spectra_a;

    for(unsigned int i=0; i<spectra_a.size(); i++) {
        spectra_c.sa_power_list[i] =\
                                    spectra_a.sa_power_list[i] *\
                                    spectra_b[i];
    }

    return spectra_c;
}

SingleSpectrum operator* (SingleSpectrum& spectra_a, double scalar) {

    auto spectra_b = spectra_a;

    for(unsigned int i=0; i<spectra_b.sa_power_list.size(); i++) {
        spectra_b.sa_power_list[i] *= scalar;
        spectra_b.uncertainties[i] *= scalar;
    }

    return spectra_b;
}

void SingleSpectrum::dBmToWatts() {
    if( current_units != Units::dBm ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nSpectra must be in units of dBm.";
        throw std::invalid_argument(err_mesg);
    }

    for (auto& power : sa_power_list ) {
        power = dbm_to_watts(power);
    }

    current_units = Units::Watts;
}

/*!
 * \brief Convert from usings of watts to units of watts above noise (i.e. excess power)
 */
void SingleSpectrum::WattsToExcessPower() {

    if( current_units != Units::Watts ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += ": Spectra is not in units of Watts.";
        throw std::invalid_argument(err_mesg);
    }

    double noise_power = power_per_bin( noise_temperature, bin_width() );

    //Yes, these functions really should be handled by the +=
    //and *= operators- unfortunately it is not so easy to implement
    //them inside a member function
    double mean_val = mean();

    for(uint i=0; i<sa_power_list.size(); i++) {
        sa_power_list.at(i) *= ( noise_power/ mean_val );
    }

    for(uint i=0; i<sa_power_list.size(); i++) {
        sa_power_list.at(i) += ( -noise_power );
    }

    PopulateUncertainties( 32 );

    current_units = Units::ExcessPower;
}


//double SingleSpectrum::kszv_power_per_bin( double freq_mhz ) {
//    double freq_ghz = freq_mhz/1e3;

//    double max_power = max_ksvz_power(effective_volume, b_field, freq_ghz);
//    double off_center_power = Q*lorentzian( center_frequency, freq_mhz, Q );

//    return max_power*off_center_power;
//}
//double lorentzian (double f0, double omega, double Q )
void SingleSpectrum::LorentzianWeight() {
    if( current_units != Units::ExcessPower ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nSpectra must be in units of excess power.";
        throw std::invalid_argument(err_mesg);
    }

    for(uint i=0; i < size(); i++) {

        double frequency = bin_mid_freq(i);

        sa_power_list.at(i) /= lorentzian( center_frequency, frequency, Q );
        uncertainties.at(i) /= lorentzian( center_frequency, frequency, Q );
    }

}

void SingleSpectrum::KSVZWeight() {

    if( current_units != Units::ExcessPower ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nSpectra must be in units of excess power.";
        throw std::invalid_argument(err_mesg);
    }

    for(uint i=0; i < size(); i++) {

        double frequency = bin_mid_freq(i);

        sa_power_list.at(i) /= max_ksvz_power( effective_volume, b_field, frequency, Q );
        uncertainties.at(i) /= max_ksvz_power( effective_volume, b_field, frequency, Q );
    }

    current_units = Units::AxionPower;

}

std::string SingleSpectrum::units() {
    switch( current_units ) {
    case Units::dBm:
        return "dBm";
        break;
    case Units::ExcessPower:
        return "Excess Power in Cavity (Watts)";
        break;
    case Units::AxionPower:
        return "Power Deposited by Axion";
        break;
    case Units::Watts:
        return "Watts";
        break;
    case Units::ExclLimit90:
        return "G a gamma gamma ( GeV ^ -1 )";
    default:
        return "";
        break;
    }
}

void SingleSpectrum::PopulateUncertainties( uint rebin_size ) {

    uncertainties.clear();

    double noise_power = power_per_bin( noise_temperature, bin_width() );
    double uniform_uncertainty = noise_power / sqrt( number_of_averages*rebin_size );

    uncertainties = std::vector<double> ( size() , uniform_uncertainty);
}

void SingleSpectrum::InitialBin ( uint bin_points ) {

    if( current_units != Units::Watts ) {
        std::string err_mesg = "Spectra needs to be in Watts before initial binning.";
        throw std::invalid_argument(err_mesg);
    }

    // choose bin size (approximate)
    // bin_window ~ # points / (0.5 * axion width / (frequency span / # points))
    // bins overlap in half bin size increments
    uint bin_window = static_cast<uint>(bin_points / 2);

    double prev_sum = 0; // sum of previous FFT points
    double curr_sum = 0; // running sum of current FFT points

    std::vector<double> rebinned_power_list;
    std::vector<double> rebinned_uncertainties;

    for ( uint fft_i = 0 ; fft_i < size() ; fft_i ++ ) {
        curr_sum += sa_power_list[fft_i];

        if ((fft_i % bin_window) == bin_window - 1) {
            // average over all the points in the current and previous sums
            if (fft_i != bin_window - 1) {
                double new_power = (curr_sum + prev_sum) / (bin_window*2);

                rebinned_power_list.push_back(new_power);
                rebinned_uncertainties.push_back(new_power);
            }
            prev_sum = curr_sum;
            curr_sum = 0;
        }
    }

    sa_power_list = rebinned_power_list;
    uncertainties = rebinned_uncertainties;
}

double SingleSpectrum::bin_width() {
    return frequency_span/static_cast<double>(size());
}

double SingleSpectrum::bin_start_freq(uint idx) {
    double freq_start = center_frequency - 0.5*frequency_span;
    double dub_idx = static_cast<double>(idx);
    double dub_size = static_cast<double>( size() );

    return freq_start+dub_idx*frequency_span/dub_size;
}

uint SingleSpectrum::bin_at_frequency(double frequency) {

    double min_freq = center_frequency - frequency_span/2.0;
    double max_freq = center_frequency + frequency_span/2.0;

    if( frequency < min_freq || frequency > max_freq ) {
        std::string err_mesg = "Requested frequency of ";
        err_mesg += boost::lexical_cast<std::string>(frequency);
        err_mesg += " is outside of spectrum range: ";
        err_mesg += boost::lexical_cast<std::string>(min_freq);
        err_mesg += "to " + boost::lexical_cast<std::string>(max_freq);

        throw std::out_of_range(err_mesg);
    }

    double bin_number = (frequency - min_freq)/frequency_span;
    bin_number *= static_cast<double>(size());
    bin_number = floor(bin_number);
    return static_cast<uint>( bin_number );
}

double SingleSpectrum::bin_mid_freq(uint idx) {

    return bin_start_freq( idx )+0.5*bin_width();
}


// rebin the spectrum by making bins of npoints, keeping the most conservative
// power and uncertainty
void SingleSpectrum::rebin( uint points_per_bin ) {

//    uint new_size = static_cast<uint>( floor( size()/ points_per_bin ) );

    std::vector<double> nu_power_list;
    std::vector<double> nu_uncertainties;

    double min_power= *std::min_element(sa_power_list.begin(), sa_power_list.end());
    double min_uncertainty = *std::min_element(uncertainties.begin(), uncertainties.end());

    double highest_power = min_power;
    double highest_uncertainty = min_uncertainty;

    for ( uint i = 0 ; i < size() ; i ++ ) {
        highest_power = std::max(highest_power, sa_power_list[i]);
        highest_uncertainty = std::max(highest_uncertainty, uncertainties[i]);

        if ((i % points_per_bin) == points_per_bin - 1) {

            nu_power_list.push_back( highest_power );
            nu_uncertainties.push_back( highest_uncertainty );

            highest_power = min_power;
            highest_uncertainty = min_uncertainty;
        }
    }

    sa_power_list.clear();
    uncertainties.clear();

    sa_power_list = nu_power_list;
    uncertainties = nu_uncertainties;
}

void SingleSpectrum::chop_bins( uint start_chop, uint end_chop ) {

    uint distance = sa_power_list.size() - end_chop;

    frequency_span *= (size() - start_chop - end_chop)/size();

    auto power_start = sa_power_list.begin() + start_chop;
    auto power_stop = sa_power_list.begin() + distance;

    sa_power_list = std::vector<double> (power_start, power_stop);

    auto uncertainty_start = uncertainties.begin() + start_chop;
    auto uncertainty_stop = uncertainties.begin() + distance;

    uncertainties = std::vector<double> (uncertainty_start, uncertainty_stop);

}

uint SingleSpectrum::num_lines(std::string raw_data) {
    return std::count(raw_data.begin(), raw_data.end(), '\n');
}

template <typename T>
bool inline check_map ( std::map<std::string, T> check_map , std::vector<std::string> check_keys ) {

    for (const auto& key : check_keys ) {
        if( check_map.count( key ) ) {
            continue;
        } else {
            std::cerr << "Could not find " << key << " in map. " <<std::endl;
            return false;
        }
    }

    return true;
}

void SingleSpectrum::FillFromHeader(std::map<std::string, double> header) {
    std::vector<std::string> search_terms = {"sa_span",
                                             "fft_length",
                                             "effective_volume",
                                             "bfield",
                                             "noise_temperature",
                                             "sa_averages",
                                             "Q",
                                             "actual_center_freq",
                                             "fitted_hwhm",
                                             "bfield"
                                            };

    if( !check_map ( header, search_terms ) ) {
        std::string err_mesg = "Insufficient information to build spectra.";
        throw std::invalid_argument(err_mesg);
        return;
    }

    center_frequency = header["actual_center_freq"];
    frequency_span = header["sa_span"];
    effective_volume = header["effective_volume"];
    noise_temperature = header["noise_temperature"];
    Q = header["Q"];
    number_of_averages = header["sa_averages"];
    fft_points = header["fft_length"];
    b_field = header["bfield"];
}

void SingleSpectrum::ParseRawData(std::string raw_data) {
    uint lines = num_lines(raw_data);
    sa_power_list.reserve( lines );

    std::istringstream data_stream(raw_data);

    std::map<std::string, double> header;

    for (uint i = 0; i < 11; i++) {
        std::string input;
        std::getline(data_stream, input);

        std::vector<std::string> strs;
        boost::split(strs, input, boost::is_any_of(";"));
        if (strs.size() >= 2) {
            std::string name = strs.at(0);
            double val = boost::lexical_cast<double>(strs.at(1));

            header[name] = val;
        }
    }

    FillFromHeader( header );

    for (uint i = 12; i < lines; i++) {

        std::string input;
        std::getline(data_stream, input);
        double val = boost::lexical_cast<double>(input);

        sa_power_list.push_back(val);

    }

}

double SingleSpectrum::min_freq() {
    return ( center_frequency - 0.5*frequency_span);
}

double SingleSpectrum::max_freq() {
    return ( center_frequency + 0.5*frequency_span);
}

double SingleSpectrum::sum(std::vector<double>& data_list, double exponent) {

    double tot = 0;

    for ( uint i = 0 ; i < data_list.size(); i ++) {
        tot += pow( data_list[i], exponent );
    }

    return tot;
}

double SingleSpectrum::mean(std::vector<double>& data_list) {
    //compute mean value of data set
    double sum_x=sum(data_list,1.0);
    double n=data_list.size();
    return sum_x/n;
}

double SingleSpectrum::mean() {
    return mean(sa_power_list);
}

double SingleSpectrum::std_dev(std::vector<double> &data_list) {

    //compute mean value of data set
    double sum_x=sum(data_list,1.0);
    double n=data_list.size();
    double mean = sum_x/n;

    //compute variance taking into account Bessel's correction i.e. n/(n-1)
    double sum_x2=sum(data_list,2.0);
    double sigma_sqr=sum_x2/(n-1.0)-n/(n-1.0)*pow(mean,2.0);

    //return square root of variance
    return sqrt(sigma_sqr);
}

double SingleSpectrum::std_dev() {
    return std_dev(sa_power_list);
}

double SingleSpectrum::norm() {
    return sqrt(sum(sa_power_list,2.0));
}
