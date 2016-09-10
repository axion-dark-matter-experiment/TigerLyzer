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

SingleSpectrum::SingleSpectrum(std::string raw_data) {
    //Load relevent parameters from string
    ParseRawData(raw_data);
    //convert from natives units of dBm to absolute power in watts
    dBmToWatts();
    //perform initial binning using a bin width of 72 points.
//    InitialBin( 72 );
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

    for(unsigned int i=0; i<spectrum.size(); i++) {
        stream << spectrum.sa_power_list[i] << std::endl;
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

/*****************************************************************************/
/**********Functions from legacy code that are being updated *****************/


#define ALPHA 7.2973525376e-3 // who knows?
#define H 6.58211899e-16*2*M_PI //h in eV s
#define G_KSVZ 0.97 // eV
#define KB 1.3806488e-23 //Watts / Hz / K


double axion_width ( double frequency ) {
    return frequency*10.0e-6/2.0;
}

//Don't ask, I don't know
double KSVZ_axion_coupling( double frequency ) {
    //compute mass in GeV^-1
    double mass_ev = frequency*H*1e6;
    //compute coupling in GeV^-1
    double axion_coupling = 1e-7*(mass_ev/0.62)*(ALPHA*G_KSVZ/M_PI);

    return axion_coupling;
}

//Aply to uncertainties and powers
//double freq_mhz = bin_mid_frequency();

double estimate_G_2 ( double freq_mhz ) {
    return pow( KSVZ_axion_coupling( freq_mhz ), 2.0);
}


/*!
 * \brief Lorentz Line Shape
 * \param f0, the center frequency
 * \param frequency of interest
 * \param Q, quality factor
 * \return value of lorenzian at omega
 */
double lorentzian (double f0, double omega, double Q ) {

    double gamma=omega/(2.0*Q);
    return pow(gamma,2.0)/( (omega-f0)*(omega-f0)+pow(gamma,2.0));

}

//get the axion photon conversion power for a Q of 1
//veff is in cm^3
//B is in Tesla
//f is in GHz
//assumes dark matter density is 0.45 GeV/cm^3
double max_ksvz_power( double effective_volume, double b_field, double frequency) {
    return 2.2e-23*pow(b_field, 2.0)*effective_volume*frequency;
}

//double noise_power=kB*runs[onrun].noise_temperature *(thespectrum.getBinWidth()*1e6);

double power_per_bin( double noise_temperature, double bin_width ) {
    return KB*noise_temperature*bin_width*1e6;
}

double inline dbm_to_watts ( double power_dbm ) {
    return pow10(power_dbm / 10.0);
}


void plot ( SingleSpectrum& spec, std::string plot_title ) {

    Gnuplot gp;

    std::vector<double> x_vals;

    double min_freq = spec.center_frequency - spec.frequency_span/2.0;
    double span = spec.frequency_span;
    double spectrum_points = static_cast<double> (spec.size());

    for( uint i = 0 ; i < spec.size() ; i ++ ) {
        double frequency = min_freq + span*i/spectrum_points;
        x_vals.push_back(frequency);
    }

    gp << "set title '"+plot_title+"'\n";
    gp << "plot '-' using 1:2 with lines title 'Power'\n";
    gp.send( boost::make_tuple( x_vals, spec.sa_power_list ) );
}

void plot ( SingleSpectrum& spec, uint num_plot_points, std::string plot_title ) {

    Gnuplot gp;

    if( num_plot_points > spec.size() ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nRequested number of points (";
        err_mesg += boost::lexical_cast<std::string>( num_plot_points );
        err_mesg += ") exceeds size of spectrum (";
        err_mesg += boost::lexical_cast<std::string>( spec.size() ) +")";
        throw std::out_of_range(err_mesg);
    }

    uint pivot = static_cast<uint> (spec.size() / num_plot_points);

    double min_freq = spec.center_frequency - spec.frequency_span/2.0;
    double span = spec.frequency_span;
    double spectrum_points = static_cast<double> (spec.size());

    std::vector<uint> pivot_indices;
    for( uint i = 0 ; i < spec.size() ; i ++ ) {
        if ( i%pivot == 0 ) {
            pivot_indices.push_back(i);
        }
    }

    std::vector<double> x_vals;
    for ( const auto& idx : pivot_indices ) {
        double frequency = min_freq + span*idx/spectrum_points;
        x_vals.push_back(frequency);
    }

    std::vector<double> y_vals;
    for ( const auto& idx : pivot_indices ) {
        y_vals.push_back( spec.sa_power_list.at(idx) );
    }

    std::vector<double> delta_y_vals;
    for ( const auto& idx : pivot_indices ) {
        delta_y_vals.push_back( spec.uncertainties.at(idx) );
    }

    gp << "set title '"+plot_title+"'\n";
    gp << "plot '-' using 1:2:3 with yerror title 'Power'\n";
    gp.send( boost::make_tuple(x_vals, y_vals, delta_y_vals) );
}

void SingleSpectrum::dBmToWatts() {
    if( current_units == Units::Watts ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nSpectra is already in units of Watts.";
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

    PopulateUncertainties();

    current_units = Units::ExcessPower;
}


double SingleSpectrum::kszv_power_per_bin( double freq_mhz ) {
    double freq_ghz = freq_mhz/1e3;

    double max_power = max_ksvz_power(effective_volume, b_field, freq_ghz);
    double off_center_power = Q*lorentzian( center_frequency, freq_mhz, Q );

    return max_power*off_center_power;
}

void SingleSpectrum::KSVZWeight() {

    if( current_units != Units::ExcessPower ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nSpectra must be in units of excess power.";
        throw std::invalid_argument(err_mesg);
    }

    for(uint i=0; i < size(); i++) {

        double frequency = bin_mid_freq(i);

        sa_power_list.at(i) /= kszv_power_per_bin( frequency );
        uncertainties.at(i) /= kszv_power_per_bin( frequency );
    }


}

//#undef ALPHA
//#undef H
//#undef G_KSVZ
//#undef KB

/*****************************************************************************/
/*****************************************************************************/

std::string SingleSpectrum::units() {
    switch( current_units ) {
    case Units::dBm:
        return "dBm";
        break;
    case Units::ExcessPower:
        return "Excess Power";
        break;
    case Units::Watts:
        return "Watts";
        break;
    default:
        return "";
        break;
    }
}

void SingleSpectrum::PopulateUncertainties() {

    uncertainties.clear();

    double noise_power = power_per_bin( noise_temperature, bin_width() );
    double uniform_uncertainty = noise_power / sqrt( number_of_averages * bin_width() );

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

    double prev_sum = 0; // sum of previous 36 FFT points
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
        highest_power = fmax(highest_power, sa_power_list[i]);
        highest_uncertainty = fmax(highest_uncertainty, uncertainties[i]);

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

double SingleSpectrum::sum(std::vector<double>& data_list,double exponent) {

    double tot=0;

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

/*
 *     //now make a grand spectrum
    double total_span=freq_max-freq_min;
    unsigned int total_bins=(unsigned int)(floor(total_span/bin_size))-1;

    Spectrum grand_g_prediction(total_bins);

    grand_g_prediction.freq_start=freq_min;
    grand_g_prediction.freq_span=total_span;

    grand_g_prediction.zero();

    for(unsigned int i=0; i<grand_g_prediction.length; i++) {

        double f=grand_g_prediction.getBinMidFreq(i);

        for(unsigned int ons=0; ons<g_predictions.size(); ons++) {

            unsigned int binno=g_predictions[ons].getBinAtFreq(f);

            if(binno==INT_MAX) continue;

            if(grand_g_prediction.uncertainty[i]==0) {

                grand_g_prediction.power[i]=g_predictions[ons].power[binno];
                grand_g_prediction.uncertainty[i]=g_predictions[ons].uncertainty[binno];

            } else {

                double a=grand_g_prediction.power[i];
                double b=g_predictions[ons].power[binno];
                double sa=grand_g_prediction.uncertainty[i];
                double sb=g_predictions[ons].uncertainty[binno];
                double taua=1.0/(sa*sa);
                double taub=1.0/(sb*sb);
                grand_g_prediction.power[i]=(taua*a+taub*b)/(taua+taub);
                grand_g_prediction.uncertainty[i]=sqrt(1.0/(taua+taub));

            }

        }
    }
 */

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

std::map<double, double> FreqIndexMap ( SingleSpectrum& spec) {

    std::map<double, double> freq_to_power;

    for( uint i = 0 ; i < spec.size() ; i++ ) {
        freq_to_power[ spec.bin_mid_freq(i) ] = i;
    }

    return freq_to_power;
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

//                    grand_spectrum.sa_power_list.at(i) = running_average (overlap_power);
//                    grand_spectrum.uncertainties.at(i) = running_average.uncertainity(overlap_uncertainity);

                    grand_spectrum.sa_power_list.at(i) = overlap_power_weight( current_power,\
                                                         overlap_power,\
                                                         current_uncertainty,\
                                                         overlap_uncertainity );

                    grand_spectrum.sa_power_list.at(i) = overlap_uncertainity_weight(\
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

SingleSpectrum Spectrum::Limits() {

    auto g_spectrum = GrandSpectrum();

    for (uint i = 0; i < g_spectrum.size() ; i++ ) {

        double g_power = positive_part( g_spectrum.sa_power_list[i]);

        g_spectrum.sa_power_list[i] = sqrt( g_power + 2.0*g_spectrum.uncertainties[i] )\
                                      * KSVZ_axion_coupling(g_spectrum.bin_mid_freq(i));

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
