#ifndef SINGLESPECTRUM_H
#define SINGLESPECTRUM_H

// C System-Headers
//
// C++ System headers
#include <vector>      //vector
#include <string>      //string
#include <fstream>     //iss* ofstream
#include <iostream>    //cout
// Boost Headers
//
// Miscellaneous Headers
//
//Project Specific Headers
#include "spectrum.h"


/*!
 * \brief Class to hold a single power spectrum and its associated parameters, such
 * as center frequency, frequency span, Q etc.
 *
 *Class is designed to 'look and feel' like a \f$ \mathcal{R}^n \f$ vector.
 * Spectra allows scalar multiplication and addition, along with 'vector'
 * multiplication and addition.
 */
class SingleSpectrum {

  public:

    SingleSpectrum(std::string raw_data);
    SingleSpectrum(uint size);
    SingleSpectrum(uint size, double min_freq, double max_freq);
    ~SingleSpectrum();

    SingleSpectrum &operator*=(double scalar);
    SingleSpectrum &operator+=(double scalar);

    friend bool operator== (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b);
    friend bool operator!= (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b);

    /*!
     * \brief Operator for printing spectrum contents to std::cout or similar
     * Only power data will be accessed and printed.
     *
     * \param stream object
     * \param spectrum to be printed
     * \return new stream object
     */
    friend std::ostream& operator<< (std::ostream& stream, SingleSpectrum& spectrum);
    /*!
     * \brief Operator for saving spectra to files
     * Unlike the ostream version, both the header and power data are accessed.
     *
     * \param stream object
     * \param spectrum to be saved
     * \return new stream object
     */
    friend std::ofstream& operator<< (std::ofstream& stream, SingleSpectrum& spectrum);

    friend SingleSpectrum operator* (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b);
    friend SingleSpectrum operator+ (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b);
    friend SingleSpectrum operator- (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b);

    friend SingleSpectrum operator* (SingleSpectrum& spectra_a, std::vector<double>& spectra_b);
    friend SingleSpectrum operator+ (SingleSpectrum& spectra_a, std::vector<double>& spectra_b);
    friend SingleSpectrum operator- (SingleSpectrum& spectra_a, std::vector<double>& spectra_b);

    friend SingleSpectrum operator* (SingleSpectrum& spectra_a, double scalar);
    friend SingleSpectrum operator+ (SingleSpectrum& spectra_a, double scalar);
    friend SingleSpectrum operator- (SingleSpectrum& spectra_a, double scalar);

    friend void plot ( SingleSpectrum& spec, std::string plot_title );
    friend void plot ( SingleSpectrum& spec, uint num_plot_points, std::string plot_title );

    friend void GaussianFilter ( SingleSpectrum& spec, uint radius );
    friend void UnsharpMask ( SingleSpectrum& spec, uint radius );

    //Unit conversions and other transformations that apply to the entire spectrum
    void InitialBin (uint bin_points = 72);
    void dBmToWatts();
    void WattsToExcessPower();
    void KSVZWeight();
    void LorentzianWeight();
    friend SingleSpectrum Spectrum::GrandSpectrum();
    friend SingleSpectrum Spectrum::Limits();
    friend SingleSpectrum Spectrum::GSquaredPrediction();

    double norm();
    double mean();
    double std_dev();

    double min_freq();
    double max_freq();

    uint size();

    std::string units();

    void rebin(uint points_per_bin);
    void chop_bins(uint start_chop, uint end_chop);
    double bin_start_freq(uint idx);
    double bin_mid_freq(uint idx);
    double bin_width();

    uint bin_at_frequency(double frequency);

  private:

    Units current_units = Units::dBm;

    void ParseRawData(std::string raw);
    uint num_lines(std:: string raw_data);

    void FillFromHeader(std::map<std::string, double> header);
    void PopulateUncertainties(uint rebin_size);

    double sum( std::vector<double>& data_list , double exponent = 1.0 );
    double mean( std::vector<double>& data_list );
    double std_dev( std::vector<double>& data_list );
    double kszv_power_per_bin( double freq_mhz );

    std::vector<double> sa_power_list;
    std::vector<double> uncertainties;

    double center_frequency = 0.0; //MHz
    double frequency_span = 0.0; //MHz
    double effective_volume = 0.0; // cm^3
    double noise_temperature = 0.0; // in Kelvin
    double Q = 0.0; //Quality Factor
    double b_field = 0.0; //Magnetic field in Tesla

    uint number_of_averages = 0; //Number of averages taken by instrument
    uint fft_points = 0; //Number of time-series points used to make FFT
};

#endif // SINGLESPECTRUM_H
