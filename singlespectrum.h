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
 * Class is designed to 'look and feel' like a \f$ \mathcal{R}^n \f$ vector.
 * Spectra allows scalar multiplication and addition, along with 'vector'
 * multiplication and addition.
 *
 * Further each SingleSpectrum tracks its own uncertainties and current units.
 */
class SingleSpectrum {

  public:

    /*!
     * \brief Build a new SingleSpectrum from a string of raw data collected from the experiment.
     *
     * Electric Tiger collects one data set per cavity length and saves the results as a plain text
     *  (.csv) file. The first few lines of each data file are the "header" and contain experiment parameters
     *  such as cavity length, peak center frequency, Quality Factor etc. Each entry in the header file is
     * has the form "parameter;value".
     *
     * At the end of the header file is a token, usually "@". Beyond this token
     * is the collected power spectrum. Each value in the spectrum is a floating point value followed by a newline.
     * Below is a sample data file:
     *
     * \f{verbatim}{
     *   sa_span;10
     *   fft_length;131072
     *   effective_volume;0.1
     *   bfield;1.54
     *   noise_temperature;400
     *   sa_averages;256
     *   Q;128.4913492063492
     *   actual_center_freq;4037.3840399002493
     *   fitted_hwhm;15.71072319201995
     *   cavity_length;7.693
     *   @
     *   -116.0639877
     *   -116.0132904
     *   -115.6643753
     *   -116.1841431
     *   -115.1884918
     *   etc...
     * \f}
     *
     * It should be noted that the header string must contain particular enteries on the constructor
     * will throw an exception. In particular the keys
     * \f{verbatim}{
        "sa_span",
        "fft_length",
        "effective_volume",
        "bfield",
        "noise_temperature",
        "sa_averages",
        "Q",
        "actual_center_freq",
        "fitted_hwhm",
        "bfield"
     * \f}
     *
     * Must all be present.
     *
     *
     * \throws std::invalid_argument
     * Thrown if the header is missing certain keys or is formatted incorrectly.
     *
     * \param raw_data
     * A std::string containing the -entire- contents of a data file.
     */
    SingleSpectrum(std::string raw_data);
    /*!
     * \brief Construct a blank ( all power values and uncertainties = 0 ) SingleSpectrum
     * with a particular number of enteries.
     *
     * Note that blank spectra do will start with all physical parameters, such as min. and max.
     *  frequency set to zero.
     *
     * \param size
     * The number of enteries that the blank spectrum should have.
     */
    SingleSpectrum(uint size);

    /*!
     * \brief Construct a blank ( all power values and uncertainties = 0 ) SingleSpectrum
     * with a particular number of enteries and a particular mim. and max. frequency.
     *
     * \param size
     * The number of enteries that the blank spectrum should have.
     *
     * \param min_freq
     * The minimum frequency of the spectrum in MHz.
     *
     * \param max_freq
     * The maximum frequency of the spectrum in MHz.
     */
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

    friend SingleSpectrum Spectrum::GrandSpectrum();
    friend SingleSpectrum Spectrum::Limits();
    friend SingleSpectrum Spectrum::GSquaredPrediction();

    /*!
     * \brief Perform initial binning of a raw power spectrum and initializes spectrum uncertainties.
     *
     * \param bin_points
     * The number of points that should be put into a single bin.
     */
    void InitialBin (uint bin_points = 32);

    /*!
     * \brief Covert from units of dBm to Watts.
     *
     * \throws std::invalid_argument
     *Thrown if spectrum is not in units of dBm.
     */
    void dBmToWatts();

    /*!
     * \brief Convert from units of Watts to ExcessPower, i.e.
     * Watts above power due to noise.
     *
     * \throws std::invalid_argument
     * Thrown if spectrum is not in units of Watts
     */
    void WattsToExcessPower();

    /*!
     * \brief Convert spectrum from units of Excess power to units of \f$ g_{\gamma\gamma} \f$
     *
     * See Ed Daw's Thesis, pg. 113, eq. 5.3
     *
     * \throws std::invalid_argument
     * Thrown if spectrum is not in units of ExcessPower.
     */
    void KSVZWeight();

    /*!
     * \brief Weight each point in the spectrum by how far away it is from the center frequency.
     *
     * Points are weighted as follows; Let \f$ f_0 \f$ be the center frequency of the spectrum,
     *  \f$ f(n) \f$ be the center frequency of the n'th bin (for \f$ 0 \leq n \leq \text{size} \f$,
     * and Q be the quality factor. Then \f$ \forall i \ \frac{P_i}{L(f,f_0,Q)} \f$ Where
     * \f$ L(f,f_0,Q) \f$ represents the Lorentzian function.
     */
    void LorentzianWeight();

    /*!
     * \brief Compute the \f$ L_2 \f$ norm of the current spectrum =
     * \f$ \sqrt{\sum_{i=1}^n | P_i |^2} \f$ Then
     */
    double norm();
    /*!
     * \brief Compute the arithmetic mean of the current spectrum
     */
    double mean();
    /*!
     * \brief Compute the standard deviation of the current spectrum
     */
    double std_dev();

    /*!
     * \brief Get the minimum (smallest) frequency stored in the current spectrum
     */
    double min_freq();
    /*!
     * \brief Get the maximum (largest) frequency stored in the current spectrum
     */
    double max_freq();

    /*!
     * \brief similar to std::vector::size(), get the number of points in the current spectrum.
     *
     * Note that this corresponds to the number of power spectrum points but does not count
     * the number of uncertainty values stored.
     */
    uint size();

    /*!
     * \brief Get the current units that the Spectrum is in.
     *
     * Possible units are:
     *
     * dBm\n
     * Watts\n
     * ExcessPower = Power above noise (in Watts)\n
     * AxionPower = Power deposited or subtracted by axion signal in units of fraction of
     * expected power as predicted by KSVZ theory\n
     * ExclLimit90 = Not a true unit of power, 90% confidence limit on \f$ g_{a\gamma\gamma} ( GeV^{-1} ) \f$
     *
     * \return
     * string expressing the current units the spectrum is in
     */
    std::string units();
    /*!
     * \brief Rebin the current spectra, averaging all power spectrum points and uncertainties
     *
     * \param points_per_bin
     * The number of points that should be averaged into a signal bin
     */
    void rebin(uint points_per_bin);
    /*!
     * \brief Elminate bins from the head and tail of the current spectrum.
     *
     * \param start_chop
     * All points with indices 0 \f$ \leq \f$ start_chop will be removed.
     *
     * \param end_chop
     * All points with indices end_chop \f$ \leq \f$ size() will be removed.
     */
    void chop_bins(uint start_chop, uint end_chop);

    /*!
     * \brief Get the frequency corresponding to the left-hand side of a bin
     * \param idx
     * The index of the requested bin
     * \return
     */
    double bin_start_freq(uint idx);
    /*!
     * \brief  Get the frequency corresponding to the center of a bin
     * \param idx
     * The index of the requested bin
     * \return
     */
    double bin_mid_freq(uint idx);
    /*!
     * \brief Get the frequency corresponding to the right-hand side of a bin
     * \param idx
     * The index of the requested bin
     * \return
     */
    double bin_width();

    /*!
     * \brief Get the index of the bin correspond to a given frequency
     *
     * \param frequency
     * \return
     *
     * \throws std::out_of_range
     * Thrown if the requested frequency is not in the current spectrum, that is
     * if frequency \f$ \notin \f$ [ min_freq(), max_freq() ]
     */
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
