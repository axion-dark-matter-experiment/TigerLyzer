#ifndef SPECTRUM_H
#define SPECTRUM_H

//C System-Headers
//
//C++ System headers
#include <vector>
#include <string>
#include <map>
#include <iostream>

//Boost Headers
//

//Miscellaneous Headers
//

//Project Specific Headers
//

enum class Units {dBm, Watts, ExcessPower, AxionPower, ExclLimit90};

class SingleSpectrum;

/*!
 * \brief Container class designed to hold all the individual spectra collected
 * in a data run.
 *
 * Class is designed to perform batch operations on multiple SingleSpectrum at once
 * e.g. calling unit conversion functions.
 *
 * Additionally class is designed to generated Grand Spectra and Exclusion Limits,
 * operations that require many individual spectra.
 */
class Spectrum {
  public:
    Spectrum();
    ~Spectrum();

    /*!
     * \brief Similar to std::vector::push_back()- insert a SingleSpectrum
     * at the back of the Spectrum class.
     *
     * \param spec
     * The SingleSpectrum class to be added.
     */
    Spectrum &operator+=(SingleSpectrum& spec);

    /*!
     * \brief Remove a SingleSpectrum class that has already been emplaced.
     *
     * \param spec
     * SingleSpectrum object to be remove- if no such object is present this
     *  function does nothing.
     */
    Spectrum &operator-=(SingleSpectrum& spec);

    /*!
     * \brief Combine all currently loaded spectra to form a Grand Spectrum.
     *
     * Power values and uncertainties for overlapping spectra are added
     *  using the rules for summing of normally distributed random variables,
     *
     * Note that the currently loaded spectra are not altered in any way
     *  when this function is called.
     *
     * \return
     * A Grand Spectrum with an appropiate min. and max. frequency.
     */
    SingleSpectrum GrandSpectrum();

    /*!
     * \brief Combine all currently
     * \return
     */
    SingleSpectrum Limits();
    SingleSpectrum GSquaredPrediction();

    /*!
     * \brief Call SingleSpectrum::dBmToWatts() on all loaded spectra.
     */

    void dBmToWatts();
    /*!
     * \brief Call SingleSpectrum::WattsToExcessPower() on all loaded spectra.
     */
    void WattsToExcessPower();

    /*!
     * \brief Call SingleSpectrum::KSVZWeight() on all loaded spectra.
     */
    void KSVZWeight();

    /*!
     * \brief Call SingleSpectrum::LorentzianWeight() on all loaded spectra.
     */
    void LorentzianWeight();

    /*!
     * \brief Similar to std::vector::size()- get the number of
     * elements in the Spectrum.
     *
     * \return
     * The total number of elements in the Spectrum
     */
    uint size();

    /*!
     * \brief Reset each power value and uncertainty of the Spectrum
     * to zero.
     */
    void clear();

    /*!
     * \brief Similiar to std::vector::at()- return the SingleSpectrum
     * at a particular index position.
     *
     * \param idx
     * \return
     */
    SingleSpectrum at(uint idx);

  private:

    SingleSpectrum BlankGrandSpectrum();

    double spectrum_weight(const SingleSpectrum& spec);
    std::vector<SingleSpectrum> spectra;

};

#endif // SPECTRUM_H
