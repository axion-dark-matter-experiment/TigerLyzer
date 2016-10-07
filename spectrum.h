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

enum class Units {dBm, Watts, ExcessPower, AxionPower};

class SingleSpectrum;

/*!
 * \brief Composite spectrum made up of individual spectra
 */
class Spectrum {
  public:
    Spectrum();
    ~Spectrum();

    Spectrum &operator+=(SingleSpectrum& spec);
    Spectrum &operator-=(SingleSpectrum& spec);

    /*!
     * \brief Use currently loaded spectra to build a Grand Spectra
     * \return Grand Spectra
     */
    SingleSpectrum GrandSpectrum();
    SingleSpectrum Limits();
    SingleSpectrum GSquaredPrediction();

    /*!
     * \brief Convert all currently loaded spectra from units of dBm to watts
     *
     * Will fail if loaded spectra are not in units of dBm
     */
    void dBmToWatts();
    /*!
     * \brief Convert all loaded spectra from units of Watts to units of
     * watts above power due to thermal noise.
     *
     * Will fail if all spectra are not in units of Watts
     */
    void WattsToExcessPower();
    /*!
     * \brief Convert all loaded spectra from units of excess power to
     * power added (or subtracted) due to axions in the cavity, as predicted
     * by KSVZ theory.
     *
     * Will fail if spectra are not in units of ExcessPower
     */
    void KSVZWeight();
    void LorentzianWeight();

    uint size();
    void clear();

    SingleSpectrum at(uint idx);

  private:

    SingleSpectrum BlankGrandSpectrum();

    double spectrum_weight(const SingleSpectrum& spec);
    std::vector<SingleSpectrum> spectra;

};

#endif // SPECTRUM_H
