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

class Spectrum {
  public:
    Spectrum();
    ~Spectrum();

    Spectrum &operator+=(SingleSpectrum& spec);
    Spectrum &operator-=(SingleSpectrum& spec);

    SingleSpectrum GrandSpectrum();
    SingleSpectrum Limits();
    SingleSpectrum GSquaredPrediction();

    void dBmToWatts();
    void WattsToExcessPower();
    void KSVZWeight();

    uint size();
    void clear();

    SingleSpectrum at(uint idx);

  private:

    SingleSpectrum BlankGrandSpectrum();

    double spectrum_weight(const SingleSpectrum& spec);
    std::vector<SingleSpectrum> spectra;

};

#endif // SPECTRUM_H
