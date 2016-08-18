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

class SingleSpectrum {
  public:

    SingleSpectrum(std::string raw_data);
    ~SingleSpectrum();

    SingleSpectrum &operator*=(double scalar);
    SingleSpectrum &operator+=(double scalar);

    friend std::ostream& operator<< (std::ostream& stream, SingleSpectrum& spectrum);

    friend SingleSpectrum operator* (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b);
    friend SingleSpectrum operator+ (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b);
    friend SingleSpectrum operator- (SingleSpectrum& spectra_a, SingleSpectrum& spectra_b);

    friend SingleSpectrum operator* (SingleSpectrum& spectra_a, std::vector<double>& spectra_b);
    friend SingleSpectrum operator+ (SingleSpectrum& spectra_a, std::vector<double>& spectra_b);
    friend SingleSpectrum operator- (SingleSpectrum& spectra_a, std::vector<double>& spectra_b);


    void Normalize();

    double std_dev();
    double mean();
    uint size();

  private:
    void ParseRawData(std::string raw);
    uint NumLines(std:: string raw_data);

    double StdDev(std::vector<double> &data_list);
    std::vector<double> Normalize(std::vector<double> &data_list);

    std::vector<double> sa_power_list;
    std::map<std::string,double> header;

};

class Spectrum {
  public:
    Spectrum();
    ~Spectrum();

    void push_back(SingleSpectrum spectrum);
    void clear();

  private:
    std::vector<SingleSpectrum> spectra;

};

#endif // SPECTRUM_H
