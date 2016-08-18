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
#include <fstream>     //iss*
#include <chrono>      // timing functions
#include <cmath>       //sqrt, abs
#include <iostream>    //cout
#include <typeinfo>    //typeid
#include <algorithm>   // transform, find, count
#include <functional>  // plus/minus
#include <utility>     //std::make_pair
#include <map>         //std::map
#include <mutex> //protect against concurrent access when using (unordered) parallel for loops

// Boost Headers
#include <boost/algorithm/string.hpp>  //split() and is_any_of for parsing .csv files
#include <boost/lexical_cast.hpp>  //lexical cast (unsurprisingly)
#include <dirent.h>

// Miscellaneous Headers
#include <omp.h>  //OpenMP pragmas

SingleSpectrum::SingleSpectrum(std::string raw_data) {
    ParseRawData(raw_data);

    //    for (const auto& val: header) {
    //        std::cout<< val.first << ";" << val.second << std::endl;
    //    }
}

SingleSpectrum::~SingleSpectrum() {
    header.clear();
    sa_power_list.clear();
}

uint SingleSpectrum::size() {
    return sa_power_list.size();
}

SingleSpectrum &SingleSpectrum::operator*=(double scalar) {

    for(unsigned int i=0; i<sa_power_list.size(); i++) {
        sa_power_list[i]*=scalar;
//        uncertainty[i]*=d;
    }
    return *this;
}

SingleSpectrum &SingleSpectrum::operator+=(double scalar) {

    for(unsigned int i=0; i<sa_power_list.size(); i++) {
        sa_power_list[i]+=scalar;
    }
    return *this;
}

std::ostream& operator << (std::ostream& stream, SingleSpectrum& spectrum) {

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

    for(unsigned int i=0; i<spectra_a.size(); i++) {
        spectra_c.sa_power_list[i] =\
                                    spectra_a.sa_power_list[i] +\
                                    spectra_b.sa_power_list[i];
    }

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

    for(unsigned int i=0; i<spectra_a.size(); i++) {
        spectra_c.sa_power_list[i] =\
                                    spectra_a.sa_power_list[i] +\
                                    spectra_b[i];
    }

    return spectra_c;
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

    for(unsigned int i=0; i<spectra_a.size(); i++) {
        spectra_c.sa_power_list[i] =\
                                    spectra_a.sa_power_list[i] -\
                                    spectra_b.sa_power_list[i];
    }

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

    for(unsigned int i=0; i<spectra_a.size(); i++) {
        spectra_c.sa_power_list[i] =\
                                    spectra_a.sa_power_list[i] -\
                                    spectra_b[i];
    }

    return spectra_c;
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
        std::string err_mesg = "Spectra and vector are not the same size ";
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

uint SingleSpectrum::NumLines(std::string raw_data) {
    return std::count(raw_data.begin(), raw_data.end(), '\n');
}

double inline line_to_val(std::istringstream& stream) {
    std::string input;
    std::getline(stream, input);
    return boost::lexical_cast<double>(input);
}

void SingleSpectrum::ParseRawData(std::string raw_data) {
    uint lines = NumLines(raw_data);

    std::istringstream data_stream(raw_data);

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

    for (uint i = 12; i < lines; i++) {
        std::string input;
        std::getline(data_stream, input);
        double val = boost::lexical_cast<double>(input);

        sa_power_list.push_back(val);
    }
}
