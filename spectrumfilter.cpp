//Header for this file
#include "spectrumfilter.h"
//C System-Headers
#include <sys/ioctl.h>
#include <fcntl.h>//fopen(),fclose()
#include <unistd.h>//read(), write()
//C++ System headers
#include <vector>//vector
#include <string>//string
#include <fstream>//iss*
#include <chrono>// timing functions
#include <cmath>//sqrt, abs
#include <iostream>//cout
#include <typeinfo>//typeid
#include <algorithm> // transform, find
#include <functional> // plus/minus
#include <utility>//std::make_pair
#include <map>//std::map
//Boost Headers
#include <boost/algorithm/string.hpp>//split() and is_any_of for parsing .csv files
#include <boost/lexical_cast.hpp>//lexical cast (unsurprisingly)
//Miscellaneous Headers
#include <omp.h>//OpenMP pragmas
//Project Specific Headers
#include "singlespectrum.h"

//sum all enteries in a vector, with optional parameter of raising each entry to a power
//used by Standard Deviation function
inline double sum(std::vector<double>& data_list,double exponent) {
    double tot=0;
    for (auto& val : data_list) {
        tot+=pow(val,exponent);
    }
    return tot;
}

inline double norm(std::vector<double>& data_list ) {
    return sqrt(sum(data_list,2.0));
}

//define a gaussian function with standard deviation sigma and a mean value of zero
//used in the construction of gaussian kernels
inline double gaussian(double x, double sigma) {
    return 1/(sqrt(M_PI_2)*sigma)*exp( -0.5 *pow(x/sigma,2.0));
}

std::vector<double> Normalize(std::vector<double>& data_list) {
    double norm_factor=sqrt(sum(data_list,2));

    for(unsigned int i = 0; i<data_list.size(); i++) {
        data_list.at(i)=data_list.at(i)/norm_factor;
    }
    return data_list;
}

//generate a gaussian kernel of radius 'r', suitable for convolutions
//kernel will have a standard deviation of r/2.
std::vector<double> GaussKernel(int r) {

    double sigma = static_cast<double>(r)/2.0;
    std::vector<double> vals;
    for( int i = -r; i<= r ; i ++) {
        vals.push_back(gaussian(i,sigma));
    }

    //normalize kernel before returning
    return Normalize(vals);
}

template <typename T>
std::vector<T> LinearConvolve( std::vector<T>& signal, std::vector<T>& kernel) {

    int kernel_size = kernel.size();
    int half_k_size = (kernel_size - 1 )/2;

    int signal_size = signal.size();
    int signal_max_index = signal_size - 1;

    std::vector<T> output( signal_size , 0);

    #pragma omp parallel for
    for ( int i = 0 ; i < signal_size ; i++ ) {

        float conv_elem = 0.0f;

        int k_max = ( ( i + half_k_size ) > signal_max_index )?( signal_max_index + half_k_size - i ):(kernel_size);
        int k_min = ( ( i - half_k_size ) < 0 )?( half_k_size - i ):(0);

        for ( int j = k_min ; j < k_max ; j++ ) {

            conv_elem += signal[ i + j - half_k_size]*kernel[ j ];
        }

        for ( int j = 0 ; j < k_min ; j++ ) {

            conv_elem += signal[ -i - j + half_k_size ]*kernel[ j ];
        }

        for ( int j = k_max ; j < kernel_size ; j++ ) {

            conv_elem += signal[ 2*signal_max_index - i - j + half_k_size ]*kernel[ j ];

        }

        output[i] = conv_elem;

    }

    return output;
}

//Convolve the input list 'data_list' with a gaussian kernel with user defined radius
//serves as a low-pass filter that surpresses noise.
std::vector<double> GaussBlur(std::vector<double>& data_list, uint radius) {

    auto gauss_matrix = GaussKernel( radius );
    return LinearConvolve( data_list, gauss_matrix );
}

std::vector<double> Unsharp(std::vector<double>& data_list, uint radius) {

    auto gauss_matrix = GaussKernel( radius );
    auto blurred_mat = LinearConvolve( data_list, gauss_matrix );


    //Even though our Gaussian kernel was normalized we cannot expect the
    //norm of our convolved matrix to be equal to the norm of the original
    //matrix. We need to compensate for this by multiplying the convolved
    //matrix by the ratio- original_matrix_norm/convolved_matrix_norm;
    double norm_factor = (norm(data_list)/norm(blurred_mat));

    std::transform(blurred_mat.begin(),\
                   blurred_mat.end(),\
                   blurred_mat.begin(),\
                   std::bind1st(std::multiplies<double>(),\
                                norm_factor));

    data_list.resize(blurred_mat.size());

    std::transform(data_list.begin(),\
                   data_list.end(),\
                   blurred_mat.begin(),\
                   data_list.begin(),\
                   std::minus<double>());

    return data_list;
}

uint AutoOptimize(SingleSpectrum& spec) {

    double target = 1.0f/sqrt( spec.size() );
    std::vector<double> deltas;

    for ( uint i = 1; i < 35; i++) {
        auto copy = spec;

        UnsharpMask( copy, i );
        double aim = copy.mean()/copy.std_dev();

        double delta = target - aim;

        deltas.push_back( delta );

    }

    auto result = std::min_element(std::begin(deltas), std::end(deltas));

    auto i = std::distance(std::begin(deltas), result);
    return static_cast<uint>( i ) + 1;

}


void GaussianFilter( SingleSpectrum& spec, uint radius ) {
    spec.sa_power_list = GaussBlur(spec.sa_power_list, radius);
}

void UnsharpMask( SingleSpectrum& spec, uint radius ) {
    spec.sa_power_list = Unsharp(spec.sa_power_list, radius);
}
