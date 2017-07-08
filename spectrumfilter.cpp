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
#include <limits>// std::numeric_limits<double>::max()
#include <utility>// std::pair
//Boost Headers
#include <boost/algorithm/string.hpp>//split() and is_any_of for parsing .csv files
#include <boost/lexical_cast.hpp>//lexical cast (unsurprisingly)
//Miscellaneous Headers
#include <omp.h>//OpenMP pragmas
//Project Specific Headers
#include "singlespectrum.h"

//sum all enteries in a vector, with optional parameter of raising each entry to a power
//used by Standard Deviation function
double sum(std::vector<double>& data_list,double exponent) {

    double tot = 0;
    for ( auto& val : data_list ) {
        tot+= std::pow( val, exponent );
    }

    return tot;
}

inline double norm(std::vector<double>& data_list ) {
    return std::sqrt(sum(data_list,2.0));
}

//define a gaussian function with standard deviation sigma and a mean value of zero
//used in the construction of gaussian kernels
inline double gaussian(double x, double sigma) {
    return 1.0/(std::sqrt(M_PI_2)*sigma)*std::exp( -0.5 *std::pow(x/sigma,2.0));
}

std::vector<double> Normalize(std::vector<double>& data_list) {
    double norm_factor=std::sqrt(sum(data_list,2));

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
    return Normalize( vals );
}

//generate a gaussian kernel of radius 'r', suitable for convolutions
//kernel will have a standard deviation set by user
std::vector<double> GaussKernel( int r, double sigma ) {

    std::vector<double> vals;
    for( int i = -r; i<= r ; i ++) {
        vals.push_back(gaussian(i,sigma));
    }

    //normalize kernel before returning
    return Normalize( vals );
}


std::vector<double> UnsharpKernel( int radius, double sigma ) {

    std::vector<double> kernel( 2*radius + 1 );

    for( int i = -radius; i <= -1 ; i ++) {
        kernel.push_back( -1.0*gaussian( i, sigma ) );
    }

    kernel.push_back( 1.0 - gaussian( 0, sigma ) );

    for( int i = 1; i <= radius ; i ++) {
        kernel.push_back( -1.0*gaussian( i, sigma ) );
    }

    //normalize kernel before returning
    return Normalize( kernel );
}

inline double special_sinc( const double x, double f_t ) {
    if( x == 0.0 ) {
        return 2.0*f_t;
    } else {
        return std::sin(2*M_PI*f_t*x)/(M_PI*x);
    }
}

std::vector< double > sinc_kernel( int radius, double cutoff_frequency, double sample_frequency ) {

    std::vector<double> vals( 2*radius + 1 );
    double f_t = cutoff_frequency/sample_frequency;

    for( int i = -radius; i <= radius ; i ++) {

        double i_f = static_cast<double>(i);
        vals.push_back( special_sinc( i_f, f_t ) );
    }

    //normalize kernel before returning
    return Normalize( vals );

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

std::vector<double> Unsharp( std::vector<double>& data_list, uint radius, double sigma ) {

    auto gauss_matrix = GaussKernel( radius, sigma );
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

std::vector<double> SincFilter( std::vector<double>& data_list, uint radius, double cutoff_frequency, double sample_frequency ) {

    auto sinc_matrix = sinc_kernel( radius, cutoff_frequency, sample_frequency );
    auto sharpened_signal = LinearConvolve( data_list, sinc_matrix );


    //Even though our Gaussian kernel was normalized we cannot expect the
    //norm of our convolved matrix to be equal to the norm of the original
    //matrix. We need to compensate for this by multiplying the convolved
    //matrix by the ratio- original_matrix_norm/convolved_matrix_norm;
    double norm_factor = ( norm(data_list)/norm(sharpened_signal) );

    std::transform(sharpened_signal.begin(),\
                   sharpened_signal.end(),\
                   sharpened_signal.begin(),\
                   std::bind1st(std::multiplies<double>(),\
                                norm_factor));

    return sharpened_signal;
}

//std::vector<double> Unsharp( std::vector<double>& data_list, uint radius, double sigma ) {

//    auto unsharp_matrix = UnsharpKernel( radius, sigma );
//    auto sharpened_signal = LinearConvolve( data_list, unsharp_matrix );


//    //Even though our Gaussian kernel was normalized we cannot expect the
//    //norm of our convolved matrix to be equal to the norm of the original
//    //matrix. We need to compensate for this by multiplying the convolved
//    //matrix by the ratio- original_matrix_norm/convolved_matrix_norm;
//    double norm_factor = ( norm(data_list)/norm(sharpened_signal) );

//    std::transform(sharpened_signal.begin(),\
//                   sharpened_signal.end(),\
//                   sharpened_signal.begin(),\
//                   std::bind1st(std::multiplies<double>(),\
//                                norm_factor));

//    return sharpened_signal;
//}

double mean( std::vector< double >& data_list ) {
    //compute mean value of data set
    double sum_x = sum( data_list, 1.0 );
    double n = static_cast<double>( data_list.size() );
    return sum_x/n;
}

double std_dev ( std::vector<double> &data_list ) {

    //compute mean value of data set
    double sum_x = sum( data_list,1.0 );
    double n = static_cast<double>( data_list.size() );
    double mean = sum_x/n;

    //compute variance taking into account Bessel's correction i.e. n/(n-1)
    double sum_x2 = sum( data_list, 2.0 );
    double sigma_sqr = sum_x2/( n - 1.0 )-n/( n-1.0 )*std::pow( mean, 2.0 );

    //return square root of variance
    return std::sqrt(sigma_sqr);
}

std::pair< uint, double > AutoOptimize( SingleSpectrum& spec, uint max_radius, double sample_frequency ) {

    double target = 1.0/sqrt( static_cast<double>( spec.size() ) );
    std::vector<double> spec_data = spec.sa_power_list;
    Normalize( spec_data );

    double smallest_delta = std::numeric_limits<double>::max();
    std::pair< uint, double > optimal_parameters { 0, 0.0 };

    double max_frequency = sample_frequency/2.0;

    for ( uint radius_it = 1; radius_it < max_radius; radius_it++ ) {

        for( uint sigma_it = 0; sigma_it < max_frequency; sigma_it++ ) {

            double sigma_f = static_cast<double>( sigma_it ) + 0.01;

            std::cout << "Trying parameters: ("<< radius_it << "," << sigma_f <<")" << std::endl;

            auto test_signal = Unsharp( spec_data, radius_it, sigma_f );

            double aim = mean( test_signal )/std_dev( test_signal );

            double delta = std::abs( target - aim );

            std::cout << "Ratio of delta to target: " << aim/target << std::endl;

            if( delta < smallest_delta ) {

                smallest_delta = delta;
                optimal_parameters = std::pair< uint, double > { radius_it, sigma_f };
                std::cout << "New parameters were smaller!" << std::endl;
            }

        }

    }

    return optimal_parameters;

}


//std::pair< uint, double > AutoOptimize( SingleSpectrum& spec, uint max_radius, uint max_sigma ) {

//    double target = 1.0/sqrt( static_cast<double>( spec.size() ) );
//    std::vector<double> spec_data = spec.sa_power_list;
//    Normalize( spec_data );

//    double smallest_delta = std::numeric_limits<double>::max();
//    std::pair< uint, double > optimal_parameters { 0, 0.0 };

//    for ( uint radius_it = 1; radius_it < max_radius; radius_it++ ) {

//        for( uint sigma_it = 0; sigma_it < max_sigma*10; sigma_it++ ) {

//            double sigma_f = static_cast<double>( sigma_it )/10.0 + 0.01;

//            std::cout << "Trying parameters: ("<< radius_it << "," << sigma_f <<")" << std::endl;

//            auto test_signal = Unsharp( spec_data, radius_it, sigma_f );

//            double aim = mean( test_signal )/std_dev( test_signal );

//            double delta = std::abs( target - aim );

//            std::cout << "Ratio of delta to target: " << aim/target << std::endl;

//            if( delta < smallest_delta ) {

//                smallest_delta = delta;
//                optimal_parameters = std::pair< uint, double > { radius_it, sigma_f };
//                std::cout << "New parameters were smaller!" << std::endl;
//            }

//        }

//    }

//    return optimal_parameters;

//}


void GaussianFilter( SingleSpectrum& spec, uint radius ) {
    spec.sa_power_list = GaussBlur(spec.sa_power_list, radius);
}

void UnsharpMask( SingleSpectrum& spec, uint radius, double sigma ) {
    spec.sa_power_list = Unsharp( spec.sa_power_list, radius, sigma );
}
