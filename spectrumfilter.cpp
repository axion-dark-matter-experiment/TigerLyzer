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

//Note we are deliberately passing data_list by reference so that we can make a copy
//We can either choose to pad data_list itself with zeros, or copy its contents
//and pad the copy. Experiments suggest padding data_list is much faster than copying
//within the function body
std::vector<double> Convolve(std::vector<double> data_list,const std::vector<double>& kernel) {
    int n = data_list.size();
    int k = kernel.size();

    //Copy data_list into deque to speed up front/back insertions that will be needed
    //for zero padding
//    std::deque<double> data_deque (data_list.begin(),data_list.end());
    //Pad data_list with zeros at front and back

//    for (int i=1; i<=k; i++) {
//        data_deque.push_front(0.0);
//        data_deque.push_back(0.0);
//    }

    //reserver space for zero padding
    data_list.reserve ( n + 2*k );

    //Pad data_list with zeros at front and back
    //Yes, a deque would be faster
    //For now we are only focusing on vector operations
    //Try not to worry about the k*O(n) complexity operations
    //that we are incuring...
    for (int i = 0; i < k; i++) {
        data_list.insert(data_list.begin(), 0.0);
        data_list.push_back(0.0);
    }

    std::vector<double> convolved_list;
    convolved_list.reserve( n - k );

    #pragma omp parallel for ordered
    for(int i = 2*k; i< n+k; i++) {

        double conv_elem=0.0;

        for(int j = 0; j<k ; j++) {
            conv_elem += data_list.at(i-j)*kernel.at(j);
        }

        #pragma omp ordered
        convolved_list.push_back( conv_elem );
    }
    return convolved_list;
}

//Convolve the input list 'data_list' with a gaussian kernel with user defined radius
//serves as a low-pass filter that surpresses noise.
std::vector<double> GaussBlur(std::vector<double>& data_list, uint radius) {

    auto gauss_matrix = GaussKernel( radius );
    return Convolve(data_list,gauss_matrix);
}

std::vector<double> Unsharp(std::vector<double>& data_list, uint radius) {

    auto gauss_matrix = GaussKernel( radius );
    auto blurred_mat = Convolve( data_list, gauss_matrix );

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


void GaussianFilter( SingleSpectrum& spec, uint radius ) {
    spec.sa_power_list = GaussBlur(spec.sa_power_list, radius);
}

void UnsharpMask( SingleSpectrum& spec, uint radius ) {
    spec.sa_power_list = Unsharp(spec.sa_power_list, radius);
}
