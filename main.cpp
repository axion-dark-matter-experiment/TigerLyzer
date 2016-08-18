#include "flatfileinterface.h"
#include "spectrum.h"

#include <iostream>
#include <mutex>

#include <chrono>

int main(int argc, char *argv[]) {

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0 ; i < 1 ; i++) {
        auto Reader = FlatFileReader("/home/bephillips2/workspace/Electric_Tiger_Control_Code/data/44_06_00_13.08.2016/");

        std::vector<SingleSpectrum> spectra;

        std::mutex guard;

//        #pragma omp parallel for
        for( uint j = 0 ; j < Reader.size() ; j ++) {

            guard.lock();
            std::cout << "Loading spectrum "<< j << std::endl;
            spectra.push_back( SingleSpectrum( Reader.at(j) ) );
            guard.unlock();
        }

        std::vector<double> dubs(spectra[0].size(), 115.0);

        auto new_spectra = spectra[0] - dubs;
        std::cout<< new_spectra;

    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = end - start;
    auto time_taken = fp_ms.count();

    std::cout<<"Took "<<time_taken<<" ms."<<std::endl;
}

