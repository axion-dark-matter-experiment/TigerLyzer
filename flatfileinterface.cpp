//Header for this file
#include "flatfileinterface.h"

//C System-Headers
#include <termios.h>  /* POSIX terminal control definitions */
#include <sys/ioctl.h>
#include <fcntl.h>//fopen(),fclose()
#include <unistd.h>//read(), write()
#include <stdio.h>

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
#include <mutex>

//Boost Headers
#include <boost/algorithm/string.hpp>//split() and is_any_of for parsing .csv files
#include <boost/lexical_cast.hpp>//lexical cast (unsurprisingly)
#include <dirent.h>

//Miscellaneous Headers
#include <omp.h>//OpenMP pragmas

FlatFileReader::FlatFileReader(std::string dir_name) {

    std::vector<std::string> file_list = EnumerateFiles(dir_name, "SA_F");

    std::mutex guard;

    #pragma omp parallel for
    for ( uint i = 0 ; i < file_list.size() ; i ++) {

        guard.lock();
        raw_data_list.push_back( FastRead( file_list.at(i) ) );
        guard.unlock();
    }
}

FlatFileReader::~FlatFileReader() {
    raw_data_list.clear();
}

std::vector<std::string> FlatFileReader::EnumerateFiles(std::string dir_name, std::string sift_term) {

    DIR *dir;
    struct dirent *ent;
    const char* c_dir_name = dir_name.c_str();

    std::vector<std::string> file_names;

    if ((dir = opendir (c_dir_name)) != NULL) {

        while ((ent = readdir (dir)) != NULL) {

            std::string file_name = std::string (ent->d_name);

            if( file_name.find(sift_term) != std::string::npos ) {
                file_names.push_back( dir_name+file_name );
            }
        }

        return file_names;
        closedir (dir);
    } else {
        /* could not open directory */
        perror ("");
        return file_names;
    }
}

std::string FlatFileReader::FastRead( std::string file_name ) {
    const char* c_file_name = file_name.c_str();

    std::ifstream file_stream(c_file_name);
    std::stringstream buffer;
    buffer << file_stream.rdbuf();

    return buffer.str();
}

uint FlatFileReader::size() {
    return raw_data_list.size();
}

bool FlatFileReader::has(uint index) {
    return (index <= raw_data_list.size())?true:false;
}

std::string FlatFileReader::at(uint index) {
    if( index >= raw_data_list.size() ) {
        std::string err_mesg = "Requested index of ";
        err_mesg += boost::lexical_cast<std::string>(index);
        err_mesg +=" is greater than the number of loaded files ";
        err_mesg +="("+boost::lexical_cast<std::string>(raw_data_list.size() - 1)+")";
        throw std::out_of_range(err_mesg);
    }

    return raw_data_list.at(index);
}
