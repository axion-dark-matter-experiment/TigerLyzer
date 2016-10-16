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

FlatFileReader::FlatFileReader(std::string dir_name, std::string sift_term ) {

    std::vector<std::string> file_list = EnumerateFiles(dir_name, sift_term);

    std::mutex guard;

    //Yes, we are reading from disk in parallel
    //It is important to note if this approach is used on a non-RAID
    //hard disk there will be a signifigant performance DECREASE.

    #pragma omp parallel for
    for ( uint i = 0 ; i < file_list.size() ; i ++) {

        std::lock_guard<std::mutex> lock (guard);
        raw_data_list.push_back( FastRead( file_list.at(i) ) );
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
        std::string err_mesg = __FUNCTION__;
        err_mesg += ": Could not open directory.";
        throw std::invalid_argument(err_mesg);
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

FlatFileSaver::FlatFileSaver(std::string dir_name) {
    save_file_path = dir_name;
}

template <typename T>
void FlatFileSaver::load(std::vector<T> vec) {

    for (const auto& val : vec ) {
        try {
            std::string num = boost::lexical_cast<std::string>(val);
            power_list.push_back(num);
        } catch (const boost::bad_lexical_cast& err) {
            std::cout << "Could not add element! Reason:" << err.what() << std::endl;
            return;
        }
    }
}

template <typename T>
void FlatFileSaver::load(std::map<std::string, T> header) {

    for (const auto& key_val : header ) {
        try {
            std::string val = boost::lexical_cast<std::string>(key_val.second);
            header_map [key_val.first] = val;
        } catch (const boost::bad_lexical_cast& err) {
            std::cout << "Could not add element! Reason:" << err.what() << std::endl;
            return;
        }
    }
}

std::string FlatFileSaver::cat() {
    std::string total;

    for ( const auto& val : power_list ) {
        total += val + "\n";
    }

    return total;
}

bool FlatFileSaver::dump() {
    const char* c_path=save_file_path.c_str();

    if ( power_list.size() <= 1 ) {
        std::cout << "Nothing to write to disk." << std::endl;
        return false;
    }
    const char* output = cat().c_str();	//convert from std::string to character array

    FILE* config;
    config = std::fopen(c_path, "w");//open the text file specified by config_file_name and write data to file

    if (config != NULL) {
        std::fputs(output, config);
        std::fclose(config);	//close config file
        return true;
    } else {
        //if config if not opened properly, exit function
        std::cout<<"Failed to write to file"<<std::endl;
        return false;
    }
}
