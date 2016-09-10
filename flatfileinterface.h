#ifndef FLATFILEINTERFACE_H
#define FLATFILEINTERFACE_H

//C++ System headers
#include <vector>//vector
#include <string>//string

//C System-Headers
//
//C++ System headers
#include <vector>
#include <string>
#include <map>
//Boost Headers
//
//Miscellaneous Headers
//

class FlatFileReader {

  public:
    FlatFileReader(std::string dir_name);
    ~FlatFileReader();

    std::string at(uint index);

    uint size();
    bool has(uint index);

  private:
    std::vector<std::string> raw_data_list;

    std::vector<std::string> EnumerateFiles(std::string dir_name, std::string sift_term);
    uint GetFileLines(std:: string file_name);

    std::string FastRead( std::string file_name);

};

class FlatFileParser : public FlatFileReader {

  public:
    FlatFileParser( std::string raw_data );
    ~FlatFileParser();

    std::vector<double> GetPowerList();
    std::map<std::string,double> GetHeader();

  private:

    void ParseRawData(std::string raw);
    uint NumLines(std:: string raw_data);

};

class FlatFileSaver {

  public:
    FlatFileSaver(std::string dir_name);
    ~FlatFileSaver();

    template <typename T>
    void load(std::map<std::string, T> header);

    template <typename T>
    void load(std::vector<T> vec);

    bool dump();

private:

    std::string cat();

    std::string save_file_path;

    std::map<std::string, std::string> header_map;
    std::vector<std::string> power_list;
};

#endif // FLATFILEINTERFACE_H
