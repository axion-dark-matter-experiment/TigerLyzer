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

/*!
 * \brief Object that handles basic file IO operations such as enumerating
 * files in a folder, opening files and loading file contents into strings.
 *
 * Upon initialization a FlatFileReader will search through a choosen directory
 *  pick out data files and load each file into a std::string. It should be noted
 *  that files are loaded from disk -in parallel-. Thus files are loaded into the
 *  classes without regard to ordered. Further parallel file IO will only result
 * in a performance increase on systems equipped with RAID or SSD's. Systems that
 * make use of non-RAID hard disks will see a signifigant performance -decrease-.
 */
class FlatFileReader {

  public:
    /*!
     * \brief Initialize a new Reader
     *
     * Upon construction the Reader will find all data files in a choosen
     * directory and load said files into individual strings.
     *
     *
     * \param dir_name
     * File path to the directory containing data collected by Electric Tiger
     *
     * \param sift_term
     * A short string used to filter data files from other files found in
     * dir_name. Only files containing sift_term in their names will be
     * loaded.
     *
     * Usually we expect data to be saved in the format "NAME_X.csv" for
     * NAME being some term common to all data files and X corresponding to an index.
     *
     * For example it may be the case that files are saved with names like
     * "SA_F0.csv", "SA_F1.csv" , "SA_F2.csv" etc, so an appropiate sift term
     * would be "SA_F" as this string is common to all data files.
     */
    FlatFileReader(std::string dir_name, std::string sift_term);
    ~FlatFileReader();

    /*!
     * \brief Return the raw file data at a certain index position.
     *
     * Function behaves identically to the at method of std::vector
     *
     * \param index
     * The position of the raw file data. Note that if index > size of class
     * this function will throw an out of bounds error.
     * \return
     * Raw data, as it was loaded from the file
     */
    std::string at(uint index);

    /*!
     * \brief Similar to std::vector::size(), get the number of data sets currently
     * loaded.
     *
     * \return Number of loaded data sets
     */
    uint size();

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
