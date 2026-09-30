#include "Reader.h"
#include <argparse/argparse.hpp>
#include <filesystem>

int main(int argc, char const *argv[])
{
    argparse::ArgumentParser parser("weightor");

    parser.add_argument("-i", "--input")
        .required()
        .help("Input file directory.");
    try {
      parser.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
      std::cerr << err.what() << std::endl;
      std::cerr << parser;
      return 1;
    }

    std::string iFilePath(parser.get<std::string>("-i"));
    std::filesystem::path directory = iFilePath;
    
    double TotalPOT = 0;

    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
      if (entry.is_regular_file() && entry.path().extension()==".root") {
	std::cout << entry.path() << std::endl;	
	Reader<float> reader(entry.path(), "cafmaker");
	double POT = reader.POT();
	std::cout << "Got POT of: " << POT << std::endl;
	TotalPOT += POT;
      }
    }

    std::cout << "Total POT:" << TotalPOT << std::endl;
    return 0;
}
