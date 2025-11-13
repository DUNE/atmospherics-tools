#ifndef SPLINECALCULATOR_H
#define SPLINECALCULATOR_H

#include <ROOT/RDataFrame.hxx>
#include "TSpline.h"
#include "TH1D.h"
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <variant>

class SplineCalculator {
public:
    // Static method to enable multithreading - call this BEFORE creating any SplineCalculator
    static void EnableMultiThreading();
    
    SplineCalculator(const std::string& fileName, const std::string& treeName);

    void addBinningAxis(const std::string& variable, const std::vector<double>& edges);

    void addSystematic(const std::string& systName,
                       const std::vector<double>& systParamNodes,
                       const std::string& vectorWeightBranch,
                       int nominalIndex);

    void run();

    void writeSplines(const std::string& outputFileName) const;

    const std::map<std::string, std::map<int, std::unique_ptr<TSpline3>>>& getSplines() const;

private:
    // Template function to handle different numeric types
    template<typename T>
    auto createTypedBinningLambda(const std::vector<long long>& strides, 
                                  const std::vector<size_t>& axis_nbins) const;
    
    // Type detection helper
    std::string getColumnType(const std::string& columnName);
    
    struct Systematic {
        std::string name;
        std::vector<double> paramNodes;
        std::string vectorWeightBranch;
        int nominalIndex;
    };

    struct BinningAxis {
        std::string variable;
        std::vector<double> edges;
    };

    ROOT::RDataFrame df;
    std::vector<BinningAxis> binningAxes;
    std::vector<Systematic> systematics;
    std::map<std::string, std::map<int, std::unique_ptr<TSpline3>>> all_splines_map;
};

#endif // SPLINECALCULATOR_H