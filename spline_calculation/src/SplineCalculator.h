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
#include "MultiHistoNDAction.hxx"

struct Systematic {
    std::string name;
    std::vector<double> paramNodes;
    std::string vectorWeightBranch;
    int nominalIndex;
};

/// Overload the << operator for std::ostream to allow easy printing of Systematic structs.
inline std::ostream& operator<<(std::ostream& os, const Systematic& syst) {
    os << "Systematic:\n"
       << "  Name: " << syst.name << "\n"
       << "  Vector Weight Branch: " << syst.vectorWeightBranch << "\n"
       << "  Nominal Index: " << syst.nominalIndex << "\n"
       << "  Parameter Nodes: [";
    for (size_t i = 0; i < syst.paramNodes.size(); ++i) {
        os << syst.paramNodes[i] << (i < syst.paramNodes.size() - 1 ? ", " : "");
    }
    os << "]";
    return os;
}



class SplineCalculator {
public:
    
    SplineCalculator(const std::string& fileName, const std::string& treeName);

    void addBinningAxis(const std::string& variable, const std::vector<double>& edges);

    void addSystematic(const std::string& systName,
                       const std::vector<double>& systParamNodes,
                       const std::string& vectorWeightBranch,
                       int nominalIndex);
    void addSystematic(const Systematic& syst);

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