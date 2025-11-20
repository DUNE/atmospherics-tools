#ifndef SPLINECALCULATOR_H
#define SPLINECALCULATOR_H

#include <ROOT/RDataFrame.hxx>
#include "TSpline.h"
#include "TGraph.h"
#include "TH1D.h"
#include <vector>
#include <string>
#include <map>
#include <memory>
#include <variant>
#include "MultiHistoNDAction.hxx"
#include "SplineContainer.h"

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

enum class InterpolationType {
    Spline, // TSpline3
    Linear  // TGraph
};

class SplineCalculator {
public:
    
    SplineCalculator(const std::string& fileName, const std::string& treeName);
    SplineCalculator(const std::string& fileName, const std::string& treeName, const std::vector<std::string>& friendFileNames, const std::vector<std::string>& friendTreeNames);

    void setInterpolationType(InterpolationType type);

    void addBinningAxis(const std::string& variable, const std::vector<double>& edges);
    void addCategoricalBinningAxis(const std::string& variable, const std::vector<int>& categories);

    void addSystematic(const std::string& systName,
                       const std::vector<double>& systParamNodes,
                       const std::string& vectorWeightBranch,
                       int nominalIndex);
    void addSystematic(const Systematic& syst);

    template<typename F>
    void addSelection(F&& f, const std::vector<std::string>& columns) {
        df = df.Filter(std::forward<F>(f), columns, "User-defined functor/lambda selection");
    }

    /**
     * @brief Adds a new column from an external std::vector by copying its data.
     * 
     * The vector must have the same number of elements as the RDataFrame has entries.
     * Note: This version creates a copy of the vector to manage its lifetime safely.
     * For large vectors, consider using the overload that accepts a std::shared_ptr to avoid the copy.
     * 
     * @tparam T The type of the elements in the vector.
     * @param columnName The name of the new column to be created.
     * @param dataVector The std::vector containing the data for the new column.
     */
    template<typename T>
    void addColumnFromVector(const std::string& columnName, const std::vector<T>& dataVector) {        
        auto dataPtr = std::make_shared<std::vector<T>>(dataVector);
        addColumnFromVector(columnName, dataPtr);
    }

    /**
     * @brief Adds a new column from an external std::vector via a shared_ptr (most memory-efficient).
     * 
     * The vector must have the same number of elements as the RDataFrame has entries.
     * This method avoids copying the vector data, making it ideal for large datasets.
     * The RDataFrame computation graph will share ownership of the vector data.
     * 
     * @tparam T The type of the elements in the vector.
     * @param columnName The name of the new column to be created.
     * @param dataVectorPtr A std::shared_ptr to the vector containing the data.
     */
    template<typename T>
    void addColumnFromVector(const std::string& columnName, std::shared_ptr<std::vector<T>> dataVectorPtr) {
        auto nEntries = df.Count().GetValue();
        if (dataVectorPtr->size() != nEntries) {
            throw std::runtime_error("The size of the provided vector (" + std::to_string(dataVectorPtr->size()) +
                                     ") does not match the number of entries in the RDataFrame (" + std::to_string(nEntries) + ").");
        }

        if constexpr (std::is_same_v<T, bool>) {
            // Special handling for std::vector<bool> to avoid issues with its proxy reference type.
            // We must explicitly cast to bool to ensure RDataFrame gets a default-constructible type.
            df = df.Define(columnName, [dataVectorPtr](ULong64_t entry) { return static_cast<bool>((*dataVectorPtr)[entry]); }, {"rdfentry_"});
        } else {
            // Default implementation for all other vector types.
            df = df.Define(columnName, [dataVectorPtr](ULong64_t entry) { return (*dataVectorPtr)[entry]; }, {"rdfentry_"});
        }
    }

    void setEventWeightColumn(const std::string& columnName);

    ULong64_t getNEntries();

    void run();

    void writeSplines(const std::string& outputFileName) const;

    const std::map<std::string, std::map<int, SplineVariant>>& getSplines() const;
    std::unique_ptr<SplineContainer> getSplineContainer(const std::string& systName) const;

private:
    // Template function to handle different numeric types
    template<typename T>
    auto createTypedBinningLambda(const std::vector<long long>& strides, 
                                  const std::vector<size_t>& axis_nbins) const;
    
    // Type detection helper
    std::string getColumnType(const std::string& columnName);

    static ROOT::RDataFrame createDataFrameWithFriends(const std::string& fileName, const std::string& treeName,
                                                       const std::vector<std::string>& friendFileNames, const std::vector<std::string>& friendTreeNames);

    struct BinningAxis {
        std::string variable;
        std::vector<double> edges;
    };
    std::map<std::string, std::map<int, float>> categorical_maps;

    InterpolationType interpolationType = InterpolationType::Spline;
    ROOT::RDF::RNode df;
    std::vector<BinningAxis> binningAxes;
    std::string eventWeightColumn;
    std::vector<Systematic> systematics;
    std::map<std::string, std::map<int, SplineVariant>> all_splines_map;
};

#endif // SPLINECALCULATOR_H