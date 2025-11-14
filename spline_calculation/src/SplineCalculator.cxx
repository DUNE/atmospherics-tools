#include "SplineCalculator.h"
#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "TFile.h"
#include "ROOT/RDFHelpers.hxx"


std::string SplineCalculator::getColumnType(const std::string& columnName) {
    // Use RDataFrame's GetColumnType method to detect actual type
    std::string rawType = df.GetColumnType(columnName);
    
    // Print debug information
    std::cout << "Column '" << columnName << "' raw type: '" << rawType << "'" << std::endl;
    
    return rawType;
}

using WeightVec_t = ROOT::RVec<double>;
using ResultPtr_t = ROOT::RDF::RResultPtr<std::vector<std::shared_ptr<THnT<double>>>>;

template<typename T>
auto SplineCalculator::createTypedBinningLambda(const std::vector<long long>& strides, 
                                               const std::vector<size_t>& axis_nbins) const {
    return [axes = this->binningAxes, strides, axis_nbins](const ROOT::RVec<T>& vals) -> long long {
        long long global_idx = 0;
        
        for (size_t i = 0; i < axes.size(); ++i) {
            const auto& axis = axes[i];
            const double val = static_cast<double>(vals[i]); // Only convert when needed for comparison

            if (val < axis.edges.front() || val >= axis.edges.back()) {
                return -1LL;
            }

            auto it = std::lower_bound(axis.edges.begin(), axis.edges.end(), val);
            int local_idx = std::distance(axis.edges.begin(), it) - 1;
            global_idx += local_idx * strides[i];
        }
        return global_idx;
    };
}

SplineCalculator::SplineCalculator(const std::string& fileName, const std::string& treeName)
    : df(treeName, fileName) {
    auto colNames = df.GetColumnNames();

    // 3. Loop over the list and print each name
    std::cout << "Column Names:" << std::endl;
    for (const auto& name : colNames) {
        std::cout << "- " << name << std::endl;
    }
    // Don't enable MT after RDataFrame construction - this causes threading issues
    // ROOT::EnableImplicitMT();
}

void SplineCalculator::addBinningAxis(const std::string& variable, const std::vector<double>& edges) {
    binningAxes.push_back({variable, edges});
}

void SplineCalculator::addSystematic(const std::string& systName,
                                     const std::vector<double>& systParamNodes,
                                     const std::string& vectorWeightBranch,
                                     int nominalIndex) {
    addSystematic({systName, systParamNodes, vectorWeightBranch, nominalIndex});
}

void SplineCalculator::addSystematic(const Systematic& syst) {
    systematics.push_back(syst);
}

void SplineCalculator::run() {
    if (binningAxes.empty()) {
        throw std::runtime_error("Binning must be set before running.");
    }
    if (systematics.empty()) {
        throw std::runtime_error("At least one systematic must be added before running.");
    }

    uint naxes = binningAxes.size();

    // Collect all binning variable names
    std::vector<std::string> binning_vars;
    binning_vars.reserve(naxes);

    std::vector<std::vector<double>> binning_edges;
    binning_edges.reserve(naxes);

    std::vector<int> binning_nbins;
    binning_nbins.reserve(naxes);



    for (const auto& axis : binningAxes) {
        binning_vars.push_back(axis.variable);
        binning_edges.push_back(axis.edges);
        binning_nbins.push_back(axis.edges.size() - 1);
    }

     ROOT::RDF::THnDModel model("", "", naxes, binning_nbins, binning_edges);
     std::vector<ResultPtr_t> all_hists;

    for (const auto& syst : systematics) {
        MultiHistoNDHelper helper(model, syst.paramNodes.size(), binning_vars);
        std::vector<std::string> col_list = {syst.vectorWeightBranch};
        col_list.insert(col_list.end(), binning_vars.begin(), binning_vars.end());

        ResultPtr_t hists = df.Book<WeightVec_t, float, float, float>(std::move(helper), col_list);
        all_hists.push_back(hists);
    }

    std::cout << "Booked all" << std::endl;
    // ROOT::RDF::SaveGraph(df, "./mydot.dot");

    for (auto mv_hists : all_hists) {
        for (const auto& hist : *mv_hists) {
            hist->Print("all");
        }
    }
    
}

void SplineCalculator::writeSplines(const std::string& outputFileName) const {
    TFile outFile(outputFileName.c_str(), "RECREATE");
    if (!outFile.IsOpen()) {
        throw std::runtime_error("Could not open output file: " + outputFileName);
    }

    for (const auto& [syst_name, splines_map] : all_splines_map) {
        if (!outFile.Get(syst_name.c_str())) {
            outFile.mkdir(syst_name.c_str());
        }
        outFile.cd(syst_name.c_str());
        for (const auto& [bin_idx, spline] : splines_map) {
            spline->Write();
        }
        outFile.cd("..");
    }
    outFile.Close();
}

const std::map<std::string, std::map<int, std::unique_ptr<TSpline3>>>& SplineCalculator::getSplines() const {
    return all_splines_map;
}