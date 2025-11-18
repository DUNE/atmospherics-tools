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
     std::shared_ptr<THnT<double>> hist = model.GetHistogram();

    // --- Pre-calculate strides for efficient global bin index calculation ---
    // The stride for an axis is the product of the number of bins of all preceding axes.
    std::vector<long long> strides(naxes);
    strides[0] = 1;
    for (uint i = 1; i < naxes; ++i) {
        strides[i] = strides[i - 1] * hist->GetAxis(i - 1)->GetNbins();
    }


     std::vector<ResultPtr_t> all_hists;

    ROOT::RDF::RNode df_with_bins = df.Define("global_bin_id", []() { return 0LL; }, {});

    // --- Iteratively build the global bin index ---
    for (uint i=0; i<naxes; ++i) {
        const auto& axis = binningAxes[i];
        // This lambda calculates the weighted contribution of a single axis to the global bin index.
        auto get_1d_bin_contribution_lambda = [hist_axis = hist->GetAxis(i), stride = strides[i]](const float& value, const long long& prev_id) -> long long {
            // return stride - 1;
            // FindBin returns a 1-based index. We need a 0-based index for the global formula.
            const int bin = hist_axis->FindBin(value);
            // Handle underflow/overflow, which FindBin places in bins 0 and nbins+1.
            return static_cast<long long>(bin) * stride + prev_id;
        };

        // For the first axis, create the "global_bin_id" column.
        df_with_bins = df_with_bins.Redefine("global_bin_id", get_1d_bin_contribution_lambda, {axis.variable, "global_bin_id"});
    }

    // auto df_with_bins = df.Define("global_bin_id", binning_lambda, binning_vars);

    for (const auto& syst : systematics) {
        MultiHistoNDHelper helper(model, syst.paramNodes.size());
        ResultPtr_t hists = df_with_bins.Book<WeightVec_t, long long>(std::move(helper), {syst.vectorWeightBranch, "global_bin_id"});
        all_hists.push_back(hists);
    }

    std::cout << "Booked all" << std::endl;
    ROOT::RDF::SaveGraph(df_with_bins, "./mydot.dot");

    for (auto mv_hists : all_hists) {
        for (const auto& hist : *mv_hists) {
            hist->Print("all");
        }
    }
    
}

///MAYBE CAN GO WITH SOME LOOP LOGIC TO ITERATIVELY BUILD THE COMPUTING GRAPH

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