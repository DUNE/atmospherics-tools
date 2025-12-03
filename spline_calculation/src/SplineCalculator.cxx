#include "SplineCalculator.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "TFile.h"
#include "ROOT/RDFHelpers.hxx"
#include "TAxis.h"

//TODO: Nicely deal with category events -> implemented, to be checked
//TODO: Allow for on the fly selection -> implemented, to be checked.
//TODO: Allow for external selections (vectors) -> implemented, to be checked
//TODO: Actually build the Splines -> implemented, to be checked
//TODO: The number of valid splines should be the same for all systematics -> BUG


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
    : df(ROOT::RDataFrame(treeName, fileName)) {
    // df is now an RNode, initialized from an RDataFrame
    auto colNames = df.GetColumnNames();

    // 3. Loop over the list and print each name
    std::cout << "Column Names:" << std::endl;
    for (const auto& name : colNames) {
        std::cout << "- " << name << std::endl;
    }
}

SplineCalculator::SplineCalculator(const std::string& fileName, const std::string& treeName, const std::vector<std::string>& friendFileNames, const std::vector<std::string>& friendTreeNames) 
    : df(createDataFrameWithFriends(fileName, treeName, friendFileNames, friendTreeNames))
{
    // df is now an RNode, initialized from an RDataFrame
    auto colNames = df.GetColumnNames();

    // 3. Loop over the list and print each name
    std::cout << "Column Names:" << std::endl;
    for (const auto& name : colNames) {
        std::cout << "- " << name << std::endl;
    }
}

ROOT::RDataFrame SplineCalculator::createDataFrameWithFriends(const std::string& fileName, const std::string& treeName, const std::vector<std::string>& friendFileNames, const std::vector<std::string>& friendTreeNames) {
    if (friendFileNames.size() != friendTreeNames.size()) {
        throw std::runtime_error("Number of friend files and friend trees must be the same.");
    }
    // Using a unique_ptr for the TFile to ensure it's closed, although TFile's destructor handles this.
    auto mainFile = new TFile(fileName.c_str(), "READ");
    if (!mainFile || mainFile->IsZombie()) {
        throw std::runtime_error("Could not open main file: " + fileName);
    }
    TTree *mainTree = nullptr;
    mainFile->GetObject(treeName.c_str(), mainTree);
    if (!mainTree) {
        throw std::runtime_error("Could not find tree: " + treeName + " in file: " + fileName);
    }
    for (size_t i = 0; i < friendFileNames.size(); ++i) {
        // The TTree::AddFriend documentation states that the friend TTree will be owned by the TChain/TTree.
        // The TFile for the friend can be closed after adding.
        mainTree->AddFriend(friendTreeNames[i].c_str(), friendFileNames[i].c_str());
    }

    return ROOT::RDataFrame(*mainTree);
}

ULong64_t SplineCalculator::getNEntries() {
    return df.Count().GetValue();
}

void SplineCalculator::setInterpolationType(InterpolationType type) {
    interpolationType = type;
}

void SplineCalculator::addBinningAxis(const std::string& variable, const std::vector<double>& edges) {
    binningAxes.push_back({variable, edges});
}

void SplineCalculator::addCategoricalBinningAxis(const std::string& variable, const std::vector<int>& categories) {
    if (categories.empty()) {
        throw std::runtime_error("Categorical axis must have at least one category.");
    }

    // Create a map from category value to a sequential bin index (0, 1, 2...)
    std::map<int, float> category_map;
    std::vector<double> edges;
    for (size_t i = 0; i < categories.size(); ++i) {
        category_map[categories[i]] = i;
        edges.push_back(static_cast<double>(i) - 0.5);
    }
    edges.push_back(static_cast<double>(categories.size()) - 0.5);

    //Need to replace the '.' in the variable name to avoid issues with RDataFrame
    std::string new_col_name = variable;
    std::replace(new_col_name.begin(), new_col_name.end(), '.', '_');
    new_col_name += "_categorical_index";
    categorical_maps[new_col_name] = category_map;

    // Define a new column that maps the original category value to the sequential index
    auto mapping_lambda = [map = categorical_maps[new_col_name]](int val) {
        return map.count(val) ? map.at(val) : -1.0f; // Return -1 for values not in map (will go to underflow bin)
    };
    df = df.Define(new_col_name, mapping_lambda, {variable});

    addBinningAxis(new_col_name, edges);
}

void SplineCalculator::setEventWeightColumn(const std::string& columnName) {
    eventWeightColumn = columnName;
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

   
    ROOT::RDF::THnDModel model("", "Template", naxes, binning_nbins, binning_edges);
    hist_template_ = model.GetHistogram();

    //Set manually the titles of the axes
    for(uint i = 0; i < binningAxes.size(); ++i){
        hist_template_->GetAxis(i)->SetTitle(binningAxes[i].variable.c_str());
    }

    // --- Pre-calculate strides for efficient global bin index calculation ---
    // The stride for an axis is the product of the number of bins of all preceding axes.
    std::vector<long long> strides(naxes);
    strides[0] = 1;
    for (uint i = 1; i < naxes; ++i) {
        strides[i] = hist_template_->GetAxis(i - 1)->GetNbins() + 2; // +2 for underflow and overflow
    }


     std::vector<ResultPtr_t> all_hists;

    // df is now the head of the graph, so we use it directly.
    auto df_with_bins = df.Define("global_bin_id", []() { return 0LL; }, {});

    if (eventWeightColumn.empty()) {
        eventWeightColumn = "eventWeight"; // Default name
        df_with_bins = df_with_bins.Define(eventWeightColumn, []() { return 1.0; });
    }

    // --- Iteratively build the global bin index ---
    for (uint i=0; i<naxes; ++i) {
        const auto& axis = binningAxes[i];
        // This lambda calculates the weighted contribution of a single axis to the global bin index.
        auto get_1d_bin_contribution_lambda = [hist_axis = hist_template_->GetAxis(i), stride = strides[i]](const float& value, const long long& prev_id) -> long long {
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
        ResultPtr_t hists = df_with_bins.Book<WeightVec_t, long long, double>(std::move(helper), {syst.vectorWeightBranch, "global_bin_id", eventWeightColumn});
        all_hists.push_back(hists);
    }
    
    // --- Trigger the event loop and build splines ---
    std::cout << "Starting event loop and spline construction..." << std::endl;
    all_splines_map.clear();

    for (size_t i = 0; i < systematics.size(); ++i) {
        const auto& syst = systematics[i];
        auto& syst_hists_ptr = all_hists[i];
        // The dereference here triggers the RDataFrame event loop for this systematic
        auto& syst_hists = *syst_hists_ptr;
        if (syst_hists.empty()) continue;

        // The histograms are 1D projections onto the global bin index.
        // All histograms for a given systematic have the same binning.
        long long n_global_bins = syst_hists[0]->GetNbins();

        //Printing the aggregated counts of each global_bin_id for debugging
        if (i == 0) { // Only print for the first systematic to avoid spam
            std::cout << "--- Debug: Aggregated counts for nominal case for systematic " << syst.name << " ---" << std::endl;
            auto nominal_hist = syst_hists[syst.nominalIndex];
            for (int bin = 0; bin <= n_global_bins; ++bin) {
                double content = nominal_hist->GetBinContent(bin);
                if (content > 0) {
                    std::cout << "  global_bin_id " << (bin) << ": " << content << " entries" << std::endl;
                }
            }
            std::cout << "-------------------------------------------------" << std::endl;
        }

        for (long long bin = 0; bin <= n_global_bins; ++bin) {
            std::vector<double> x_nodes, y_values;
            x_nodes.reserve(syst.paramNodes.size());
            y_values.reserve(syst.paramNodes.size());

            for (size_t j = 0; j < syst.paramNodes.size(); ++j) {
                const double nominal_value = syst_hists[syst.nominalIndex]->GetBinContent(bin);
                const double variation_value = syst_hists[j]->GetBinContent(bin);
                double ratio = 1.0; // Default to 1.0 (no change)
                if (nominal_value != 0) {
                    ratio = variation_value / nominal_value;
                }
                x_nodes.push_back(syst.paramNodes[j]);
                y_values.push_back(ratio);
            }

            // Check for "dummy" splines: not enough points or all y_values are 1.0
            if (x_nodes.size() < 2) { // Need at least 2 points for any meaningful interpolation
                continue;
            }
            bool all_ratios_one = true;
            for (double y : y_values) {
                if (y != 1.0) {
                    all_ratios_one = false;
                    break;
                }
            }
            if (all_ratios_one) {
                continue; // No systematic effect in this bin, don't store a spline
            }

            TString graph_name = TString::Format("%s_bin%lld", syst.name.c_str(), bin);

            // Create the appropriate graph or spline object
            if (interpolationType == InterpolationType::Spline) {
                auto spline = std::make_unique<TSpline3>(graph_name, x_nodes.data(), y_values.data(), x_nodes.size());
                all_splines_map[syst.name][bin - 1] = std::move(spline);
            } else { 
                // Linear interpolation is a TGraph
                auto graph = std::make_unique<TGraph>(x_nodes.size(), x_nodes.data(), y_values.data());
                graph->SetName(graph_name);
                graph->SetTitle(graph_name);
                all_splines_map[syst.name][bin - 1] = std::move(graph);
            }
        }
    }
    std::cout << "Spline construction complete." << std::endl;
}

void SplineCalculator::writeSplines(const std::string& outputFileName) const {
    TFile outFile(outputFileName.c_str(), "RECREATE");
    if (!outFile.IsOpen()) {
        throw std::runtime_error("Could not open output file: " + outputFileName);
    }

    long long total_possible_splines = 0;
    if (!binningAxes.empty()) {
        total_possible_splines = 1;
        for (const auto& axis : binningAxes) {
            total_possible_splines *= (axis.edges.size() - 1);
        }
    }
    total_possible_splines *= systematics.size();

    long long total_splines_created = 0;
    long long written_splines = 0;

    for (const auto& [syst_name, splines_map] : all_splines_map) {
        if (!outFile.Get(syst_name.c_str())) {
            outFile.mkdir(syst_name.c_str());
        }
        outFile.cd(syst_name.c_str());
        total_splines_created += splines_map.size();
        for (const auto& [bin_idx, spline_variant] : splines_map) {
            // Only write TSpline3 objects, as requested
            if (std::holds_alternative<std::unique_ptr<TSpline3>>(spline_variant)) {
                auto& spline_ptr = std::get<std::unique_ptr<TSpline3>>(spline_variant);
                if (spline_ptr) {
                    spline_ptr->Write(spline_ptr->GetName());
                    written_splines++;
                }
            }
        }
        outFile.cd("..");
    }

    //Also write the histogram template for binning info
    outFile.cd();
    if (hist_template_) {
        hist_template_->SetName("BinningTemplate");
        hist_template_->Write("BinningTemplate");
    }

    outFile.Close();

    double percentage = (total_possible_splines > 0) ? (100.0 * total_splines_created / total_possible_splines) : 0.0;
    std::cout << "--- Spline Writing Summary ---" << std::endl;
    std::cout << "Total possible splines (bins * systematics): " << total_possible_splines << std::endl;
    std::cout << "Total non-dummy splines created: " << total_splines_created << " (" << std::fixed << std::setprecision(2) << percentage << "%)" << std::endl;
    std::cout << "Total TSpline3 objects written to file: " << written_splines << std::endl;
    std::cout << "------------------------------" << std::endl;
}

const std::map<std::string, std::map<int, SplineVariant>>& SplineCalculator::getSplines() const {
    return all_splines_map;
}

std::unique_ptr<SplineContainer> SplineCalculator::getSplineContainer(const std::string& systName) const {
    if (all_splines_map.find(systName) == all_splines_map.end()) {
        throw std::runtime_error("Systematic '" + systName + "' not found in SplineCalculator.");
    }

    if (!hist_template_) {
      throw std::runtime_error("Histogram template is not available. Have you called run()?");
    }
    auto spline_container = std::make_unique<SplineContainer>(hist_template_.get());

    // Populate the container by cloning the splines from all_splines_map
    const auto& splines_map_for_syst = all_splines_map.at(systName);
    for (const auto& pair : splines_map_for_syst) {
        const auto& bin_idx = pair.first;
        const auto& spline_variant = pair.second;
        std::visit([&spline_container, &bin_idx](const auto& spline_ptr){
            if(spline_ptr) {
                spline_container->AddSpline(bin_idx, std::unique_ptr<TGraph>(static_cast<TGraph*>(spline_ptr->Clone())));
            }
        }, spline_variant);
    }

    return spline_container;
}