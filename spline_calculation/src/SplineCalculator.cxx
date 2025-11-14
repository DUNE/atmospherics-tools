#include "SplineCalculator.h"
#include <algorithm>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "TFile.h"


std::string SplineCalculator::getColumnType(const std::string& columnName) {
    // Use RDataFrame's GetColumnType method to detect actual type
    std::string rawType = df.GetColumnType(columnName);
    
    // Print debug information
    std::cout << "Column '" << columnName << "' raw type: '" << rawType << "'" << std::endl;
    
    return rawType;
}

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
    systematics.push_back({systName, systParamNodes, vectorWeightBranch, nominalIndex});
}

void SplineCalculator::run() {
    if (binningAxes.empty()) {
        throw std::runtime_error("Binning must be set before running.");
    }
    if (systematics.empty()) {
        throw std::runtime_error("At least one systematic must be added before running.");
    }

    // Pre-calculate total number of bins to avoid dynamic histogram range calculation
    long long total_bins = 1;
    std::vector<size_t> axis_nbins;
    axis_nbins.reserve(binningAxes.size());
    
    for (const auto& axis : binningAxes) {
        size_t nbins = axis.edges.size() - 1;
        axis_nbins.push_back(nbins);
        total_bins *= nbins;
    }

    // Pre-calculate strides for faster binning
    std::vector<long long> strides(binningAxes.size());
    strides[0] = 1;
    for (size_t i = 1; i < binningAxes.size(); ++i) {
        strides[i] = strides[i-1] * axis_nbins[i-1];
    }

    // Collect all binning variable names
    std::vector<std::string> binning_vars;
    binning_vars.reserve(binningAxes.size());
    for (const auto& axis : binningAxes) {
        binning_vars.push_back(axis.variable);
    }

    // OPTIMIZATION: Cache type detection to avoid redundant calls
    std::unordered_map<std::string, std::string> type_cache;
    auto get_cached_type = [&](const std::string& col) {
        auto it = type_cache.find(col);
        if (it == type_cache.end()) {
            it = type_cache.emplace(col, getColumnType(col)).first;
        }
        return it->second;
    };
    
    // Create a lambda that works with ROOT's automatic type system
    auto robust_lambda = [axes = this->binningAxes, strides](ROOT::RVecF vals_in) -> long long {
        // ROOT will automatically convert input types to float internally
        long long global_idx = 0;
        
        for (size_t i = 0; i < axes.size() && i < vals_in.size(); ++i) {
            const auto& axis = axes[i];
            const double val = static_cast<double>(vals_in[i]); // Only convert for comparison

            if (val < axis.edges.front() || val >= axis.edges.back()) {
                return -1LL;
            }

            auto it = std::lower_bound(axis.edges.begin(), axis.edges.end(), val);
            int local_idx = std::distance(axis.edges.begin(), it) - 1;
            global_idx += local_idx * strides[i];
        }
        return global_idx;
    };
    
    ROOT::RDF::RNode df_with_bin_idx = df.Define("bin_index", robust_lambda, binning_vars);

    // Filter for valid bins
    auto filtered_df = df_with_bin_idx.Filter("bin_index >= 0");

    // Pre-allocate histogram specifications to avoid recalculation
    const int n_bins = static_cast<int>(total_bins);
    const double bin_low = -0.5;
    const double bin_high = static_cast<double>(total_bins) - 0.5;

    // MAJOR OPTIMIZATION: Group systematics by weight branch to avoid redundant column definitions
    std::unordered_map<std::string, std::vector<size_t>> branch_to_systematics;
    for (size_t i = 0; i < systematics.size(); ++i) {
        if (!systematics[i].paramNodes.empty()) {
            branch_to_systematics[systematics[i].vectorWeightBranch].push_back(i);
        }
    }
    
    auto final_df = filtered_df;
    std::vector<std::string> all_nominal_cols, all_syst_cols;
    std::vector<size_t> syst_start_indices;
    
    // Define columns only once per unique weight branch
    for (const auto& [weight_branch, syst_indices] : branch_to_systematics) {
        get_cached_type(weight_branch); // Cache the type detection
        
        // Find all unique indices we need for this branch
        std::set<size_t> required_indices;
        for (size_t syst_idx : syst_indices) {
            const auto& syst = systematics[syst_idx];
            required_indices.insert(syst.nominalIndex);
            for (size_t i = 0; i < syst.paramNodes.size(); ++i) {
                required_indices.insert(i);
            }
        }
        
        // Define columns for all required indices at once
        std::unordered_map<size_t, std::string> index_to_column;
        for (size_t idx : required_indices) {
            std::string col_name = weight_branch + "_idx_" + std::to_string(idx);
            std::string expr = "Take(" + weight_branch + ", " + std::to_string(idx) + ")";
            final_df = final_df.Define(col_name, expr);
            index_to_column[idx] = std::move(col_name);
        }
        
        // Map systematic columns to the shared weight columns
        for (size_t syst_idx : syst_indices) {
            const auto& syst = systematics[syst_idx];
            syst_start_indices.push_back(all_syst_cols.size());
            
            // Nominal column
            all_nominal_cols.push_back(index_to_column[syst.nominalIndex]);
            
            // Systematic columns
            for (size_t i = 0; i < syst.paramNodes.size(); ++i) {
                all_syst_cols.push_back(index_to_column[i]);
            }
        }
    }
    
    // OPTIMIZATION: Book all histograms without triggering computation
    std::vector<ROOT::RDF::RResultPtr<TH1D>> all_nominal_histos, all_syst_histos;
    std::vector<size_t> syst_param_counts; // Track parameter counts for each systematic
    
    size_t syst_idx = 0;
    size_t col_offset = 0;
    
    all_nominal_histos.reserve(systematics.size());
    
    for (const auto& [weight_branch, syst_indices] : branch_to_systematics) {
        for (size_t orig_syst_idx : syst_indices) {
            const auto& syst = systematics[orig_syst_idx];
            if (syst.paramNodes.empty()) continue;

            syst_param_counts.push_back(syst.paramNodes.size());

            // Book nominal histogram
            auto nominal_histo = final_df.Histo1D(
                {("nominal_" + std::to_string(syst_idx)).c_str(), "nominal", n_bins, bin_low, bin_high}, 
                "bin_index", all_nominal_cols[syst_idx]);
            all_nominal_histos.push_back(std::move(nominal_histo));

            // Book all systematic histograms for this systematic
            for (size_t i = 0; i < syst.paramNodes.size(); ++i) {
                std::string histo_name = "syst_" + std::to_string(syst_idx) + "_" + std::to_string(i);
                auto histo = final_df.Histo1D(
                    {histo_name.c_str(), histo_name.c_str(), n_bins, bin_low, bin_high}, 
                    "bin_index", all_syst_cols[col_offset + i]);
                all_syst_histos.push_back(std::move(histo));
            }
            
            col_offset += syst.paramNodes.size();
            syst_idx++;
        }
    }

    // OPTIMIZATION: Trigger computation once for all histograms using batch processing
    std::cout << "Processing " << all_nominal_histos.size() << " systematics..." << std::endl;
    
    // Access all histograms to trigger the single computation pass
    std::vector<TH1D*> nominal_hist_ptrs, syst_hist_ptrs;
    nominal_hist_ptrs.reserve(all_nominal_histos.size());
    syst_hist_ptrs.reserve(all_syst_histos.size());
    
    for (auto& histo : all_nominal_histos) {
        nominal_hist_ptrs.push_back(histo.GetPtr());
    }
    for (auto& histo : all_syst_histos) {
        syst_hist_ptrs.push_back(histo.GetPtr());
    }

    // OPTIMIZATION: Process splines with batch operations and pre-allocated containers
    size_t processed_systematics = 0;
    col_offset = 0;
    
    for (const auto& [weight_branch, syst_indices] : branch_to_systematics) {
        for (size_t orig_syst_idx : syst_indices) {
            const auto& syst = systematics[orig_syst_idx];
            if (syst.paramNodes.empty()) continue;

            TH1D* nominal_hist = nominal_hist_ptrs[processed_systematics];
            const size_t param_count = syst_param_counts[processed_systematics];
            
            // Pre-allocate map and vectors for this systematic
            std::map<int, std::unique_ptr<TSpline3>> splines_for_this_syst;
            std::vector<double> avg_responses;
            avg_responses.reserve(param_count);
            
            // Pre-allocate title string buffer
            std::string spline_title;
            spline_title.reserve(64);

            // Get systematic histogram pointers for this systematic
            std::vector<TH1D*> current_syst_hist_ptrs;
            current_syst_hist_ptrs.reserve(param_count);
            for (size_t i = 0; i < param_count; ++i) {
                current_syst_hist_ptrs.push_back(syst_hist_ptrs[col_offset + i]);
            }

            // OPTIMIZATION: Vectorized bin processing with minimal allocations
            const int max_bin = nominal_hist->GetNbinsX();
            for (int bin = 1; bin <= max_bin; ++bin) {
                const double nominal_sum = nominal_hist->GetBinContent(bin);
                if (nominal_sum == 0.0) continue;

                // Fast vector reuse instead of reallocation
                avg_responses.clear();
                const double inv_nominal = 1.0 / nominal_sum; // Pre-compute division
                
                // Vectorized response calculation
                for (TH1D* syst_hist : current_syst_hist_ptrs) {
                    avg_responses.push_back(syst_hist->GetBinContent(bin) * inv_nominal);
                }

                // Fast bin index calculation
                const int actual_bin_idx = static_cast<int>(nominal_hist->GetBinCenter(bin));

                // Optimized spline title generation
                spline_title = "spline_";
                spline_title += syst.name;
                spline_title += "_bin_";
                spline_title += std::to_string(actual_bin_idx);
                
                auto spline = std::make_unique<TSpline3>(
                    spline_title.c_str(),
                    const_cast<double*>(syst.paramNodes.data()),
                    avg_responses.data(),
                    static_cast<int>(param_count)
                );
                splines_for_this_syst[actual_bin_idx] = std::move(spline);
            }

            all_splines_map[syst.name] = std::move(splines_for_this_syst);
            col_offset += param_count;
            processed_systematics++;
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