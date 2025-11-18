#include "SplineCalculator.h"
#include <iostream>
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include <ROOT/RLogger.hxx>
 


std::vector<Systematic> getSystematicsFromRootFile(const std::string& fileName, const std::string& treeName) {
    std::vector<Systematic> systematics;
    TFile file(fileName.c_str(), "READ");
    caf::SRGlobal* srGlobal = nullptr;
    TTree* tree = nullptr;
    file.GetObject(treeName.c_str(), tree);
    if (!tree) {
        throw std::runtime_error("Could not find tree: " + treeName);
    }

    tree->SetBranchAddress("SRGlobal", &srGlobal);
    uint nread = tree->GetEntry(0); // Read the first entry to get SRGlobal

    if (!srGlobal) {
        throw std::runtime_error("SRGlobal branch is null.");
    }

    std::cout << "Found " << srGlobal->wgts.params.size() << " systematic parameters." << std::endl;
    
    for (const auto& param : srGlobal->wgts.params) {
        Systematic syst;
        syst.name = param.name;
        syst.nominalIndex = (param.nshifts - 1) / 2; // Assuming nominal is in the middle
        // For simplicity, we create dummy parameter nodes here
        for (int i = 0; i < param.nshifts; ++i) {
            syst.paramNodes.push_back(static_cast<double>(i - syst.nominalIndex));
        }
        syst.vectorWeightBranch = param.name; // Assuming branch name matches systematic name
        systematics.emplace_back(std::move(syst));
    }

    file.Close();
    return systematics;
}

int main() {
        // Enable multithreading before creating any RDataFrame objects
        // SplineCalculator::EnableMultiThreading();
        ROOT::EnableImplicitMT();

        // this increases RDF's verbosity level as long as the `verbosity` variable is in scope
        auto verbosity = ROOT::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), ROOT::ELogLevel::kInfo);

        std::vector<Systematic> systematics = getSystematicsFromRootFile("/Users/pgranger/test-ff/atmospherics-tools/spline_calculation/nusyst_new_sum.root", "SRGlobal");
        for (const auto& syst : systematics) {
            std::cout << syst << std::endl;
        }

        // 1. Create a calculator instance
        SplineCalculator calculator("/Users/pgranger/test-ff/atmospherics-tools/spline_calculation/nusyst_new_sum.root", "SystWeights");
        // // // 2. Configure binning (unchanged)
        calculator.addBinningAxis("Ecalo", {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 5.0, 10.0});
        calculator.addBinningAxis("Pmu_x", {-1, -0.5, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1});
        calculator.addBinningAxis("Pmu_x", {-1, -0.5, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1});

        // // 3. Add a systematic from a vector branch
        // // The parameter nodes should correspond to the elements in the vector branch

        for (const auto& syst : systematics) {
            calculator.addSystematic(syst);
        }

        uint nentries = calculator.getNEntries();
        std::cout << "Number of entries: " << nentries << std::endl;

        std::vector<bool> dummy_column(nentries, true);
        for(uint i=0; i<nentries; ++i) {
            dummy_column[i] = (i % 2 == 0); // Example condition: even indices
        }
        calculator.addColumnFromVector("dummy_selection", dummy_column);

        // calculator.addSelection([](float Ecalo, float Pmu_x) {
        //     return (Ecalo > 0.2) && (Pmu_x > 0.0);
        // }, {"Ecalo", "Pmu_x"});
        calculator.addSelection([](bool sel) {
            return sel;
        }, {"dummy_selection"});

        // // // 4. Run the calculation
        calculator.run();

        // // // 5. Write the results to a file
        calculator.writeSplines("splines_vector_branch.root");

        std::cout << "Successfully created splines in splines_vector_branch.root" << std::endl;

    return 0;
}
