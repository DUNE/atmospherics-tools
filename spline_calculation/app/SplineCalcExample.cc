#include "SplineCalculator.h"
#include <iostream>
#include "duneanaobj/StandardRecord/SRGlobal.h"

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
                     // The index of the nominal weight (0.0 sigma)


        // ROOT::RDataFrame df("SystWeights", "/Users/pgranger/test-ff/atmospherics-tools/spline_calculation/nusyst_new_sum.root");
        
        // // 1. Create and book histograms for the "MaCCRES" systematic
        // MultiHistoNDHelper<double, 2> ma_helper{"ma_hists", // Name
        //            "A THn with 2 dimensions",        // Title
        //             7,                  // Number of histograms
        //            {5, 5},                     // NBins
        //            {0, 0},            // Axes min values
        //            {10., 10.}};               // Axes max values

        // // By providing the column types as template arguments, we avoid the need for the JIT compiler
        // // and gInterpreter. The compiler generates the action code at compile time.
        // // NOTE: Replace these types with the actual types of your columns.
        // using WeightVec_t = ROOT::RVec<double>;
        // // The column order must match the Exec method: weights first, then binning variables.
        // auto ma_hists = df.Book<WeightVec_t, float, float>(std::move(ma_helper), {"MaCCRES", "Ecalo", "Pmu_x"});

        // // 2. Create and book histograms for the "MvCCRES" systematic
        // MultiHistoNDHelper<double, 2> mv_helper{"mv_hists", // Name
        //            "A THn with 2 dimensions for MvCCRES",        // Title
        //             7,                  // Number of histograms
        //            {5, 5},                     // NBins
        //            {0, 0},            // Axes min values
        //            {10., 10.}};               // Axes max values
        // auto mv_hists = df.Book<WeightVec_t, float, float>(std::move(mv_helper), {"MvCCRES", "Ecalo", "Pmu_x"});

        // // The event loop is triggered when you access the results.
        // // Let's print the contents of the first histogram from each set.
        // std::cout << "--- MaCCRES Histograms ---" << std::endl;
        // for (const auto& hist : *ma_hists) {
        //     hist->Print("all");
        // }
        
        // std::cout << "\n--- MvCCRES Histograms ---" << std::endl;
        // for (const auto& hist : *mv_hists) {
        //     hist->Print("all");
        // }

        
        // // calculator.addBinningAxis("reco_lepton_angle", {0.0, 0.5, 1.0, 2.5, 3.14});

        // // // 4. Run the calculation
        calculator.run();

        // // // 5. Write the results to a file
        // // calculator.writeSplines("splines_vector_branch.root");

        // // std::cout << "Successfully created splines in splines_vector_branch.root" << std::endl;

    return 0;
}
