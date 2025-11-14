#include "SplineCalculator.h"
#include <iostream>

int main() {
    try {
        // Enable multithreading before creating any RDataFrame objects
        // SplineCalculator::EnableMultiThreading();
        ROOT::EnableImplicitMT();
        ROOT::RDataFrame df("SystWeights", "/Users/pgranger/test-ff/atmospherics-tools/spline_calculation/nusyst_new_sum.root");
        
        // 1. Create and book histograms for the "MaCCRES" systematic
        MultiHistoNDHelper<double, 2> ma_helper{"ma_hists", // Name
                   "A THn with 2 dimensions",        // Title
                    7,                  // Number of histograms
                   {5, 5},                     // NBins
                   {0, 0},            // Axes min values
                   {10., 10.}};               // Axes max values

        // By providing the column types as template arguments, we avoid the need for the JIT compiler
        // and gInterpreter. The compiler generates the action code at compile time.
        // NOTE: Replace these types with the actual types of your columns.
        using WeightVec_t = ROOT::RVec<double>;
        // The column order must match the Exec method: weights first, then binning variables.
        auto ma_hists = df.Book<WeightVec_t, float, float>(std::move(ma_helper), {"MaCCRES", "Ecalo", "Pmu_x"});

        // 2. Create and book histograms for the "MvCCRES" systematic
        MultiHistoNDHelper<double, 2> mv_helper{"mv_hists", // Name
                   "A THn with 2 dimensions for MvCCRES",        // Title
                    7,                  // Number of histograms
                   {5, 5},                     // NBins
                   {0, 0},            // Axes min values
                   {10., 10.}};               // Axes max values
        auto mv_hists = df.Book<WeightVec_t, float, float>(std::move(mv_helper), {"MvCCRES", "Ecalo", "Pmu_x"});

        // The event loop is triggered when you access the results.
        // Let's print the contents of the first histogram from each set.
        std::cout << "--- MaCCRES Histograms ---" << std::endl;
        for (const auto& hist : *ma_hists) {
            hist->Print("all");
        }
        
        std::cout << "\n--- MvCCRES Histograms ---" << std::endl;
        for (const auto& hist : *mv_hists) {
            hist->Print("all");
        }

        // 1. Create a calculator instance
        SplineCalculator calculator("/Users/pgranger/test-ff/atmospherics-tools/spline_calculation/nusyst_new_sum.root", "SystWeights");
        // // // 2. Configure binning (unchanged)
        calculator.addBinningAxis("Ecalo", {0.0, 1.0, 2.0, 5.0, 10.0});
        calculator.addBinningAxis("Pmu_x", {-1, -0.5, 0, 0.5, 1});
        // // calculator.addBinningAxis("reco_lepton_angle", {0.0, 0.5, 1.0, 2.5, 3.14});

        // // // 3. Add a systematic from a vector branch
        // // // The parameter nodes should correspond to the elements in the vector branch
        // calculator.addSystematic("MaCCRES",
        //                          {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
        //                          "MaCCRES",               // The branch with the vector of weights
        //                          2);                          // The index of the nominal weight (0.0 sigma)
        // calculator.addSystematic("ZExpA1CCQE",
        //                          {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
        //                          "ZExpA1CCQE",               // The branch with the vector of weights
        //                          2);                          // The index of the nominal weight (0.0 sigma)
        // calculator.addSystematic("ZExpA2CCQE",
        //                          {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
        //                          "ZExpA2CCQE",               // The branch with the vector of weights
        //                          2);                          // The index of the nominal weight (0.0 sigma)
        // calculator.addSystematic("ZExpA3CCQE",
        //                          {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
        //                          "ZExpA3CCQE",               // The branch with the vector of weights
        //                          2);                          // The index of the nominal weight (0.0 sigma)
        // calculator.addSystematic("ZExpA4CCQE",
        //                          {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
        //                          "ZExpA4CCQE",               // The branch with the vector of weights
        //                          2);                          // The index of the nominal weight (0.0 sigma)

        // // // 4. Run the calculation
        // calculator.run();

        // // // 5. Write the results to a file
        // // calculator.writeSplines("splines_vector_branch.root");

        // // std::cout << "Successfully created splines in splines_vector_branch.root" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
