#include "SplineCalculator.h"
#include <iostream>

int main() {
    try {
        // Enable multithreading before creating any RDataFrame objects
        SplineCalculator::EnableMultiThreading();
        
        // 1. Create a calculator instance
        SplineCalculator calculator("/Users/pgranger/test-ff/atmospherics-tools/spline_calculation/nusyst_new_sum.root", "SystWeights");
        // // 2. Configure binning (unchanged)
        calculator.addBinningAxis("Ecalo", {0.0, 1.0, 2.0, 5.0, 10.0});
        calculator.addBinningAxis("Pmu_x", {-1, -0.5, 0, 0.5, 1});
        // calculator.addBinningAxis("reco_lepton_angle", {0.0, 0.5, 1.0, 2.5, 3.14});

        // // 3. Add a systematic from a vector branch
        // // The parameter nodes should correspond to the elements in the vector branch
        calculator.addSystematic("MaCCRES",
                                 {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
                                 "MaCCRES",               // The branch with the vector of weights
                                 2);                          // The index of the nominal weight (0.0 sigma)
        calculator.addSystematic("ZExpA1CCQE",
                                 {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
                                 "ZExpA1CCQE",               // The branch with the vector of weights
                                 2);                          // The index of the nominal weight (0.0 sigma)
        calculator.addSystematic("ZExpA2CCQE",
                                 {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
                                 "ZExpA2CCQE",               // The branch with the vector of weights
                                 2);                          // The index of the nominal weight (0.0 sigma)
        calculator.addSystematic("ZExpA3CCQE",
                                 {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
                                 "ZExpA3CCQE",               // The branch with the vector of weights
                                 2);                          // The index of the nominal weight (0.0 sigma)
        calculator.addSystematic("ZExpA4CCQE",
                                 {-2.0, -1.0, 0.0, 1.0, 2.0}, // 5 parameter nodes
                                 "ZExpA4CCQE",               // The branch with the vector of weights
                                 2);                          // The index of the nominal weight (0.0 sigma)

        // // 4. Run the calculation
        calculator.run();

        // // 5. Write the results to a file
        // calculator.writeSplines("splines_vector_branch.root");

        // std::cout << "Successfully created splines in splines_vector_branch.root" << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
