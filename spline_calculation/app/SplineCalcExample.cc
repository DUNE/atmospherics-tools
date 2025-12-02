#include <iostream>
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include <ROOT/RLogger.hxx>
#include "yaml-cpp/yaml.h"
#include <string>
#include "SplineCalculator.h"


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

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <config.yaml>" << std::endl;
        return 1;
    }

    std::string config_path = argv[1];
    YAML::Node config = YAML::LoadFile(config_path);

    // --- Configure ROOT based on YAML ---
    if (config["multithreading"] && config["multithreading"].as<bool>()) {
        ROOT::EnableImplicitMT();
        std::cout << "Multithreading enabled." << std::endl;
    }

    auto logLevel = ROOT::ELogLevel::kInfo; // Default
    if (config["log_level"]) {
        std::string levelStr = config["log_level"].as<std::string>();
        if (levelStr == "Debug") logLevel = ROOT::ELogLevel::kDebug;
        else if (levelStr == "Info") logLevel = ROOT::ELogLevel::kInfo;
        else if (levelStr == "Warning") logLevel = ROOT::ELogLevel::kWarning;
        else if (levelStr == "Error") logLevel = ROOT::ELogLevel::kError;
    }
    auto verbosity = ROOT::RLogScopedVerbosity(ROOT::Detail::RDF::RDFLogChannel(), logLevel);

    // --- Load Data ---
    const auto& input_data = config["input_data"];
    SplineCalculator calculator(
        input_data["main_file"].as<std::string>(),
        input_data["main_tree"].as<std::string>(),
        input_data["friend_files"].as<std::vector<std::string>>(),
        input_data["friend_trees"].as<std::vector<std::string>>()
    );

    // --- Configure Binning ---
    std::cout << "Configuring binning axes..." << std::endl;
    for (const auto& axis_node : config["binning_axes"]) {
        std::string type = axis_node["type"].as<std::string>();
        std::string variable = axis_node["variable"].as<std::string>();
        if (type == "categorical") {
            auto categories = axis_node["categories"].as<std::vector<int>>();
            calculator.addCategoricalBinningAxis(variable, categories);
            std::cout << "  Added categorical axis for '" << variable << "'" << std::endl;
        } else if (type == "continuous") {
            auto edges = axis_node["edges"].as<std::vector<double>>();
            calculator.addBinningAxis(variable, edges);
            std::cout << "  Added continuous axis for '" << variable << "'" << std::endl;
        }
    }

    // --- Configure Interpolation Type ---
    if (config["interpolation_type"]) {
        std::string type = config["interpolation_type"].as<std::string>();
        if (type == "Linear") {
            calculator.setInterpolationType(InterpolationType::Linear);
            std::cout << "Set interpolation type to Linear." << std::endl;
        } else {
            calculator.setInterpolationType(InterpolationType::Spline);
            std::cout << "Set interpolation type to Spline." << std::endl;
        }
    }

    // --- Load and Add Systematics ---
    std::cout << "Loading systematics..." << std::endl;
    const auto& syst_config = config["systematics"];
    std::vector<Systematic> systematics = getSystematicsFromRootFile(
        syst_config["file"].as<std::string>(),
        syst_config["tree"].as<std::string>()
    );
    for (const auto& syst : systematics) {
        calculator.addSystematic(syst);
    }
    std::cout << "Added " << systematics.size() << " systematics." << std::endl;

    // --- Selections from Config ---
    if (config["selections"]) {
        std::cout << "Applying selections from config..." << std::endl;
        for (const auto& selection_node : config["selections"]) {
            auto selection_string = selection_node.as<std::string>();
            calculator.addSelection(selection_string);
            std::cout << "  Applied selection: \"" << selection_string << "\"" << std::endl;
        }
    }

    // --- Run and Write Output ---
    std::cout << "Running calculation..." << std::endl;
    calculator.run();

    const std::string output_file = config["output_file"].as<std::string>();
    std::cout << "Writing splines to " << output_file << "..." << std::endl;
    calculator.writeSplines(output_file);

    std::cout << "Successfully created splines in " << output_file << std::endl;

    return 0;
}
