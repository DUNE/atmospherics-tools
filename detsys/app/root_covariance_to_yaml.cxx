#include <TFile.h>
#include <TMatrixT.h>
#include <TVectorT.h>

#include <yaml-cpp/yaml.h>

#include <map>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <algorithm>
#include <fstream>

#include <TString.h>



struct KinematicCut {
    std::string Name;
    double LowerBound;
    double UpperBound;
};


const std::map<std::string, std::string> KinematicNameMap = {
    {"kNuERec",     "RecoNeutrinoEnergy"},
    {"kNuRecoCosZ", "RecoCosineZ"},
    {"kNuTrueE",    "TrueNeutrinoEnergy"},
    {"kNuTrueCosZ", "TrueCosineZ"}
};

bool TryParseDouble(const std::string& text, double& value)
{
    char* end = nullptr;
    value = std::strtod(text.c_str(), &end);

    return end != text.c_str()
        && *end == '\0'
        && std::isfinite(value);
}

bool TryParseFlavorCut(
    const std::string& parameterName,
    KinematicCut& cut)
{
    std::vector<std::string> tokens;
    std::stringstream stream(parameterName);
    std::string token;

    while (std::getline(stream, token, '_')) {
        tokens.push_back(token);
    }

    for (const auto& part : tokens) {
        if (part == "NuE") {
            cut = {"SelFlavour", -0.5, 0.5};
            return true;
        }

        if (part == "NuMu") {
            cut = {"SelFlavour", 0.5, 1.5};
            return true;
        }
        if (part == "NC") {
            cut = {"SelFlavour", 1.5, 2.5};
            return true;
        }
    }

    return false;
}

std::vector<KinematicCut> ParseKinematicCuts(
    const std::string& parameterName)
{
    std::vector<KinematicCut> cuts;
    KinematicCut flavorCut;

    if (TryParseFlavorCut(parameterName, flavorCut)) {
        cuts.push_back(flavorCut);
    }
    std::vector<std::string> tokens;
    std::stringstream stream(parameterName);
    std::string token;

    while (std::getline(stream, token, '_')) {
        tokens.push_back(token);
    }
    for (std::size_t i = 0; i < tokens.size(); ++i) {
        const auto mapping = KinematicNameMap.find(tokens[i]);
        if (mapping == KinematicNameMap.end()) {
            continue;
        }
        if (i + 2 >= tokens.size()) {
            throw std::runtime_error(
                "Missing kinematic bounds after '" + tokens[i] +
                "' in parameter: " + parameterName
            );
        }

        double lower = 0.0;
        double upper = 0.0;

        if (!TryParseDouble(tokens[i + 1], lower) ||
            !TryParseDouble(tokens[i + 2], upper)) {
            throw std::runtime_error(
                "Invalid kinematic bounds after '" + tokens[i] +
                "' in parameter: " + parameterName
            );
        }

        if (lower > upper) {
            throw std::runtime_error(
                "Lower kinematic bound exceeds upper bound in: " +
                parameterName
            );
        }
        cuts.push_back({
            mapping->second,
            lower,
            upper
        });
        i += 2;
    }

    return cuts;
}

struct Config {
    std::vector<KinematicCut> KinematicCuts;
    std::vector<std::string> SampleNames;
    std::vector<std::string> DetID;

    std::string ParameterGroup = "Osc";
    std::string Type = "Oscillation";

    bool FlatPrior = false;

    double Generated = 1.0;
    double PreFitValue = 1.0;
    double StepScaleMCMC = 1.0;
};

template <typename T>
T* GetRequiredObject(TFile& file, const char* key) {
    T* obj = nullptr;
    file.GetObject(key, obj);

    if (!obj) {
        throw std::runtime_error(
            std::string("Missing or incompatible ROOT object: ") + key
        );
    }

    return obj;
}

void ConvertRootToYaml(
    const std::string& inputFileName,
    const std::string& outputFileName,
    const Config& config
) {
    TFile inputFile(inputFileName.c_str(), "READ");

    if (inputFile.IsZombie()) {
        throw std::runtime_error(
            "Could not open ROOT file: " + inputFileName
        );
    }

    auto* physicsParameter = GetRequiredObject<TString>(
        inputFile, "detector_paremeter_name"
    );

    const std::string physicsName = physicsParameter->Data();

    auto* names = GetRequiredObject<std::vector<std::string>>(
        inputFile, "param_names"
    );

    auto* covariance = GetRequiredObject<TMatrixT<double>>(
        inputFile, "covariance"
    );

   auto* lowerBounds = GetRequiredObject<TVectorT<double>>(
        inputFile, "param_lb"
    );

    auto* upperBounds = GetRequiredObject<TVectorT<double>>(
        inputFile, "param_ub"
    );

    const int n = static_cast<int>(names->size());

    if (covariance->GetNrows() != n ||
        covariance->GetNcols() != n) {
        throw std::runtime_error(
            "Covariance matrix dimensions do not match param_names."
        );
    }

    if (lowerBounds->GetNrows() != n ||
        upperBounds->GetNrows() != n) {
        throw std::runtime_error(
            "Parameter bounds dimensions do not match param_names."
        );
    }

   std::vector<double> errors(n);

    for (int i = 0; i < n; ++i) {
        const double variance = (*covariance)(i, i);

        if (!std::isfinite(variance) || variance < -1e-12) {
            throw std::runtime_error(
                "Invalid covariance diagonal at index " +
                std::to_string(i)
            );
        }

        errors[i] = std::sqrt(std::max(0.0, variance));
    }

    YAML::Node document;
    document["Systematics"] = YAML::Node(YAML::NodeType::Sequence);

    std::vector<std::string> fullNames;
    fullNames.reserve(names->size());

    for (const auto& name : *names) {
        fullNames.push_back(physicsName + "_" + name);
    }

    for (int i = 0; i < n; ++i) {
        //const std::string& name = names->at(i);
        const std::string name = fullNames[i];
        YAML::Node systematic;

        systematic["Names"]["FancyName"] = name;
        systematic["Names"]["ParameterName"] = name;

        if (!config.DetID.empty()) {
            for (const auto& id : config.DetID) {
                systematic["DetID"].push_back(id);
            }
        } else {
            for (const auto& sample : config.SampleNames) {
                systematic["SampleNames"].push_back(sample);
            }
        }

        systematic["Error"] = errors[i];
        systematic["FlatPrior"] = config.FlatPrior;


        const auto cuts = ParseKinematicCuts(name);

        for (const auto& cut : cuts) {
            YAML::Node bounds(YAML::NodeType::Sequence);
            bounds.push_back(cut.LowerBound);
    	    bounds.push_back(cut.UpperBound);

    	    YAML::Node cutEntry;
    	    cutEntry[cut.Name] = bounds;

    	    systematic["KinematicCuts"].push_back(cutEntry);
	}	


        systematic["ParameterBounds"].push_back((*lowerBounds)[i]);
        systematic["ParameterBounds"].push_back((*upperBounds)[i]);

        systematic["ParameterGroup"] = config.ParameterGroup;

        systematic["ParameterValues"]["Generated"] =
            config.Generated;
        systematic["ParameterValues"]["PreFitValue"] =
            config.PreFitValue;

        systematic["StepScale"]["MCMC"] = config.StepScaleMCMC;

        systematic["Type"] = config.Type;

        YAML::Node correlations(YAML::NodeType::Sequence);

        for (int j = 0; j < n; ++j) {
            double rho = 0.0;
            const double denominator = errors[i] * errors[j];

            if (denominator > 0.0) {
                rho = (*covariance)(i, j) / denominator;
            }

            if (!std::isfinite(rho)) {
                throw std::runtime_error(
                    "Invalid correlation for parameters " +
                    name + " and " + fullNames[j]
                );
            }

            YAML::Node entry;
            entry[fullNames[j]] = rho;
            correlations.push_back(entry);
        }

        systematic["Correlations"] = correlations;

        YAML::Node entry;
        entry["Systematic"] = systematic;

        document["Systematics"].push_back(entry);
    }

    YAML::Emitter emitter;
    emitter << document;

    if (!emitter.good()) {
        throw std::runtime_error(
            std::string("YAML serialization failed: ") +
            emitter.GetLastError()
        );
    }

    const std::string output = emitter.c_str();


   std::ofstream outputFile(outputFileName);

    if (!outputFile) {
        throw std::runtime_error(
            "Could not create output file: " + outputFileName
        );
    }

    outputFile << "---\n" << output << '\n';
    outputFile.close();

    inputFile.Close();

    std::cout << "Converted " << n << " parameters.\n"
              << "Output written to: " << outputFileName << '\n';
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr
            << "Usage: " << argv[0]
            << " input.root output.yaml\n";
        return 1;
    }

    try {
        Config config;
        config.SampleNames = {"ATM"};
        config.ParameterGroup = "DetNorm";
        config.Type = "Norm"; //To change, do not know to what
        config.FlatPrior = false;
        config.Generated = 1.0;
        config.PreFitValue = 1.0;
        config.StepScaleMCMC = 1.0;

       ConvertRootToYaml(argv[1], argv[2], config);
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
