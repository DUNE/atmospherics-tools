#include <TFile.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TObject.h>
#include <TH2D.h>
#include <TAxis.h>

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

struct MergeConfig {
    std::vector<std::vector<int>> rows;    // Y bins
    std::vector<std::vector<int>> columns; // X bins
};

struct AxisMapping {
    std::vector<double> edges;
    std::vector<int> oldToNew;
};

AxisMapping makeAxisMapping(
    const TAxis* axis,
    const std::vector<std::vector<int>>& groups)
{
    const int n = axis->GetNbins();

    std::map<int, int> groupEnd;
    std::set<int> occupied;

    for (const auto& group : groups) {
        if (group.size() < 2) {
            throw std::runtime_error(
                "Each merge group must contain at least two bins.");
        }

        for (std::size_t i = 0; i < group.size(); ++i) {
            const int bin = group[i];

            if (bin < 1 || bin > n) {
                throw std::runtime_error(
                    "Merge bin index is outside the valid range.");
            }

            if (i > 0 && bin != group[i - 1] + 1) {
                throw std::runtime_error(
                    "Merge groups must contain consecutive bins.");
            }

            if (!occupied.insert(bin).second) {
                throw std::runtime_error(
                    "Overlapping merge groups are not allowed.");
            }
        }

        groupEnd[group.front()] = group.back();
    }

    AxisMapping result;
    result.oldToNew.resize(n + 2, 0);

    result.edges.push_back(axis->GetBinLowEdge(1));

    int newBin = 0;
    int oldBin = 1;

    while (oldBin <= n) {
        int lastBin = oldBin;

        const auto it = groupEnd.find(oldBin);
        if (it != groupEnd.end()) {
            lastBin = it->second;
        }

        ++newBin;

        for (int b = oldBin; b <= lastBin; ++b) {
            result.oldToNew[b] = newBin;
        }

        result.edges.push_back(axis->GetBinUpEdge(lastBin));

        oldBin = lastBin + 1;
    }

    result.oldToNew[0] = 0;
    result.oldToNew[n + 1] = newBin + 1;

    return result;
}

std::unique_ptr<TH2D> mergeHistogram(
    const TH2D* source,
    const MergeConfig& config)
{
    const AxisMapping xMap =
        makeAxisMapping(source->GetXaxis(), config.columns);

    const AxisMapping yMap =
        makeAxisMapping(source->GetYaxis(), config.rows);

    const int nxNew = static_cast<int>(xMap.edges.size()) - 1;
    const int nyNew = static_cast<int>(yMap.edges.size()) - 1;

    auto result = std::make_unique<TH2D>(
        source->GetName(),
        source->GetTitle(),
        nxNew, xMap.edges.data(),
        nyNew, yMap.edges.data());

    result->SetDirectory(nullptr);

    result->GetXaxis()->SetTitle(source->GetXaxis()->GetTitle());
    result->GetYaxis()->SetTitle(source->GetYaxis()->GetTitle());
    result->GetZaxis()->SetTitle(source->GetZaxis()->GetTitle());
    result->SetLineColor(source->GetLineColor());
    result->SetLineStyle(source->GetLineStyle());
    result->SetLineWidth(source->GetLineWidth());
    result->SetFillColor(source->GetFillColor());
    result->SetFillStyle(source->GetFillStyle());
    result->SetMarkerColor(source->GetMarkerColor());
    result->SetMarkerStyle(source->GetMarkerStyle());
    result->SetMarkerSize(source->GetMarkerSize());
    //result->SetStats(source->GetStats());
    result->SetOption(source->GetOption());
    for (int x = 0; x <= source->GetNbinsX() + 1; ++x) {
        for (int y = 0; y <= source->GetNbinsY() + 1; ++y) {
            const int newX = xMap.oldToNew[x];
            const int newY = yMap.oldToNew[y];

            const double content = source->GetBinContent(x, y);
            const double error = source->GetBinError(x, y);

            const int targetBin =
                result->GetBin(newX, newY);

            result->SetBinContent(
                targetBin,
                result->GetBinContent(targetBin) + content);

            const double previousError =
                result->GetBinError(targetBin);

            result->SetBinError(
                targetBin,
                std::sqrt(previousError * previousError
                          + error * error));
        }
    }
    result->SetEntries(source->GetEntries());

    return result;
}

std::map<std::string, MergeConfig> loadConfig(
    const std::string& yamlFile,
    std::string& inputFile,
    std::string& outputFile)
{
    const YAML::Node root = YAML::LoadFile(yamlFile);

    inputFile = root["input_file"].as<std::string>();
    outputFile = root["output_file"].as<std::string>();

    std::map<std::string, MergeConfig> configs;

    for (const auto& hist : root["histograms"]) {
        const std::string name = hist["name"].as<std::string>();

        if (configs.count(name)) {
            throw std::runtime_error(
                "Duplicate histogram configuration: " + name);
        }

        MergeConfig config;

        auto readGroups = [](const YAML::Node& node,
                             std::vector<std::vector<int>>& groups) {
            if (!node) {
                return;
            }

            for (const auto& groupNode : node) {
                std::vector<int> group;

                for (const auto& bin : groupNode) {
                    group.push_back(bin.as<int>());
                }

                groups.push_back(std::move(group));
            }
        };

        readGroups(hist["merge_rows"], config.rows);
        readGroups(hist["merge_columns"], config.columns);

        configs.emplace(name, std::move(config));
    }

    return configs;
}

void copyDirectory(
    TDirectory* inputDir,
    TDirectory* outputDir,
    const std::string& prefix,
    const std::map<std::string, MergeConfig>& configs)
{
    TIter nextKey(inputDir->GetListOfKeys());
    TKey* key = nullptr;
    std::set<std::string> processed;

    while ((key = static_cast<TKey*>(nextKey()))) {
        const std::string name = key->GetName();

        if (!processed.insert(name).second) {
            continue;
        }

        const std::string path =
            prefix.empty() ? name : prefix + "/" + name;

        std::unique_ptr<TObject> object(key->ReadObj());

        if (!object) {
            throw std::runtime_error(
                "Could not read ROOT object: " + path);
        }

        if (auto* inputSubdir = dynamic_cast<TDirectory*>(object.get())) {
            TDirectory* outputSubdir =
                outputDir->GetDirectory(name.c_str());

            if (!outputSubdir) {
                outputSubdir =
                    outputDir->mkdir(name.c_str());

                if (!outputSubdir) {
                    throw std::runtime_error(
                        "Could not create output directory: " + path);
                }
            }

            copyDirectory(
                inputSubdir, outputSubdir, path, configs);

            continue;
        }

        const auto configIt = configs.find(path);

        if (configIt != configs.end()) {
            auto* hist = dynamic_cast<TH2D*>(object.get());

            if (!hist) {
                throw std::runtime_error(
                    "Configured object is not a TH2D: " + path);
            }

            auto merged = mergeHistogram(hist, configIt->second);

            outputDir->cd();
            if (merged->Write(name.c_str(), TObject::kOverwrite) <= 0) {
                throw std::runtime_error(
                    "Failed to write merged histogram: " + path);
            }

            std::cout << "Merged histogram: " << path << '\n';
        } else {
            outputDir->cd();

            if (object->Write(name.c_str(), TObject::kOverwrite) <= 0) {
                throw std::runtime_error(
                    "Failed to copy ROOT object: " + path);
            }
        }
    }
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr
            << "Usage: " << argv[0]
            << " config/bin_merging.yaml\n";
        return 1;
    }

    try {
        std::string inputFileName;
        std::string outputFileName;

        const auto configs = loadConfig(
            argv[1], inputFileName, outputFileName);

        if (inputFileName == outputFileName) {
            throw std::runtime_error(
                "Input and output files must be different.");
        }

        std::unique_ptr<TFile> inputFile(
            TFile::Open(inputFileName.c_str(), "READ"));

        if (!inputFile || inputFile->IsZombie()) {
            throw std::runtime_error(
                "Could not open input ROOT file: " + inputFileName);
        }

        std::unique_ptr<TFile> outputFile(
            TFile::Open(outputFileName.c_str(), "RECREATE"));

        if (!outputFile || outputFile->IsZombie()) {
            throw std::runtime_error(
                "Could not create output ROOT file: " + outputFileName);
        }

        copyDirectory(
            inputFile.get(), outputFile.get(), "", configs);

        outputFile->Write();
        outputFile->Close();
        inputFile->Close();

        std::cout << "Output written to: "
                  << outputFileName << '\n';

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
