#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH2D.h>
#include <TKey.h>
#include <TObject.h>
#include <TROOT.h>
#include <TStyle.h>

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

struct HistogramConfig {
    double errorThreshold;
};

using ConfigMap = std::map<std::string, HistogramConfig>;

ConfigMap loadConfig(
    const std::string& yamlFile,
    std::string& inputFile,
    std::string& outputFile,
    std::string& pdfFile)
{
    const YAML::Node root = YAML::LoadFile(yamlFile);

    inputFile = root["input_file"].as<std::string>();
    outputFile = root["output_file"].as<std::string>();

    pdfFile = root["pdf_file"]
        ? root["pdf_file"].as<std::string>()
        : "bin_error_threshold.pdf";

    if (inputFile == outputFile) {
        throw std::runtime_error(
            "Input and output ROOT files must be different.");
    }

    if (pdfFile.empty()) {
        throw std::runtime_error("PDF filename cannot be empty.");
    }

    ConfigMap configs;

    for (const auto& hist : root["histograms"]) {
        const std::string name = hist["name"].as<std::string>();
        const double threshold =
            hist["error_threshold"].as<double>();

        if (name.empty()) {
            throw std::runtime_error(
                "Histogram name cannot be empty.");
        }

        if (!std::isfinite(threshold) || threshold < 0.0) {
            throw std::runtime_error(
                "Error threshold must be finite and non-negative: "
                + name);
        }

        if (!configs.emplace(
                name, HistogramConfig{threshold}).second) {
            throw std::runtime_error(
                "Duplicate histogram configuration: " + name);
        }
    }

    if (configs.empty()) {
        throw std::runtime_error(
            "No histograms specified in the YAML configuration.");
    }

    return configs;
}

int processHistogram(TH2D* hist, double threshold)
{
    int modifiedBins = 0;

    for (int x = 1; x <= hist->GetNbinsX(); ++x) {
        for (int y = 1; y <= hist->GetNbinsY(); ++y) {
            const double content = hist->GetBinContent(x, y);
            const double error = hist->GetBinError(x, y);

            if (error > threshold){
                if (error > 0.5 || error > std::abs(content - 1.0)) {
                    hist->SetBinContent(x, y, 1.0);
                    ++modifiedBins;
                    if (error > 0.5) std::cout<< "WARNING LARGE ERROR FOR BIN x=" << x << ", y= " << y << std::endl;
                }
            }
        }
    }

    return modifiedBins;
}

TDirectory* makeOutputDirectory(
    TFile* output,
    const std::string& path)
{
    TDirectory* current = output;

    if (path.empty()) {
        return current;
    }

    std::size_t start = 0;

    while (start < path.size()) {
        const std::size_t slash = path.find('/', start);

        const std::string component = path.substr(
            start,
            slash == std::string::npos
                ? std::string::npos
                : slash - start);

        if (component.empty()) {
            throw std::runtime_error(
                "Invalid ROOT directory path: " + path);
        }

        TDirectory* next =
            current->GetDirectory(component.c_str());

        if (!next) {
            next = current->mkdir(component.c_str());
        }

        if (!next) {
            throw std::runtime_error(
                "Cannot create ROOT directory: " + path);
        }

        current = next;

        if (slash == std::string::npos) {
            break;
        }

        start = slash + 1;
    }

    return current;
}

std::pair<std::string, std::string> splitHistogramPath(
    const std::string& path)
{
    const std::size_t slash = path.find_last_of('/');

    if (slash == std::string::npos) {
        return {"", path};
    }

    if (slash == 0 || slash + 1 == path.size()) {
        throw std::runtime_error(
            "Invalid histogram path: " + path);
    }

    return {
        path.substr(0, slash),
        path.substr(slash + 1)
    };
}


void addPdfPage(
    TCanvas& canvas,
    TH2D* hist,
    const std::string& title,
    const std::string& pdfFile)
{
    const std::string originalTitle = hist->GetTitle();
    hist->SetTitle(title.c_str());

    canvas.cd();
    canvas.Clear();

    canvas.SetLogy();
    canvas.SetLogx();

    hist->Draw("COLZ TEXT");
    canvas.Modified();
    canvas.Update();
    canvas.Print(pdfFile.c_str());

    hist->SetTitle(originalTitle.c_str());
}

std::unique_ptr<TH2D> cloneHistogram(
    const TH2D* source,
    const std::string& cloneName)
{
    auto* cloned = dynamic_cast<TH2D*>(
        source->Clone(cloneName.c_str()));

    if (!cloned) {
        throw std::runtime_error(
            "Failed to clone histogram: "
            + std::string(source->GetName()));
    }

    cloned->SetDirectory(nullptr);
    return std::unique_ptr<TH2D>(cloned);
}

int main(int argc, char* argv[])
{
    if (argc != 2) {
        std::cerr
            << "Usage: " << argv[0]
            << " config/bin_error_threshold.yaml\n";
        return 1;
    }

    try {
        std::string inputFileName;
        std::string outputFileName;
        std::string pdfFileName;

        const ConfigMap configs = loadConfig(
            argv[1],
            inputFileName,
            outputFileName,
            pdfFileName);
       if (pdfFileName == inputFileName ||
            pdfFileName == outputFileName) {
            throw std::runtime_error(
                "The PDF filename must differ from both ROOT filenames.");
        }

        std::unique_ptr<TFile> input(
            TFile::Open(inputFileName.c_str(), "READ"));

        if (!input || input->IsZombie()) {
            throw std::runtime_error(
                "Cannot open input ROOT file: " + inputFileName);
        }

        std::unique_ptr<TFile> output(
            TFile::Open(outputFileName.c_str(), "RECREATE"));

        if (!output || output->IsZombie()) {
            throw std::runtime_error(
                "Cannot create output ROOT file: " + outputFileName);
        }
        TCanvas canvas("bin_error_canvas",
                       "Bin error threshold",
                       1200, 900);

        canvas.SetRightMargin(0.15);
        canvas.SetLeftMargin(0.12);
        canvas.SetBottomMargin(0.12);

        gStyle->SetOptStat(0);
        gStyle->SetPaintTextFormat("g");

        canvas.Print((pdfFileName + "[").c_str());

        for (const auto& entry : configs) {
            const std::string& path = entry.first;
            const double threshold = entry.second.errorThreshold;

            std::cout << "Processing " << path
                      << " (threshold = " << threshold << ")\n";
          TObject* object = input->Get(path.c_str());

            if (!object) {
                throw std::runtime_error(
                    "Configured histogram not found: " + path);
            }

            auto* source = dynamic_cast<TH2D*>(object);

            if (!source) {
                throw std::runtime_error(
                    "Configured object is not a TH2D: " + path);
            }
           auto before = cloneHistogram(
                source, std::string("before_") + source->GetName());

            auto after = cloneHistogram(
                source, source->GetName());

            const int modifiedBins =
                processHistogram(after.get(), threshold);

            std::cout << "  Modified bins: "
                      << modifiedBins << '\n';
          addPdfPage(
                canvas,
                before.get(),
                path + " (before modification)",
                pdfFileName);

            addPdfPage(
                canvas,
                after.get(),
                path + " (after modification)",
                pdfFileName);

           const auto pathParts = splitHistogramPath(path);

            TDirectory* outputDir = makeOutputDirectory(
                output.get(), pathParts.first);

            outputDir->cd();

            if (after->Write(
                    pathParts.second.c_str(),
                    TObject::kOverwrite) <= 0) {
                throw std::runtime_error(
                    "Failed to write histogram: " + path);
            }
        }
       canvas.Print((pdfFileName + "]").c_str());

        output->Write();
        output->Close();
        input->Close();

        std::cout << "\nOutput ROOT file: " << outputFileName
                  << "\nPDF report: " << pdfFileName << '\n';

    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
