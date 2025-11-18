#ifndef SPLINECONTAINER_H
#define SPLINECONTAINER_H
#include "THn.h"
#include "TGraph.h"
#include "TSpline.h"
#include <vector>
#include <map>
#include <variant>
#include <memory>
#include <string>

// Forward declaration for the variant type
using SplineVariant = std::variant<std::unique_ptr<TGraph>, std::unique_ptr<TSpline3>>;

class SplineContainer {
public:
    SplineContainer() : fHistogram(nullptr) {} // For ROOT I/O
    SplineContainer(const THnT<double>* modelHist);
    virtual ~SplineContainer();

    // Method to add a spline for a given global bin index
    void AddSpline(int globalBin, SplineVariant&& spline);
    void AddSpline(int globalBin, std::unique_ptr<TGraph> spline);

    // The main user-facing methods
    const SplineVariant* GetSpline(const std::vector<double>& coords) const;
    const SplineVariant* GetSplineByBin(int globalBin) const;

    // Helper to calculate global bin index from coordinates
    int GetGlobalBin(const std::vector<double>& coords) const;

private:
    std::unique_ptr<THnT<double>> fHistogram; // Owns the binning definition
    std::map<int, SplineVariant> fSplineMap;
};

#endif // SPLINECONTAINER_H