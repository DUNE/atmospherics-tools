#include "SplineContainer.h"
#include <stdexcept>

SplineContainer::SplineContainer(const THnT<double>* modelHist) {
    if (modelHist) {
        fHistogram.reset(static_cast<THnT<double>*>(modelHist->Clone()));
        fHistogram->Reset(); // We only need the structure, not the content
    }
}

SplineContainer::~SplineContainer() = default;

void SplineContainer::AddSpline(int globalBin, SplineVariant&& spline) {
    fSplineMap[globalBin] = std::move(spline); // This should work if spline is SplineVariant
}

void SplineContainer::AddSpline(int globalBin, std::unique_ptr<TGraph> spline) {
    fSplineMap[globalBin] = std::move(spline);
}

const SplineVariant* SplineContainer::GetSpline(const std::vector<double>& coords) const {
    int globalBin = GetGlobalBin(coords);
    return GetSplineByBin(globalBin);
}

const SplineVariant* SplineContainer::GetSplineByBin(int globalBin) const {
    auto it = fSplineMap.find(globalBin);
    if (it != fSplineMap.end()) {
        return &it->second;
    }
    return nullptr; // Not found
}

int SplineContainer::GetGlobalBin(const std::vector<double>& coords) const {
    if (!fHistogram) {
        throw std::runtime_error("SplineContainer is not properly initialized with a histogram model.");
    }
    if (coords.size() != static_cast<size_t>(fHistogram->GetNdimensions())) {
        throw std::runtime_error("Coordinate dimensions do not match number of axes.");
    }

    // THn::GetBin returns a 1-based global bin index. We use 0-based indices internally.
    return fHistogram->GetBin(coords.data()) - 1;
}