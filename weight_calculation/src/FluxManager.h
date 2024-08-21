#pragma once
#include "Flux.h"
#include <string>
#include <map>

enum Flavour{
    NuE = 12,
    NuMu = 14,
    NuTau = 16,
    NuEBar = -12,
    NuMuBar = -14,
    NuTauBar = -16
};

static std::map<Flavour, std::string> flavours = {
    {Flavour::NuE, "nue"},
    {Flavour::NuMu, "numu"},
    {Flavour::NuTau, "nutau"},
    {Flavour::NuEBar, "nuebar"},
    {Flavour::NuMuBar, "numubar"},
    {Flavour::NuTauBar, "nutaubar"},
};

class FluxManager
{
private:
    std::map<Flavour, Flux> _fluxes;
    const std::map<Flavour, std::string> _fnames;

    void Process();
public:
    FluxManager(const std::map<Flavour, std::string> &fnames);
    ~FluxManager();
    double GetFlux(Flavour fl, double energy, double costh, double phi) const;
};


