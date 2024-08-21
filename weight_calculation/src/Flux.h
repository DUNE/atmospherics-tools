#pragma once
#include <iostream>
#include "TFile.h"
#include "TH3D.h"
#include "TH2D.h"

class Flux
{
private:
    TH3D *_flux = nullptr;
    TH2D *_flux2D = nullptr;
    bool _is2D = false;
    void Open(std::string filename);

public:
    Flux(std::string filename);
    Flux() = default;
    ~Flux();
    double GetFlux(double energy, double costh, double phi) const;
};
