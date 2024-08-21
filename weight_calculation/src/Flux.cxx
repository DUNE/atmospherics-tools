#include "Flux.h"

Flux::Flux(std::string filename)
{
    this->Open(filename);
}

Flux::~Flux()
{
}

void Flux::Open(std::string filename){
    std::cout << "Going to open flux file: " << filename << std::endl;

    TFile *f = new TFile(filename.c_str(), "READ");
    if(f == nullptr){
        std::cout << "Could not open file " << filename << std::endl;
        abort();
    }

    TH3D *flux;
    f->GetObject("flux", flux);

    if(flux == nullptr){
        std::cout << "Could not get the flux histogram" << std::endl;
        abort();
    }

    _flux = static_cast<TH3D*>(flux->Clone());
    _flux->SetDirectory(0);
    
    if(_flux->GetZaxis()->GetNbins() == 1){
        _is2D = true;
        _flux2D = static_cast<TH2D*>(flux->Project3D("yx"));
        _flux2D->SetDirectory(0);
    }
    
    f->Close();
}

double Flux::GetFlux(double energy, double costh, double phi) const{
    if(_is2D){ //no binning in phi, bilinear interpolation only so using the 2D hist
        return _flux2D->Interpolate(energy, costh);
    }
    else {
        return _flux->Interpolate(energy, costh, phi);
    }
}
