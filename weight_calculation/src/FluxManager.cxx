#include "FluxManager.h"

FluxManager::FluxManager(const std::map<Flavour, std::string> &fnames)
: _fnames(fnames)
{
    Process();
}

FluxManager::~FluxManager()
{
}

void FluxManager::Process(){
    for(auto const &[fl, fname] : _fnames){
        _fluxes[fl] = Flux(fname);
    } 
}

double FluxManager::GetFlux(Flavour fl, double energy, double costh, double phi) const{
    if(_fluxes.count(fl) == 0){
        std::cout << "No flux available for flavour: " << flavours[fl] << std::endl;
        abort();
    }

    return _fluxes.at(fl).GetFlux(energy, costh, phi);
}