#include "Calculator.h"

template<typename T>
Calculator<T>::Calculator(const FluxManager& mgr, Reader<T>& rdr, Writer& wrt, float exposure_scaling)
: _mgr(mgr), _rdr(rdr), _wrt(wrt), _exposure_scaling(exposure_scaling)
{
}

template<typename T>
Calculator<T>::~Calculator()
{
}

template<>
Data<double> Calculator<float>::FetchData(){
    return float2double(_rdr.GetData());
}

template<>
Data<double> Calculator<double>::FetchData(){
    return _rdr.GetData();
}

template<typename T>
void Calculator<T>::Process(){
    double POT = _rdr.POT();
    std::cout << "Got POT of: " << POT << std::endl;
    while(_rdr.GetEntry()){
        Data<double> data = FetchData();
        double Enu = data.ev;
        double costh = data.NuMomY/sqrt(pow(data.NuMomX, 2) + pow(data.NuMomY, 2) + pow(data.NuMomZ, 2));
        double phi = 60; //To be modified in the case of 3D flux.

        Flavour nuE = (data.nuPDG > 0) ? Flavour::NuE : Flavour::NuEBar;
        Flavour nuMu = (data.nuPDG > 0) ? Flavour::NuMu : Flavour::NuMuBar;

        double nuE_flux = _mgr.GetFlux(nuE, Enu, costh, phi);
        double nuMu_flux = _mgr.GetFlux(nuMu, Enu, costh, phi);

        //We assume here that the numu flux was used for the generation
        
        double xsec_w = data.weight/nuMu_flux/POT*_exposure_scaling;
        double nuE_w = nuE_flux;
        double nuMu_w = nuMu_flux;

        double osc_from_e_w = 0;
        double osc_from_mu_w = 0;
        double final_oscillated_w = 0;

        if(_osc_calculator && _earth){
            std::tie(osc_from_e_w, osc_from_mu_w) = GetOscWeight(Enu, costh, (Flavour)data.nuPDG);
            final_oscillated_w = xsec_w*(osc_from_e_w*nuE_w + osc_from_mu_w*nuMu_w);
        }

        data.genie_weight = data.weight;
        data.flux_nue = nuE_flux;
        data.flux_numu = nuMu_flux;
        data.xsec = xsec_w;
        data.nue_w = nuE_w;
        data.numu_w = nuMu_w;
        data.osc_from_e_w = osc_from_e_w;
        data.osc_from_mu_w = osc_from_mu_w;
        data.final_oscillated_w = final_oscillated_w;

        _wrt.Fill(data);
    }
}

template <typename T>
void Calculator<T>::SetOscCalculator(osc::PMNS *pmns, osc::EarthModel *earth, float prodh, float rdet){
    if(!pmns || !earth){
        std::cerr << "Invalid oscillation calculator provided!" << std::endl;
        abort();
    }
    _osc_calculator = pmns;
    _earth = earth;
    _prodh = prodh;
    _rdet = rdet;
}

template <typename T>
std::tuple<double, double> Calculator<T>::GetOscWeight(double E, double costh, Flavour final){
    if(!_osc_calculator){
        return {0, 0};
    }

    int is_anti = (final < 0) ? -1 : 1;

    std::list<double> L;
    std::list<double> Ne;
    _earth->LineProfile(_prodh, costh, _rdet, L, Ne);
    _osc_calculator->Reset();
    _osc_calculator->PropMatter(L, E, Ne, is_anti);

    int final_flavour_osclib = (final*is_anti - 12)/2;
    return {_osc_calculator->P(0, final_flavour_osclib), _osc_calculator->P(1, final_flavour_osclib)};
}

template class Calculator<float>;
