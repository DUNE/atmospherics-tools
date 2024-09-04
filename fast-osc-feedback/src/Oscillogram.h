#pragma once
#include <vector>
#include "PremModel.h"
#include "PMNS_Fast.h"

struct OscPars
{
    double dm21;
    double dm31;
    double th12;
    double th13;
    double th23;
    double dcp;
};


class Oscillogram {
    public:
        Oscillogram(const std::vector<double> &Ebins, const std::vector<double> &Czbins);

        std::vector<double> Compute(OscPars pars);

    private:
        const std::vector<double> _Ebins;
        const std::vector<double> _Czbins;
        OscProb::PremModel _prem;
        OscProb::PMNS_Fast _pmns;

};