#include "Oscillogram.h"

using namespace OscProb;

Oscillogram::Oscillogram(const std::vector<double> &Ebins, const std::vector<double> &Czbins, const std::string &premPath) :
_Ebins(Ebins),
_Czbins(Czbins),
_prem(PremModel(premPath))
{
}

std::vector<double> Oscillogram::Compute(OscPars pars){
    _pmns.SetDm(2, pars.dm21);
    _pmns.SetDm(3, pars.dm31);
    _pmns.SetAngle(1, 2, pars.th12);
    _pmns.SetAngle(1, 3, pars.th13);
    _pmns.SetAngle(2, 3, pars.th23);
    _pmns.SetDelta(1, 3, pars.dcp);

    _pmns.SetAvgProbPrec(0.01); //Setting the precision for the averaging


    int N_Ebins = _Ebins.size() - 1;
    int N_Czbins = _Czbins.size() - 1;
    int osc_size = N_Ebins*N_Czbins*2*3*2;
    std::vector<double> oscillogram(osc_size);

    for(uint i = 0; i < N_Czbins; i++){
        double cosT = 0.5*(_Czbins[i] + _Czbins[i + 1]);
        _prem.FillPath(cosT);
        _pmns.SetPath(_prem.GetNuPath());

        for(uint j = 0; j < N_Ebins; j++){
            double E = 0.5*(_Ebins[j] + _Ebins[j + 1]);
            double dE = _Ebins[j + 1] - _Ebins[j];

            int base_idx = 2*3*2*(i*N_Ebins + j);

            _pmns.SetIsNuBar(false);
            const std::vector<std::vector<double>> &mat = _pmns.AvgProbMatrix(2, 3, E, dE);
            _pmns.SetIsNuBar(true);
            const std::vector<std::vector<double>> &mat_anti = _pmns.AvgProbMatrix(2, 3, E, dE);

            oscillogram[base_idx] = mat[0][0];
            oscillogram[base_idx + 1] = mat[0][1];
            oscillogram[base_idx + 2] = mat[0][2];
            oscillogram[base_idx + 3] = mat[1][0];
            oscillogram[base_idx + 4] = mat[1][1];
            oscillogram[base_idx + 5] = mat[1][2];
            
            oscillogram[base_idx + 6] = mat_anti[0][0];
            oscillogram[base_idx + 7] = mat_anti[0][1];
            oscillogram[base_idx + 8] = mat_anti[0][2];
            oscillogram[base_idx + 9] = mat_anti[1][0];
            oscillogram[base_idx + 10] = mat_anti[1][1];
            oscillogram[base_idx + 11] = mat_anti[1][2];
        }
    }

    return oscillogram;
}