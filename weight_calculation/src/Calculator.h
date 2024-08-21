#pragma once
#include <string>
#include "FluxManager.h"
#include "Reader.h"
#include "Writer.h"
#include "OscLib/PMNS.h"
#include "OscLib/EarthModel.h"

template<typename T>
class Calculator
{
private:
    const FluxManager& _mgr;
    Reader<T>& _rdr;
    Writer& _wrt;
    float _exposure_scaling;
    float _prodh;
    float _rdet;
    osc::PMNS* _osc_calculator = nullptr;
    osc::EarthModel* _earth = nullptr;

    Data<double> FetchData();
    std::tuple<double, double> GetOscWeight(double E, double costh, Flavour final);

    

public:
    Calculator(const FluxManager& mgr, Reader<T>& rdr, Writer& wrt, float exposure_scaling);
    ~Calculator();
    void Process();
    void SetOscCalculator(osc::PMNS *pmns, osc::EarthModel *earth, float prodh, float rdet);
};
