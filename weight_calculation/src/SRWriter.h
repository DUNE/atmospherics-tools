#pragma once
#include "Writer.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"

class SRWriter : public Writer
{
private:
    void SetupTree();
    void Create(std::string name);
    caf::StandardRecord *_sr = nullptr;
public:
    SRWriter(std::string name, caf::StandardRecord *sr);
    ~SRWriter();
    void Fill();
};
