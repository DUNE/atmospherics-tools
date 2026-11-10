#pragma once
#include "Writer.h"

class DirWriter : public Writer
{
private:
    bool _write_sample_id = true;
    void SetupTree();
    void Create(std::string name);
public:
    DirWriter(std::string name, bool write_sample_id = true);
    ~DirWriter();
};
