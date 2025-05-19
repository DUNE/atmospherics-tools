#include "Reader.h"
#include <argparse/argparse.hpp>
#include <filesystem>

using FLOAT_T = float;

int main(int argc, char const *argv[])
{
    argparse::ArgumentParser parser("ComparisonScript");

    parser.add_argument("-i", "--input")
        .required()
        .help("Input file to process.");

    parser.add_argument("-o", "--output")
        .required()
        .help("Output file to save.");

    parser.add_argument("--fluxdir")
        .default_value("/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/TaskForce_Flux/atmos/Honda_interp/")
        .help("Location of the flux files");

    parser.add_argument("--nue")
        .default_value("honda_2d_homestake_2015_nue.root")
        .help("Specifies the name of the nue flux inside the flux directory");

    parser.add_argument("--nuebar")
        .default_value("honda_2d_homestake_2015_nuebar.root")
        .help("Specifies the name of the nuebar flux inside the flux directory");

    parser.add_argument("--numu")
        .default_value("honda_2d_homestake_2015_numu.root")
        .help("Specifies the name of the numu flux inside the flux directory");

    parser.add_argument("--numubar")
        .default_value("honda_2d_homestake_2015_numubar.root")
        .help("Specifies the name of the numubar flux inside the flux directory");

    parser.add_argument("--ref_flux")
        .default_value("/cvmfs/dune.osgstorage.org/pnfs/fnal.gov/usr/dune/persistent/stash/TaskForce_Flux/atmos/Honda_interp/honda_2d_homestake_2015_numu.root")
        .help("Specifies the FULL PATH to the reference flux used to weight the sample at generation. THE DEFAULT VALUE SHOULD BE FINE.");

    parser.add_argument("-e", "--exposure")
        .default_value(400.f)
        .scan<'g', float>()
        .help("Target exposure in kton.yr");
    
    parser.add_argument("-m", "--mass")
        .default_value(1.71958f) //Default value for 1x2x6 production
        .scan<'g', float>()
        .help("Detector mass in kton USED FOR THE SIMULATION.");

    parser.add_argument("--prodh")
        .default_value(20.f)
        .scan<'g', float>()
        .help("Production height of the neutrinos in the atmosphere (km).");

    parser.add_argument("--deth")
        .default_value(-0.3f)
        .scan<'g', float>()
        .help("Detector height above ground (should be negative for burried detectors) (km).");

    parser.add_argument("--earthmodel")
        .default_value("PREM")
        .help("Earth model to use for the oscillation calculation. Avail: [\"PREM\", \"STACEY\"]");

    parser.add_argument("--pmns")
        .nargs(6)
        .default_value(std::vector<float>{0.583, 0.737, 0.150, 4.05, 7.41e-5, 2.51e-3})
        .scan<'g', float>()
        .help("List of 6 oscillation parameters to use under the form: th12 th23 th13 dcp dms12 dms23 (rad, eV^2)\nDefault values are taken from NuFIT 5.2 w. SK atm.: 0.583 0.737 0.150 4.05 7.41e-5 2.51e-3");

    try {
        parser.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string ifilename(parser.get<std::string>("-i"));

    Reader<FLOAT_T> reader(ifilename, "cafmaker");

    for (int i=0;i<reader.GetNentries();i++) {
      reader.GetEntry(i);
      const Data<FLOAT_T> Data = reader.GetData();
      std::cout << i << " " << Data.mode << " " << Data.erec << " " << Data.NuMomX << std::endl;
    }
    
    return 0;
}
