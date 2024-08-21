#include "FluxManager.h"
#include "Reader.h"
#include "DirWriter.h"
#include "Writer.h"
#include "Calculator.h"
#include <argparse/argparse.hpp>
#include <filesystem>

int main(int argc, char const *argv[])
{
    argparse::ArgumentParser parser("weightor");

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
    //TODO: Create a distinction between the numu flux used to generate the sample, and the one people want to use to reweight.
    // They are assumed to be the same for now!

    parser.add_argument("--numubar")
        .default_value("honda_2d_homestake_2015_numubar.root")
        .help("Specifies the name of the numubar flux inside the flux directory");

    parser.add_argument("-e", "--exposure")
        .default_value(400.f)
        .scan<'g', float>()
        .help("Target exposure in kton.yr");
    
    parser.add_argument("-m", "--mass")
        .default_value(2.5f) //TODO: Pick the exact default value that matches the production
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

    std::filesystem::path fluxdir(parser.get<std::string>("--fluxdir"));
    std::filesystem::path nue_file(parser.get<std::string>("--nue"));
    std::filesystem::path nuebar_file(parser.get<std::string>("--nuebar"));
    std::filesystem::path numu_file(parser.get<std::string>("--numu"));
    std::filesystem::path numubar_file(parser.get<std::string>("--numubar"));

    nue_file = fluxdir / nue_file;
    nuebar_file = fluxdir / nuebar_file;
    numu_file = fluxdir / numu_file;
    numubar_file = fluxdir / numubar_file;

    std::map<Flavour, std::string> fluxes = {
        {Flavour::NuE, nue_file},
        {Flavour::NuMu, numu_file},
        {Flavour::NuEBar, nuebar_file},
        {Flavour::NuMuBar, numubar_file},
    };

    std::string ifilename(parser.get<std::string>("-i"));
    std::string ofilename(parser.get<std::string>("-o"));

    float det_mass = parser.get<float>("-m"); // kton
    float target_exposure = parser.get<float>("-e"); // kton.yr
    int year = 3600*24*365; // seconds/yr

    float exposure_scaling = target_exposure*year/det_mass; // seconds

    float prodh = parser.get<float>("--prodh");
    float deth = parser.get<float>("--deth");
    std::string earth_model = parser.get<std::string>("earthmodel");
    std::vector<float> pmns_params = parser.get<std::vector<float>>("--pmns");

    FluxManager manager(fluxes);
    Reader<float> reader(ifilename, "cafmaker");

    DirWriter writer(ofilename);

    Calculator<float> calc(manager, reader, writer, exposure_scaling);

    osc::EarthModel earth(earth_model.c_str(), 0.02);
    osc::PMNS pmns(pmns_params[0],
        pmns_params[1],
        pmns_params[2],
        pmns_params[3],
        pmns_params[4],
        pmns_params[5]);

    calc.SetOscCalculator(&pmns, &earth, prodh, deth);
    calc.Process();

    writer.GetFile()->cd();
    writer.AddTree(reader.GetTree()->CloneTree(-1, "fast"));
    writer.AddTree(reader.GetGenieTree()->CloneTree(-1, "fast"));
    writer.AddTree(reader.GetGlobalTree()->CloneTree(-1, "fast"));

    writer.Write();

    return 0;
}
