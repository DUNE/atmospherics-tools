#include "FluxManager.h"
#include "Reader.h"
#include "SRWriter.h"
#include "Calculator.h"
#include "progressbar.hpp"
#include <argparse/argparse.hpp>

enum Sel{
    SelNuE,
    SelNuMu,
    SelNC
};

int main(int argc, char const *argv[])
{
    argparse::ArgumentParser parser("weightor");

    parser.add_argument("-i", "--input")
        .required()
        .help("Input file to process. Must have been produced by weightor before.");

    parser.add_argument("-o", "--output")
        .required()
        .help("Output file basename (no extension).");

    parser.add_argument("--cvn_numu")
        .default_value(0.56f)
        .help("CVN threshold to select an event as numu. (applies first)");

    parser.add_argument("--cvn_nue")
        .default_value(0.55f)
        .help("CVN threshold to select an event as nue. (applies second)");

    try {
        parser.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string ifilename(parser.get<std::string>("-i"));
    std::string ofilename(parser.get<std::string>("-o"));

    float cvn_numu = parser.get<float>("--cvn_numu");
    float cvn_nue = parser.get<float>("--cvn_nue");

    std::cout << "Opening file: " << ifilename << std::endl;
    TFile *ifile = TFile::Open(ifilename.c_str(), "READ");
    if(ifile == nullptr){
        std::cout << "Could not open file: " << ifilename << std::endl;
        abort();
    }

    TTree* caf_tree = nullptr;
    ifile->GetObject("cafTree", caf_tree);
    if(caf_tree == nullptr){
        std::cout << "No tree named weights in the input file " << ifilename << std::endl;
        abort();
    }

    TTree* weights_tree = nullptr;
    ifile->GetObject("weights", weights_tree);
    if(weights_tree == nullptr){
        std::cout << "No tree named weights in the input file " << ifilename << std::endl;
        abort();
    }
    double nuE_flux, nuMu_flux, xsec_w;
    weights_tree->SetBranchAddress("flux_nue", &nuE_flux);
    weights_tree->SetBranchAddress("flux_numu", &nuMu_flux);
    weights_tree->SetBranchAddress("xsec", &xsec_w);

    caf::StandardRecord *sr = nullptr;
    caf_tree->SetBranchAddress("rec", &sr);
    
    std::map<std::tuple<Flavour, Flavour, Sel>, SRWriter*> channels;
    
    for(Flavour ifl : {Flavour::NuE, Flavour::NuMu, Flavour::NuEBar, Flavour::NuMuBar}){
        for(Flavour ofl : {Flavour::NuE, Flavour::NuMu, Flavour::NuTau, Flavour::NuEBar, Flavour::NuMuBar, Flavour::NuTauBar}){
            if(ifl*ofl < 0){//nu -> nubar osc impossible
                continue;
            }
            for(Sel sel : {Sel::SelNuE, Sel::SelNuMu, Sel::SelNC}){
                std::string sel_str;
                switch (sel)
                {
                case Sel::SelNuE:
                    sel_str = "nueselec";
                    break;
                case Sel::SelNuMu:
                    sel_str = "numuselec";
                    break;
                case Sel::SelNC:
                    sel_str = "ncselec";
                    break;
                }

                std::string fname = ofilename + flavours[ifl] + "_x_" + flavours[ofl] + "_" + sel_str + ".root";
                SRWriter *ch_writer = new SRWriter(fname, sr);
                channels.insert({{ifl, ofl, sel}, ch_writer});
            } 
        }
    }

    uint nentries = caf_tree->GetEntries();

    progressbar bar(nentries);
    bar.set_todo_char(" ");
    bar.set_done_char("█");

    for(uint i = 0; i < nentries; i++){
        caf_tree->GetEntry(i);
        weights_tree->GetEntry(i);
        bar.update();

        Flavour ofl = Flavour(sr->mc.nu[0].pdg);
        std::vector<Flavour> ifls;

        if(ofl < 0){
            ifls = {Flavour::NuEBar, Flavour::NuMuBar};
        }
        else{
            ifls = {Flavour::NuE, Flavour::NuMu};
        }

        Sel sel = Sel::SelNC; // By defaukt an event is NC if it does not pass any cvn threshold

        if(sr->common.ixn.pandora.size() != 1){
            continue;
        }

        if(sr->common.ixn.pandora[0].nuhyp.cvn.numu > cvn_numu){
            sel = Sel::SelNuMu;
        }
        else if(sr->common.ixn.pandora[0].nuhyp.cvn.nue > cvn_nue){
            sel = Sel::SelNuE;
        }

        SRWriter* channel_writer = channels[{ifls[0], ofl, sel}];
        sr->mc.nu[0].genweight = xsec_w*nuE_flux;
        channel_writer->Fill();

        channel_writer = channels[{ifls[1], ofl, sel}];
        sr->mc.nu[0].genweight = xsec_w*nuMu_flux;
        channel_writer->Fill();
    }

    std::cout << std::endl;

    for(auto channel : channels){
        // channel.second->AddNorm();
        // channel.second->AddTree(reader.GetGlobalTree());
        channel.second->Write();
        channel.second->WriteFile();
    }

    return 0;
}
