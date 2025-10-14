#include "FluxManager.h"
#include "Reader.h"
#include "SRWriter.h"
#include "Calculator.h"
#include "progressbar.hpp"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include <argparse/argparse.hpp>

enum Sel{
    SelNuE,
    SelNuMu,
    SelNC,
    SelNuEBar,
    SelNuMuBar
};

int main(int argc, char const *argv[])
{
    //Splits the provided input file into 2 new output files by applying some Michel-e tagging efficiency

    argparse::ArgumentParser parser("weightor");

    parser.add_argument("-i", "--input")
        .required()
        .help("Input file to process. Must have been produced by weightor before.");

    try {
        parser.parse_args(argc, argv);
    }
    catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string ifilename(parser.get<std::string>("-i"));

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

    TTree *genie_tree = nullptr;
    ifile->GetObject("genieEvt", genie_tree);
    if(genie_tree == nullptr){
        std::cout << "No tree named genieEvt in the input file " << ifilename << std::endl;
        abort();
    }

    genie::NtpMCEventRecord *genie_record = new genie::NtpMCEventRecord();
    genie_tree->SetBranchAddress("genie_record", &genie_record);

    caf::StandardRecord *sr = nullptr;
    caf_tree->SetBranchAddress("rec", &sr);
    
    std::map<Sel, SRWriter*> channels;

    //Get basename of ifilename
    std::string basename = ifilename.substr(0, ifilename.find_last_of(".root")) + "_";

    for(Sel sel : {Sel::SelNuMu, Sel::SelNuMuBar}){
        std::string sel_str;
        switch (sel)
        {
        case Sel::SelNuMu:
            sel_str = "michele_numuselec";
            break;
        case Sel::SelNuMuBar:
            sel_str = "michele_numubarselec";
            break;
        }

        std::string fname = basename + sel_str + ".root";
        SRWriter *ch_writer = new SRWriter(fname, sr);
        ch_writer->SetGenieTree(genie_tree->CloneTree(0));
        channels.insert({sel, ch_writer});
    } 

    uint nentries = caf_tree->GetEntries();
    progressbar bar(nentries);
    bar.set_todo_char(" ");
    bar.set_done_char("█");

    for(uint i = 0; i < nentries; i++){
        caf_tree->GetEntry(i);
        weights_tree->GetEntry(i);
        genie_tree->GetEntry(i);
        bar.update();

        Flavour ofl = Flavour(sr->mc.nu[0].pdg);

        Sel sel = Sel::SelNuMuBar;
        //You need to check here to see what you want to do with wrong ID events, meaning for instance nue events reco as numu, where should they go?
        if(ofl > 0){ //If we have a neutrino, we want to apply the efficiency
            float rand = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
            if(rand > 0.82){
                sel = Sel::SelNuMu;
            }
        }

        SRWriter* channel_writer = channels[sel]; //Getting the corresponding picked channel
        channel_writer->Fill(); //Filling the event in the selected channel
        channel_writer->GetGenieTree()->Fill();

        delete genie_record->event;
    }

    std::cout << std::endl;

    for(auto channel : channels){
        // channel.second->AddNorm();
        // channel.second->AddTree(reader.GetGlobalTree());
        channel.second->Write();
        channel.second->WriteFile();
        channel.second->GetFile()->Close();
        channel.second->_file = nullptr;
    }

    ifile->Close();

    return 0;
}
