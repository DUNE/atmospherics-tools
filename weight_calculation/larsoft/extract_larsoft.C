// Extracts the truth information needed by the weightor directly from LArSoft
// (art/ROOT) files and writes it as a flat ntuple that `weightor --flat` can read.
//
// Must run inside a LArSoft/dunesw environment (dictionaries for simb:: and
// sumdata:: classes are needed), e.g.:
//   root -l -b -q 'extract_larsoft.C("filelist.txt", "flat.root")'
// where filelist.txt holds one art file path per line.
//
// Output:
//   events : one entry per art event
//            run, subrun, eid             (art::EventAuxiliary)
//            Ev, NuMomX, NuMomY, NuMomZ   (MCTruth::GetNeutrino().Nu(), GeV, detector frame)
//            nuPDG, isCC, mode            (MCTruth::GetNeutrino())
//            genie_weight                 (GTruth::fweight, i.e. CAF mc.nu[0].genweight)
//   meta   : one entry per art subrun
//            run, subrun, pot             (POTSummary::totpot = GENIE exposure of that subrun)

#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Provenance/EventAuxiliary.h"
#include "canvas/Persistency/Provenance/SubRunAuxiliary.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "TFile.h"
#include "TTree.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void extract_larsoft(const char *filelist, const char *outname,
                     const char *module = "generator", const char *process = "GenieGen")
{
    std::vector<std::string> files;
    {
        std::ifstream in(filelist);
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line[0] != '#') files.push_back(line);
        }
    }
    if (files.empty()) {
        std::cerr << "No input files found in " << filelist << std::endl;
        exit(1);
    }

    const std::string mctBranch = std::string("simb::MCTruths_") + module + "__" + process + ".";
    const std::string gtBranch  = std::string("simb::GTruths_") + module + "__" + process + ".";
    const std::string potBranch = std::string("sumdata::POTSummary_") + module + "__" + process + ".";

    TFile fout(outname, "RECREATE");
    TTree tev("events", "Truth information extracted from LArSoft files");
    TTree tmeta("meta", "Per-subrun exposure");

    int run, subrun, eid, nuPDG, isCC, mode;
    double Ev, NuMomX, NuMomY, NuMomZ, genie_weight, pot;
    tev.Branch("run", &run);
    tev.Branch("subrun", &subrun);
    tev.Branch("eid", &eid);
    tev.Branch("Ev", &Ev);
    tev.Branch("NuMomX", &NuMomX);
    tev.Branch("NuMomY", &NuMomY);
    tev.Branch("NuMomZ", &NuMomZ);
    tev.Branch("nuPDG", &nuPDG);
    tev.Branch("isCC", &isCC);
    tev.Branch("mode", &mode);
    tev.Branch("genie_weight", &genie_weight);
    tmeta.Branch("run", &run);
    tmeta.Branch("subrun", &subrun);
    tmeta.Branch("pot", &pot);

    long nev_tot = 0;
    double pot_tot = 0;

    for (const auto &fname : files) {
        TFile *f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "ERROR: cannot open " << fname << std::endl;
            exit(1);
        }
        TTree *t = nullptr, *s = nullptr;
        f->GetObject("Events", t);
        f->GetObject("SubRuns", s);
        if (!t || !s || !t->GetBranch(mctBranch.c_str()) || !t->GetBranch(gtBranch.c_str()) ||
            !s->GetBranch(potBranch.c_str())) {
            std::cerr << "ERROR: missing Events/SubRuns trees or " << module << "/" << process
                      << " products in " << fname << std::endl;
            exit(1);
        }

        // Only read the branches we need: the reco files are large.
        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("EventAuxiliary*", 1);
        t->SetBranchStatus((mctBranch + "*").c_str(), 1);
        t->SetBranchStatus((gtBranch + "*").c_str(), 1);

        art::EventAuxiliary *aux = nullptr;
        art::Wrapper<std::vector<simb::MCTruth>> *mct = nullptr;
        art::Wrapper<std::vector<simb::GTruth>> *gt = nullptr;
        t->SetBranchAddress("EventAuxiliary", &aux);
        t->SetBranchAddress(mctBranch.c_str(), &mct);
        t->SetBranchAddress(gtBranch.c_str(), &gt);

        const long n = t->GetEntries();
        for (long i = 0; i < n; i++) {
            t->GetEntry(i);
            if (!mct->isPresent() || !gt->isPresent() || mct->product()->size() != 1 ||
                gt->product()->size() != 1) {
                std::cerr << "ERROR: expected exactly one MCTruth and one GTruth in " << fname
                          << " entry " << i << std::endl;
                exit(1);
            }
            const simb::MCNeutrino &mcnu = mct->product()->at(0).GetNeutrino();
            const simb::MCParticle &nu = mcnu.Nu();
            run = aux->run();
            subrun = aux->subRun();
            eid = aux->event();
            Ev = nu.E();
            NuMomX = nu.Px();
            NuMomY = nu.Py();
            NuMomZ = nu.Pz();
            nuPDG = nu.PdgCode();
            isCC = (mcnu.CCNC() == simb::kCC);
            mode = mcnu.Mode();
            genie_weight = gt->product()->at(0).fweight;
            tev.Fill();
        }
        nev_tot += n;

        art::SubRunAuxiliary *saux = nullptr;
        art::Wrapper<sumdata::POTSummary> *potw = nullptr;
        s->SetBranchAddress("SubRunAuxiliary", &saux);
        s->SetBranchAddress(potBranch.c_str(), &potw);
        double pot_file = 0;
        for (long i = 0; i < s->GetEntries(); i++) {
            s->GetEntry(i);
            if (!potw->isPresent()) {
                std::cerr << "ERROR: POTSummary not present in " << fname << " subrun entry " << i << std::endl;
                exit(1);
            }
            run = saux->run();
            subrun = saux->subRun();
            pot = potw->product()->totpot;
            pot_file += pot;
            tmeta.Fill();
        }
        pot_tot += pot_file;
        std::cout << fname << ": " << n << " events, exposure " << pot_file << std::endl;

        // Objects allocated by SetBranchAddress are owned (and deleted) by the trees.
        f->Close();
        delete f;
    }

    fout.cd();
    tev.Write();
    tmeta.Write();
    fout.Close();
    std::cout << "Wrote " << nev_tot << " events from " << files.size() << " files, total exposure "
              << pot_tot << " to " << outname << std::endl;
}
