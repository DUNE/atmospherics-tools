#include <iostream>
#include <string>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/Conventions/Constants.h" //for calculating event kinematics

Int_t fIsCC, fPdg;
Float_t fLep_Px, fLep_Py, fLep_Pz, fLep_E;
Float_t fNu_Px, fNu_Py, fNu_Pz, fNu_E, fR_nucleus;
Float_t fHitNuc_Px, fHitNuc_Py, fHitNuc_Pz, fHitNuc_E, fRemovalE, Q2, y, v, q0, q3;
std::string fProc;

TTree* setup_otree(){
    TTree *otree = new TTree("genie_dump", "GENIE event tree");
    otree->Branch("isCC", &fIsCC, "isCC/I");
    otree->Branch("pdg", &fPdg, "pdg/I");
    otree->Branch("lep_Px", &fLep_Px, "lep_Px/F");
    otree->Branch("lep_Py", &fLep_Py, "lep_Py/F");
    otree->Branch("lep_Pz", &fLep_Pz, "lep_Pz/F");
    otree->Branch("lep_E", &fLep_E, "lep_E/F");
    otree->Branch("nu_Px", &fNu_Px, "nu_Px/F");
    otree->Branch("nu_Py", &fNu_Py, "nu_Py/F");
    otree->Branch("nu_Pz", &fNu_Pz, "nu_Pz/F");
    otree->Branch("nu_E", &fNu_E, "nu_E/F");
    otree->Branch("hitnuc_Px", &fHitNuc_Px, "hitnuc_Px/F");
    otree->Branch("hitnuc_Py", &fHitNuc_Py, "hitnuc_Py/F");
    otree->Branch("hitnuc_Pz", &fHitNuc_Pz, "hitnuc_Pz/F");
    otree->Branch("hitnuc_E", &fHitNuc_E, "hitnuc_E/F");
    otree->Branch("removalE", &fRemovalE, "removalE/F");
    otree->Branch("R_nucleus", &fR_nucleus, "R_nucleus/F");
    otree->Branch("Q2", &Q2, "Q2/F");
    otree->Branch("y", &y, "y/F");
    otree->Branch("v", &v, "v/F");
    otree->Branch("q0", &q0, "q0/F");
    otree->Branch("q3", &q3, "q3/F");

    otree->Branch("proc", &fProc);
    return otree;
}

void extract_genie(std::string ifilename, std::string ofilename){
    TFile *ifile = new TFile(ifilename.c_str(), "READ");
    TFile *ofile = new TFile(ofilename.c_str(), "RECREATE");

    TTree *otree = setup_otree();
    
    TTree * genie_tree;
    ifile->GetObject("genieEvt", genie_tree);
    if(genie_tree == nullptr){
        std::cerr << "Error: Could not find the genieEvt tree in the input file." << std::endl;
        return;
    }

    genie::NtpMCEventRecord *fEventRecord = nullptr;
    genie_tree->SetBranchAddress("genie_record", &fEventRecord);

    uint nentries = genie_tree->GetEntries();
    uint pct = nentries/100;

    for (int i = 0; i < nentries; i++){
        if(i % pct == 0){
            std::cout << "Processing event " << i << " of " << nentries << std::endl;
        }
        genie_tree->GetEntry(i);
        const genie::Interaction *inter = fEventRecord->event->Summary();
        // get the different components making up the interaction
        const genie::InitialState &initState  = inter->InitState();
        const genie::ProcessInfo  &procInfo   = inter->ProcInfo();

        fIsCC = !procInfo.IsWeakNC();
        fProc = procInfo.AsString();

        const TLorentzVector v4_null;
        genie::GHepParticle* probe = fEventRecord->event->Probe();
        genie::GHepParticle* finallepton = fEventRecord->event->FinalStatePrimaryLepton();
        genie::GHepParticle * hitnucl = fEventRecord->event->HitNucleon();
        const TLorentzVector & k1 = ( probe ? *(probe->P4()) : v4_null );
        const TLorentzVector & k2 = ( finallepton ? *(finallepton->P4()) : v4_null );
        const TLorentzVector & k3 = ( hitnucl ? *(hitnucl->P4()) : v4_null );

        fPdg = probe->Pdg();

        fLep_Px = k2.Px();
        fLep_Py = k2.Py();
        fLep_Pz = k2.Pz();
        fLep_E = k2.E();
        fNu_Px = k1.Px();
        fNu_Py = k1.Py();
        fNu_Pz = k1.Pz();
        fNu_E = k1.E();
        fHitNuc_Px = k3.Px();
        fHitNuc_Py = k3.Py();
        fHitNuc_Pz = k3.Pz();
        fHitNuc_E = k3.E();
        if(hitnucl)
            fRemovalE = hitnucl->RemovalEnergy();
        else
            fRemovalE = 0;

        TLorentzVector q  = k1-k2;

        Q2 = -q.M2();
        v = q.Energy();
        y = v/k1.E();
        q0 = q.Energy();
        q3 = q.P();

        otree->Fill();

        delete fEventRecord->event;
    }
    
    ifile->Close();
    ofile->cd();
    otree->Write();
    ofile->Close();
}