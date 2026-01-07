#include <iostream>
#include <string>

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Conventions/Constants.h" //for calculating event kinematics

Int_t fIsCC, fPdg;
Float_t fLep_Px, fLep_Py, fLep_Pz, fLep_E;
Float_t fNu_Px, fNu_Py, fNu_Pz, fNu_E, fR_nucleus;
Float_t fHitNuc_Px, fHitNuc_Py, fHitNuc_Pz, fHitNuc_E, fRemovalE, Q2, y, v, q0, q3, fDecayAngMEC;
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
    otree->Branch("DecayAngMEC", &fDecayAngMEC, "DecayAngMEC/F");

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

        bool is_mec = inter->ProcInfo().IsMEC();
        if(is_mec){
            const int recoil_nucleon_cluster_pos = 5;
            genie::GHepParticle* nucleon_cluster = fEventRecord->event->Particle( recoil_nucleon_cluster_pos );
            if(!nucleon_cluster){
                throw("Expected to find a nucleon cluster particle in MEC event");
            }
            int cluster_pdg = nucleon_cluster->Pdg();
            if ( !genie::pdg::Is2NucleonCluster(cluster_pdg) ) {
                throw("Expected nucleon cluster to have a 2-nucleon-cluster PDG code");
            }
            int first = nucleon_cluster->FirstDaughter();
            int last = nucleon_cluster->LastDaughter();
            if ( !nucleon_cluster->HasDaughters() || (last - first) != 1 ) {
                throw("Invalid number of daughters for a two-nucleon cluster");
            }
            genie::GHepParticle* N1 = fEventRecord->event->Particle( first );
            genie::GHepParticle* N2 = fEventRecord->event->Particle( last );
            if ( !genie::pdg::IsNucleon(N1->Pdg()) || !genie::pdg::IsNucleon(N2->Pdg()) ) {
                throw("Expected daughters of nucleon cluster to be nucleons");
            }
            // Get the 4-momenta of the two outgoing nucleons
            TLorentzVector p4N1 = *N1->P4();
            TLorentzVector p4N2 = *N2->P4();

            // Boost the 4-momenta of the two nucleons from the lab frame to their
            // CM frame (which is also the rest frame of the recoiling nucleon cluster)
            TLorentzVector p4Cluster = p4N1 + p4N2;
            TVector3 boostToCM = -p4Cluster.BoostVector();

            p4N1.Boost( boostToCM );
            p4N2.Boost( boostToCM );
            q.Boost( boostToCM );
            // Use the 3-momentum transfer in the two-nucleon CM frame as the reference
            // z-axis for the altered angular distribution
            TVector3 q3 = q.Vect().Unit();

            // Determine a rotation axis and angle that will cause the 3-momentum to
            // point along the +z direction
            TVector3 zvec(0., 0., 1.);
            TVector3 rot = ( q3.Cross(zvec) ).Unit();
            double angle = zvec.Angle( q3 );

            // Handle the edge case where q3 is along -z, so the
            // cross product above vanishes
            if ( q3.Perp() == 0. && q3.Z() < 0. ) {
                rot = TVector3(0., 1., 0.);
                angle = genie::constants::kPi;
            }

            // If the rotation vector is non-null (within numerical precision) then
            // rotate the CM frame 3-momentum of nucleon #1 into a frame where q3 points along +z
            TVector3 p3N1 = p4N1.Vect();
            if ( rot.Mag() >= 1e-10 ) {
                p3N1.Rotate(angle, rot);
            }

            // We now have what we need. Compute the emission angles for nucleon #1 relative to the
            // 3-momentum transfer in the rest frame of the recoiling nucleon cluster.
            fDecayAngMEC = p3N1.Theta();
        } else {
            fDecayAngMEC = -9999.;
        }

        otree->Fill();

        delete fEventRecord->event;
    }
    
    ifile->Close();
    ofile->cd();
    otree->Write();
    ofile->Close();
}