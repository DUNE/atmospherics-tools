#include "TFile.h"
#include "TTree.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
  TFile *f = TFile::Open("/exp/dune/data/users/pgranger/xsec_systematics_atmospherics_sbn.root");
  if (!f || f->IsZombie()) {
    std::cerr << "Could not open file" << std::endl;
    return 1;
  }
  
  TTree *t_weights = (TTree*)f->Get("SystWeights");
  TTree *t_genie = (TTree*)f->Get("genieEvt");
  
  if (!t_weights || !t_genie) {
    std::cerr << "Could not find trees" << std::endl;
    return 1;
  }
  
  long long nentries = t_weights->GetEntries();
  std::cout << "Total entries: " << nentries << std::endl;
  
  // Set up weight branch
  double FrAbs_N[7];
  t_weights->SetBranchAddress("GENIEReWeight_ICARUS_v2_multisigma_FrAbs_N", FrAbs_N);
  
  // Set up genie branch
  genie::NtpMCEventRecord *GenieNtpl = nullptr;
  t_genie->SetBranchAddress("genie_record", &GenieNtpl);
  
  int nan_count = 0;
  int non_nan_count = 0;
  
  std::cout << "\n--- SCANNING FIRST 50 EVENTS ---" << std::endl;
  for (long long i = 0; i < nentries; ++i) {
    t_weights->GetEntry(i);
    t_genie->GetEntry(i);
    
    bool is_nan = std::isnan(FrAbs_N[0]);
    if (is_nan) nan_count++;
    else non_nan_count++;
    
    if (i < 30) {
      genie::EventRecord *rec = GenieNtpl->event;
      std::cout << "Entry " << std::setw(3) << i 
                << " | is_nan=" << (is_nan ? "YES" : " NO")
                << " | Target=" << rec->Summary()->InitState().Tgt().Pdg()
                << " | Probe=" << rec->Summary()->InitState().ProbePdg()
                << " | Proc=" << rec->Summary()->ProcInfo().AsString()
                << " | NumParticles=" << rec->GetEntries()
                << std::endl;
                
      // If it's a NaN, let's print particles in this event to see what FSI states look like
      if (is_nan && nan_count <= 5) {
        std::cout << "    Particles in NaN event " << i << ":" << std::endl;
        for (int p = 0; p < rec->GetEntries(); ++p) {
          genie::GHepParticle *part = rec->Particle(p);
          std::cout << "      [" << p << "] PDG=" << std::setw(6) << part->Pdg()
                    << " | Status=" << std::setw(2) << part->Status()
                    << " | Name=" << part->Name()
                    << " | Mother1=" << part->FirstMother()
                    << " | Mother2=" << part->LastMother()
                    << std::endl;
        }
      }
    }
  }
  
  std::cout << "\nSummary:" << std::endl;
  std::cout << "  NaN events: " << nan_count << " (" << (100.0 * nan_count / nentries) << "%)" << std::endl;
  std::cout << "  Valid events: " << non_nan_count << " (" << (100.0 * non_nan_count / nentries) << "%)" << std::endl;
  
  f->Close();
  return 0;
}
