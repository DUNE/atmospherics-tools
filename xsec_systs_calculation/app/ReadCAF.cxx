// std
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// ROOT
#include "TObjString.h"
#include "TChain.h"
#include "TFile.h"
// systematicstools
#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"
// duneanaobj
#include "duneanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/Flat/FlatRecord.h"

namespace cliopts {
  std::string input_filename = "";
  size_t NMax = std::numeric_limits<size_t>::max();
  size_t NSkip = 0;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-i <ghep.root>   : GENIE TChain descriptor to read events\n"
               "\t                   from. (n.b. quote wildcards).\n"
               "\t-N <NMax>        : Maximum number of events to process.\n"
               "\t-s <NSkip>       : Number of events to skip.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::input_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = systtools::str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-s") {
      cliopts::NSkip = systtools::str2T<size_t>(argv[++opt]);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {

  HandleOpts(argc, argv);
  if (!cliopts::input_filename.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }

  // Input file
  TFile *f_input = new TFile(cliopts::input_filename.c_str());
/*
  // Global tree
  TTree *t_input_global = (TTree *)f_input->Get("globalTree");
  caf::SRGlobal* srglobal = nullptr;
  t_input_global->SetBranchAddress("global", &srglobal);
  t_input_global->GetEntry(0);

  for(auto sp:srglobal->wgts.params){
    printf("- (id, name, nshifts) = (%d, %s, %d)\n", sp.id, sp.name.c_str(), sp.nshifts);
  }
*/
  // CAF tree

  TTree *t_input_caf(nullptr);
  f_input->GetObject("cafmaker/cafTree", t_input_caf);

  if(!t_input_caf){
    std::cerr << "Could not find input tree cafmaker/cafTree" << std::endl;
    throw;
  }

  size_t NEvs = t_input_caf->GetEntries();
  size_t NToRead = std::min(NEvs, cliopts::NMax);
  printf("@@ Number of CAF events = %ld\n", NEvs);

  // Proxy
  caf::StandardRecordProxy* sr = new caf::StandardRecordProxy(t_input_caf, "rec");

  // Loop over CAFTree
  for (size_t cafev_it = cliopts::NSkip; cafev_it < NToRead; ++cafev_it) {

    printf("CAFTree index: %ld\n", cafev_it);

    t_input_caf->GetEntry(cafev_it);

    size_t N_MC = sr->mc.nu.size();
    printf("- Number of MC particles = %ld\n", N_MC);
    for(size_t i_nu=0; i_nu<N_MC; i_nu++){
      const auto& nu = sr->mc.nu[i_nu];
      printf("  - %ld-th neutrino in this CAF event\n", i_nu);
      std::cout << nu.E << std::endl;
    } // END nu loop

  } // END caf event loop


}
