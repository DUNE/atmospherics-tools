// std
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// fhiclcpp
#include "fhiclcpp/ParameterSet.h"
// systematicstools
#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"
// nusystematics
#include "nusystematics/utility/GENIEUtils.hh"
#include "nusystematics/utility/enumclass2int.hh"
#include "nusystematics/utility/response_helper.hh"
// GENIE
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
// ROOT
#include "TObjString.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSpline.h"
// duneanaobj
#include "duneanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "duneanaobj/StandardRecord/StandardRecord.h"
#include "duneanaobj/StandardRecord/SRGlobal.h"
#include "duneanaobj/StandardRecord/Flat/FlatRecord.h"

namespace cliopts {
  std::string fclname = "";
  std::string input_filename = "";
  std::string output_filename = "";
  std::string envvar = "FHICL_FILE_PATH";
  std::string fhicl_key = "generated_systematic_provider_configuration";
  size_t NMax = std::numeric_limits<size_t>::max();
  size_t NSkip = 0;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-c <config.fcl>  : fhicl file to read.\n"
               "\t-k <list key>    : fhicl key to look for parameter headers,\n"
               "\t                   "
               "\"generated_systematic_provider_configuration\"\n"
               "\t                   by default.\n"
               "\t-i <ghep.root>   : GENIE TChain descriptor to read events\n"
               "\t                   from. (n.b. quote wildcards).\n"
               "\t-N <NMax>        : Maximum number of events to process.\n"
               "\t-s <NSkip>       : Number of events to skip.\n"
               "\t-o <out.root>    : File to write validation canvases to.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-c") {
      cliopts::fclname = argv[++opt];
    } else if (std::string(argv[opt]) == "-k") {
      cliopts::fhicl_key = argv[++opt];
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::input_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = systtools::str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-s") {
      cliopts::NSkip = systtools::str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-o") {
      cliopts::output_filename = argv[++opt];
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {

  genie::Messenger::Instance()->SetPrioritiesFromXmlFile("Messenger_laconic.xml"); // quiet mode

  HandleOpts(argc, argv);
  if (!cliopts::fclname.size()) {
    std::cout << "[ERROR]: Expected to be passed a -c option." << std::endl;
    SayUsage(argv);
    return 1;
  }
  if (!cliopts::input_filename.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }

  // Input file
  TFile *f_input = new TFile(cliopts::input_filename.c_str());

  // CAF tree
  TTree *t_input_caftree(nullptr);
  f_input->GetObject("cafmaker/cafTree", t_input_caftree);

  if(!t_input_caftree){
    std::cerr << "Could not find input tree cafTree" << std::endl;
    throw;
  }
  size_t NEvs = t_input_caftree->GetEntries();
  size_t NToRead = std::min(NEvs, cliopts::NMax);
  printf("@@ Number of CAF events = %ld\n", NEvs);
  caf::StandardRecord* sr = nullptr; // we will update this and write it to the output
  t_input_caftree->SetBranchAddress("rec", &sr);

  // GENIE tree
  TTree *t_input_genie(nullptr);
  f_input->GetObject("cafmaker/genieEvt", t_input_genie);

  if(!t_input_genie){
    std::cerr << "Could not find input tree cafmaker/genieEvt" << std::endl;
    throw;
  }

  genie::NtpMCEventRecord *GenieNtpl = nullptr;
  if (t_input_genie->SetBranchAddress("genie_record", &GenieNtpl) != TTree::kMatch) {
    std::cout << "[ERROR]: Failed to set branch address on ghep tree."
              << std::endl;
    return 1;
  }
  size_t NGENIEEvents = t_input_genie->GetEntries();
  printf("@@ Number of genie events = %ld\n", NGENIEEvents);
  
  // Proxy
  caf::StandardRecordProxy* srproxy = new caf::StandardRecordProxy(t_input_caftree, "rec");

  // nusyst response_helper
  nusyst::response_helper resp_helper(cliopts::fclname);

  // SRGlobal
  caf::SRGlobal srglobal = caf::SRGlobal();
  srglobal.wgts.params.clear();
  printf("@@ Writting Header\n");
  for(systtools::paramId_t pid : resp_helper.GetParameters()) {
    systtools::SystParamHeader const &hdr = resp_helper.GetHeader(pid);

    srglobal.wgts.params.emplace_back();
    // Name
    srglobal.wgts.params.back().name = hdr.prettyName;
    // TODO better save paramVariations
    if (hdr.isCorrection) {
      srglobal.wgts.params.back().nshifts = 1;
    } else {
      srglobal.wgts.params.back().nshifts = hdr.paramVariations.size();
    }
    // ParamID
    srglobal.wgts.params.back().id = pid;
  }
  printf("@@ Printing SRGlobal\n");
  for(const auto& sp:srglobal.wgts.params){
    printf("- (id, name, nshifts) = (%d, %s, %d)\n", sp.id, sp.name.c_str(), sp.nshifts);
  }

  // Loop over CAFTree
  for (size_t cafev_it = cliopts::NSkip; cafev_it < NToRead; ++cafev_it) {

    t_input_caftree->GetEntry(cafev_it);

    size_t N_MC = srproxy->mc.nu.size();
    for(size_t i_nu=0; i_nu<N_MC; i_nu++){
      auto& nu = srproxy->mc.nu[i_nu];
      auto mode = nu.mode;
      size_t genieIdx = nu.genieIdx;

      if (genieIdx > t_input_genie->GetEntries()) {
	std::cout << "Was expecting genieIdx < t_input_genie->GetEntries()" << std::endl;
	std::cout << "genieIdx:" << genieIdx << std::endl;
	std::cout << "t_input_genie->GetEntries():" << t_input_genie->GetEntries() << std::endl;
	throw;
      }
      t_input_genie->GetEntry(genieIdx);
      genie::EventRecord const &GenieGHep = *GenieNtpl->event;

      //DB: Force the assigned weight to be equal to exactly 1. so we renormalise the variation weight with respect to the flux weights
      //genie::EventRecord CopyGenieEventRecord(GenieGHep);
      //CopyGenieEventRecord.SetWeight(1.);

      // Evaluate reweights
      systtools::event_unit_response_w_cv_t resp = resp_helper.GetEventVariationAndCVResponse(GenieGHep);

      for(const auto& v: resp){
        const systtools::paramId_t& pid = v.pid;
	systtools::SystParamHeader const &hdr = resp_helper.GetHeader(pid);
	if (hdr.isCorrection) continue;
	
        const double& CVw = v.CV_response;
        const std::vector<double>& ws = v.responses;

	std::cout << "- EventID:" << cafev_it << ", Mode: " << mode << ", ParamID:" << pid << ": RW values = {";
        if (!ws.empty()) {
	  std::cout << ws[0];
          for (size_t i = 1; i < ws.size(); ++i) {
	    std::cout << ", " << ws[i];
          }
        }
	std::cout << "}" << std::endl;

        size_t nPoints = ws.size();

        if (nPoints != hdr.paramVariations.size()) {
	  std::cout << "nPoints:" << nPoints << std::endl;
	  std::cout << "hdr.paramVariations.size():" << hdr.paramVariations.size() << std::endl;
          throw;
        }
        double* XVals = new double[nPoints];
        double* YVals = new double[nPoints];
        for (int iP=0;iP<nPoints;iP++) {
          XVals[iP] = hdr.paramVariations.at(iP);
          YVals[iP] = ws.at(iP);
        }
        for (int iP=0;iP<nPoints;iP++) {
	  std::cout << "\t\tVals: " << iP << " , (X,Y) = (" << XVals[iP] << "," << YVals[iP] << ")" << std::endl;
        }

	TGraph Graph(nPoints,XVals,YVals);
	TSpline3 Spline(Form("Event_%i_Syst_%s",int(cafev_it),hdr.prettyName.c_str()),&Graph);
	
	std::cout << "\t" << Spline.GetTitle() << " : nKnots:" << Spline.GetNp() << std::endl;
	double x, y;
	for (int iP=0;iP<Spline.GetNp();iP++) {
	  Spline.GetKnot(iP,x,y);
	  std::cout << "\t\tKnot: " << iP << " , (X,Y) = (" << x << "," << y << ")" << std::endl;
	}

      } // END resp loop

      // UPDATE RECORD HERE
      // E.g., sr->mc.nu[i_nu].E = <updated value>

    } // END nu loop

    //t_output_caftree->Fill();
  } // END caf event loop

  // Finalize output

  std::cout << "@@ Closing input" << std::endl;
  f_input->Close();
}
