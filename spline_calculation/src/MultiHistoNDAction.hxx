#ifndef MULTI_HISTO_ND_ACTION_HXX
#define MULTI_HISTO_ND_ACTION_HXX

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>
#include <ROOT/RDF/RDisplay.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDF/RActionImpl.hxx> // <-- The required base class
#include <ROOT/RDF/HistoModels.hxx>
#include <TString.h> // For Form

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <iostream>


class MultiHistoNDHelper : public ROOT::Detail::RDF::RActionImpl<MultiHistoNDHelper> {
public:
   // C++17 fold expression to fill the array from a tuple
   template <typename Tuple, std::size_t... Is>
   static void fill_array_from_tuple_impl(const Tuple& t, std::array<double, sizeof...(Is)>& arr, std::index_sequence<Is...>) {
      // Use a fold expression to unpack the tuple into the array
      ((arr[Is] = static_cast<double>(std::get<Is>(t))), ...);
   }

   /// This type is a requirement for every helper.
   using Result_t = std::vector<std::shared_ptr<THnT<double>>>;
 
private:
   std::vector<std::vector<std::shared_ptr<THnT<double>>>> fHistos; // one per data processing slot
 
public:
   /// This constructor takes all the parameters necessary to build the THnTs. In addition, it requires the names of
   /// the columns which will be used.
   MultiHistoNDHelper(const ROOT::RDF::THnDModel &model, unsigned int nHistos)
   {
      const auto nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
      fHistos.resize(nSlots);
      for (auto i : ROOT::TSeqU(nSlots)) {
        for (unsigned int j = 0; j < nHistos; ++j) {
            fHistos[i].emplace_back(model.GetHistogram());
        }
      }
   }
   MultiHistoNDHelper(MultiHistoNDHelper &&) = default;
   MultiHistoNDHelper(const MultiHistoNDHelper &) = delete;
   std::shared_ptr<Result_t> GetResultPtr() const { return std::make_shared<Result_t>(fHistos[0]); }
   void Initialize() {}
   void InitTask(TTreeReader *, unsigned int) {}
   /// This is a method executed at every entry
   template <typename WeightColumnType, typename BinningColumnType, typename EventWeightColumnType>
   void Exec(unsigned int slot, const WeightColumnType& weights, BinningColumnType bin_id, EventWeightColumnType event_weight)
   { 
      const size_t nIterations = std::min(fHistos[slot].size(), static_cast<size_t>(weights.size()));
      if (nIterations == 0) {
         return;
      }

      // Get the global bin index once for this event.
      // This assumes all histograms for this slot have the same binning,
      // and we are explicitly not calculating errors.
      // const Long64_t bin = fHistos[slot][0]->GetBin(valuesArr.data());

      for (size_t i = 0; i < nIterations; ++i) {
         const double weight = static_cast<double>(weights[i]*event_weight);
         fHistos[slot][i]->AddBinContent(bin_id, weight);
      }
      return;
   }
   /// This method is called at the end of the event loop. It is used to merge all the internal THnTs which
   /// were used in each of the data processing slots.
   void Finalize()
   {
      // fHistos[0] is the main result, we merge from all other slots into it.
      for (auto slot : ROOT::TSeqU(1, fHistos.size())) {
         for (size_t i = 0; i < fHistos[0].size(); ++i) {
            fHistos[0][i]->Add(fHistos[slot][i].get());
         }
      }
   }
 
   std::string GetActionName() const {
      return "MultiHistoNDHelper";
   }
};


#endif // MULTI_HISTO_ND_ACTION_HXX