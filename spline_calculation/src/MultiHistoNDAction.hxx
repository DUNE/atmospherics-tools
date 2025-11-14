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
   // Recursive helper to fill the values array from a tuple
   template <std::size_t I = 0, typename... Tp>
   static typename std::enable_if<I == sizeof...(Tp), void>::type
   fill_array_from_tuple(const std::tuple<Tp...>&, std::array<double, sizeof...(Tp)>&) {
      // Base case: end of recursion
   }

   template <std::size_t I = 0, typename... Tp>
   static typename std::enable_if<I < sizeof...(Tp), void>::type
   fill_array_from_tuple(const std::tuple<Tp...>& t, std::array<double, sizeof...(Tp)>& arr) {
      arr[I] = static_cast<double>(std::get<I>(t));
      fill_array_from_tuple<I + 1, Tp...>(t, arr);
   }

   /// This type is a requirement for every helper.
   using Result_t = std::vector<std::shared_ptr<THnT<double>>>;
 
private:
   std::vector<std::vector<std::shared_ptr<THnT<double>>>> fHistos; // one per data processing slot
   ROOT::RDF::ColumnNames_t fcolumnList;
 
public:
   /// This constructor takes all the parameters necessary to build the THnTs. In addition, it requires the names of
   /// the columns which will be used.
   MultiHistoNDHelper(const ROOT::RDF::THnDModel &model, unsigned int nHistos, const ROOT::RDF::ColumnNames_t &columnList)
   {
      const auto nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
      fHistos.resize(nSlots);
      for (auto i : ROOT::TSeqU(nSlots)) {
        for (unsigned int j = 0; j < nHistos; ++j) {
            fHistos[i].emplace_back(model.GetHistogram());
        }
      }
      fcolumnList = columnList;
   }
   MultiHistoNDHelper(MultiHistoNDHelper &&) = default;
   MultiHistoNDHelper(const MultiHistoNDHelper &) = delete;
   std::shared_ptr<Result_t> GetResultPtr() const { return std::make_shared<Result_t>(fHistos[0]); }
   void Initialize() {}
   void InitTask(TTreeReader *, unsigned int) {}
   /// This is a method executed at every entry
   template <typename WeightColumnType, typename... BinningColumnTypes>
   void Exec(unsigned int slot, const WeightColumnType& weights, BinningColumnTypes... binning_values)
   {
      // Deduce the dimensionality from the number of binning columns provided.
      constexpr size_t NDIM = sizeof...(BinningColumnTypes);

      // Pack binning values into a tuple for easier handling
      auto binning_tuple = std::make_tuple(binning_values...);

      // Create and fill the array of doubles for THn_t::Fill
      std::array<double, sizeof...(BinningColumnTypes)> valuesArr;
      fill_array_from_tuple(binning_tuple, valuesArr);

      // Ensure we don't try to fill more histograms than we have
      const auto nHistos = fHistos[slot].size();
      const auto nWeights = weights.size();

      for (size_t i = 0; i < nHistos && i < nWeights; ++i) {
         fHistos[slot][i]->Fill(valuesArr.data(), static_cast<double>(weights[i]));
      }
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
 
   std::string GetActionName(){
      return "MultiHistoNDHelper";
   }
};


#endif // MULTI_HISTO_ND_ACTION_HXX