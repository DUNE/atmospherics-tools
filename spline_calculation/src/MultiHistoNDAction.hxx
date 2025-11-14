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

template <typename T, unsigned int NDIM>
class MultiHistoNDHelper : public ROOT::Detail::RDF::RActionImpl<MultiHistoNDHelper<T, NDIM>> {
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

   /// This is a handy, expressive shortcut.
   using THn_t = THnT<T>;
   /// This type is a requirement for every helper.
   using Result_t = std::vector<std::shared_ptr<THn_t>>;
 
private:
   std::vector<std::vector<std::shared_ptr<THn_t>>> fHistos; // one per data processing slot
 
public:
   /// This constructor takes all the parameters necessary to build the THnTs. In addition, it requires the names of
   /// the columns which will be used.
   MultiHistoNDHelper(std::string_view name, std::string_view title, uint nvar, std::array<int, NDIM> nbins, std::array<double, NDIM> xmins,
             std::array<double, NDIM> xmax)
   {
      const auto nSlots = ROOT::IsImplicitMTEnabled() ? ROOT::GetThreadPoolSize() : 1;
      for (auto i : ROOT::TSeqU(nSlots)) {
            fHistos.emplace_back();
        for (uint j = 0; j < nvar; ++j) {
            fHistos[i].emplace_back(std::make_shared<THn_t>(std::string(name).c_str(), std::string(title).c_str(),
                                                      NDIM, nbins.data(), xmins.data(), xmax.data()));
         
        }
      }
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
      // The number of binning variables must match the histogram dimensionality
      static_assert(sizeof...(BinningColumnTypes) == NDIM, "Number of binning columns must match histogram NDIM.");

      // Pack binning values into a tuple for easier handling
      auto binning_tuple = std::make_tuple(binning_values...);

      // Create and fill the array of doubles for THnT::Fill
      std::array<double, NDIM> valuesArr;
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