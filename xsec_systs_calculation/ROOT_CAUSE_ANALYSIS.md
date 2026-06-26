# Root-cause analysis — why some dials are inert / pathological

Investigation of the variation-production setup (`UpdateReweight.cxx`,
`fcl/systs_atmospherics_v2.fcl`, `genie_config/config/*`, job `.out` logs) explaining the
three issues seen in the wsbn output (`merged_new_systs_wsbn_wheader.root`). Generation tune:
**AR23_20i_00_000**; reweighting GENIE: **v3_06_02e** (standard larsoft) per `setup.sh`.

## Issue 1 — 27 dials are identically unit. Three distinct causes:

### 1a. The 16 SBN FSI dials — wrong GENIE build (need the SBN fork)
`FrG4{,LoE,M1E,M2E,HiE}_N`, `FrINCL{,LoE,M1E,M2E,HiE}_N`, `MFP{Lo,M1,M2,Hi}E_N`,
`FrKin_PiPro{Fix,Bias}_N` (the `GENIEReWeight ... SBN_v3` provider) are **SBN-specific
extensions** — Geant4/INCL-cascade and energy-binned MFP / kinematic FSI dials. They are
implemented in the **SBN GENIE fork (`v3_06_02_sbn1/sbn2`)**, but `setup.sh` sets up
**standard `genie v3_06_02e`**, which doesn't implement them → GReWeight returns 1.0.

Evidence: `genie_config/config/GSystUncertaintyTable.xml` has 1σ entries **only** for the
standard hA dials (`FrAbs/FrCEx/FrInel/FrPiProd/MFP` for π and N) — exactly the dials that
**do** vary in the output — and **none** for `FrG4/FrINCL/MFP*E/FrKin`. The fcl was clearly
ported from an SBN config (the CCQE block even references
`reweight_data_..._v3_06_02_sbn2_...`), but the build links the non-SBN GENIE.
→ **Fix:** build/run nusystematics against `genie v3_06_02_sbn2` (+ matching `genie_xsec`),
or remove these dials from the fcl for a v3_06_02e production.

### 1b. The MEC model-morph dials — alternative MEC model fails to configure
`XSecShape_CCMEC_Martini`, `XSecShape_CCMEC_Empirical`, `EnergyDependence_CCMEC`,
`FracDelta_CCMEC`, `DecayAng2MEC`, `DecayAngMECVariationResponse` need a *live alternative MEC
cross-section model* to reweight to. The job log shows it failing:
```
FATAL  ... Key: HadronTensorAlg does not exist in pools from algorithm :
       genie::MartiniEricsonChanfrayMarteauMECPXSec2024/Default
WARN   ... No Configuration available for genie::MartiniEricsonChanfrayMarteauMECPXSec2024/Default
```
Root cause: **`genie_config/config/reweight_master_config.xml` is incomplete** — it registers
only 3 algorithms:
```
genie::rew::ObservableMuonMomentum, genie::rew::ObservablePMuEnu, genie::rew::GSystUncertaintyTable
```
It omits the Martini/Empirical MEC PXSec and their hadron-tensor models. So during reweighting
the Martini MEC PXSec can't resolve its `HadronTensorAlg` (`genie::MartiniMECHadronTensorModel/
Default`) from the reweight ConfigPool, the model can't be built, and the dial silently → 1.0.
(The full `master_config.xml` *does* register these algorithms; the reweight pool doesn't pull
them in.) The "shape/fraction" MEC dials that don't need a model swap (`XSecShape_CCMEC`,
`FracPN_CCMEC`, `DecayAngMEC`) work fine.
→ **Fix:** add the MEC-model + hadron-tensor `<config>` lines (Martini*, EmpiricalMEC*,
SuSAv2MEC*, Nieves*) from `master_config.xml` into `reweight_master_config.xml`.

### 1c. Expected / benign (not bugs)
- `CCQEXSecCorr`, `Theta_Delta2Npi`, `VecFFCCQEshape`: corrections — `UpdateReweight.cxx:289`
  does `if (hdr.isCorrection) continue;` and the header is seeded to 1.0, so unit by design.
- `QEIntf_dial_5`: last QE-interference PCA component — ~zero eigen-variance, expected ≈1.
- `CCQETemplate ... q0bin4` (≈7e-5 %): the highest-q0 template bin has almost no events.

## Issue 2 — MEC q0q3 weights pile up at exactly 1000
`MECq0q3InterpWeighting_SBN_v3_SuSATo{Mar,Val}` set `WeightLimits: [0.1, 1000]` (fcl lines
1601, 1743). The (q0,q3) model-ratio interpolation extrapolates to extreme values in sparse /
edge cells — note the **coarse** `Q0Bins: [0, 0.2, 0.4, 0.6, 10]` (a single bin spans 0.6–10
GeV), `Q3ApplyMax: 2`, `UseNearestBin: false` — and those get clamped to exactly 1000. Only
~5–14 events out of 2.88 M hit the cap, but they dominate individual fine bins and produce the
"hot" cells in the q0/q3 plots. The *integrated* MEC effect is unaffected (~0.1–0.4 %).
→ **Fix:** tighten the cap, refine the q0 binning (the 0.6–10 bin is too wide), or
`UseNearestBin: true` to avoid edge extrapolation.

## Issue 3 — CCQE-template negative / huge weights (uncapped)
`CCQETemplateReweight_SBN_v3_{HFToCRPA,LFGToHF,LFGToSF}` reweight via the **ratio of two model
histograms** in (ENu, q3, q0) — e.g. `LowE/h_ENuq3q0_G21_HF-CRPA` ÷ `..._G21_HF` from
`$SBNDATA_DIR/systematics/AR23Plus/CCQETemplateReweight.root` (`RWMode: "q3q0"`). In sparse /
low-stat template cells (again coarse `q0_bin_edges` with a wide 0.7–10 bin) the denominator is
near-zero or noisy, so the ±1…±3σ-scaled ratio overshoots into large positive **and negative**
values (down to `wmin ≈ -318`). Crucially these three blocks have **no `WeightLimits`** (only
the 4th CCQE provider at fcl line 408 and the MEC providers do), so the pathological weights
pass straight through (and the analysis later clips <0 → 0, biasing those universes).
→ **Fix:** add `WeightLimits` to the HF/LFG/SF CCQETemplate blocks and/or floor/smooth the
template histograms; revisit the coarse q0 edges.

## Minor code note (`UpdateReweight.cxx`)
`sys_weights[pid]` is seeded to 1.0 once before the loop (line 212) and only overwritten for
params returned in each event's `resp` (lines 308–310) — never reset per event. It's correct
here because `GetEventVariationAndCVResponse` returns every configured param each event, but if
any provider ever returned a partial list, weights would go stale from the previous event.
Also `genieIdx = cafev_it` (line 265) assumes CAF↔genie row alignment (true for this file).
