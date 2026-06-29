# Root-cause analysis — why some dials are inert / pathological

Investigation of the variation-production setup (`UpdateReweight.cxx`,
`fcl/systs_atmospherics_v2.fcl`, `genie_config/config/*`, job `.out` logs) explaining the
three issues seen in the wsbn output (`merged_new_systs_wsbn_wheader.root`). Generation tune:
**AR23_20i_00_000** (FSI = `genie::HAIntranuke2018/Default`, i.e. hA2018). The reweighting links
a **custom local GENIE build** (`/exp/dune/app/users/pgranger/systematics-lars/{Generator,Reweight,
local_install}`), NOT the ups `v3_06_02e` (that only provides base deps + data). This matters: the
fixes are about the custom build / its data, not the ups version.

## Issue 1 — 27 dials are identically unit. Three distinct causes:

### 1a. The 16 SBN FSI dials — no GReWeight calculator in this build, and N/A for hA2018
`FrG4{,LoE,M1E,M2E,HiE}_N`, `FrINCL{,LoE,M1E,M2E,HiE}_N`, `MFP{Lo,M1,M2,Hi}E_N`,
`FrKin_PiPro{Fix,Bias}_N` (the `GENIEReWeight ... SBN_v3` provider) are Geant4/INCL-cascade and
energy-binned MFP / kinematic FSI dials. The custom GReWeight build only implements the **standard
hA INuke dials** — `grep` of `Reweight/src/RwCalculators/GReWeightINuke*` gives exactly
`kINukeTwkDial_{FrAbs,FrCEx,FrElas,FrInel,FrPiProd,MFP}_{N,pi}` and **nothing** for
`FrG4/FrINCL/MFP*E/FrKin`. Consistently, `GSystUncertaintyTable.xml` has 1σ entries only for those
standard dials (the ones that **do** vary in the output). On top of that the sample's FSI is
**hA2018** (`AR23_20i/ModelConfiguration.xml: HadronTransp-Model = genie::HAIntranuke2018/Default`),
for which Geant4/INCL-cascade reweighting is not even defined. So these 16 dials cannot produce
weights here under any config.
→ **Fix (recommended):** remove the `GENIEReWeight ... SBN_v3` provider from the fcl. The hA FSI
uncertainties are already covered by the working `ICARUS_v2` `FrAbs/FrCEx/FrInel/FrPiProd/MFP`
dials. (Getting genuine G4/INCL-cascade FSI dials would require an SBN-fork GReWeight build **and**
regenerating the sample with a cascade FSI model — a separate, major effort, not warranted.)

### 1b. The MEC model-morph dials — Martini/Empirical hadron-tensor DATA is missing at runtime
`XSecShape_CCMEC_Martini`, `XSecShape_CCMEC_Empirical`, `EnergyDependence_CCMEC`,
`FracDelta_CCMEC`, `DecayAng2MEC`, `DecayAngMECVariationResponse` need a *live alternative MEC
cross-section model*. The job log shows it failing:
```
FATAL  ... Key: HadronTensorAlg does not exist in pools from algorithm :
       genie::MartiniEricsonChanfrayMarteauMECPXSec2024/Default
WARN   ... No Configuration available for ...MartiniEricsonChanfrayMarteauMECPXSec2024/Default
```
Root cause is a **data-path bug, not a config-registration one.** The configs are fine —
`master_config.xml` *does* register Martini/Empirical (lines 227/230/272/273), the model XMLs have
proper `Default` param_sets, and `reweight_master_config.xml` is additively loaded
(`AlgConfigPool.cxx:210-212`). But `MartiniMECHadronTensorModel.xml` loads tensors from
`DataPath: data/evgen/hadron_tensors/martini` *relative to `$GENIE`*, and `job.sh:46` does
`ln -s ${ORIG_GENIE}/data genie_config/data` with `ORIG_GENIE` = the **base ups `v3_06_02e`** —
whose `data/evgen/hadron_tensors/` contains only `crpa_susav2` and `nieves`, **no `martini`**. The
Martini/Empirical tensors exist **only in the custom build** (`Generator/data/evgen/hadron_tensors/
martini/`, `local_install/GENIE-Generator/data/...`). So `MartiniMECHadronTensorModel` can't load
its `.dat` files → can't configure → the Martini PXSec can't resolve `HadronTensorAlg` → the dial
→ 1.0. (The `MECq0q3Interp` SuSAToMar/Val dials work because they use **precomputed** weight files
from `sbndata`, not the live model.)
**Fix applied (data shipping):** ship the custom `martini` hadron tensors in the tarball
(`genie_config/hadron_tensors_custom/martini`, 16M) and build `genie_config/data` in `job.sh` as a
merged tree (base subdirs symlinked from `${ORIG_GENIE}/data` + the shipped martini dir).

**VALIDATION (ran the built `UpdateReweight` on 200 events with the fix):** the Martini `FATAL`
is **gone** and the model now instantiates — but the data fix is **necessary, not sufficient**.
`XSecShape_CCMEC_Martini` / `_Empirical` still output **1.0** on the exact CC-MEC events where the
generic `XSecShape_CCMEC` correctly morphs (e.g. evt 39: XSecShape=1.667, Martini=1.0). AR23_20i's
default MEC model is `genie::SuSAv2MECPXSec` (`AR23_20i/ModelConfiguration.xml`), which the morph
gate (`GReWeightXSecMEC::CalcWeightXSecShape_Martini`, line ~144) *does* handle — yet the weight
comes back unity. So there is a **second bug inside the custom `GReWeightXSecMEC.cxx`** weight
logic (the file is WIP, e.g. line 405 `TODO: Change this line once the Martini model is available
in GENIE`, line 380 `TODO: change hard-coding`). Pinning it needs a `pDEBUG`/`RwMEC` run to see the
per-event `diff_xsec_def`/`diff_xsec_alt`.

Note these dials are **redundant** with the working `MECq0q3InterpWeighting` SuSA→Martini/Valencia
dials (precomputed weight tables from `sbndata` — those already produce real MEC-model variation).

**Second bug — FOUND & FIXED.** Instrumented `GReWeightXSecMEC` with `std::cerr` and traced it: a
`ReW`=DEBUG run showed `SetSystematic` fires for `XSecShape_CCMEC`/`FracPN_CCMEC`/`FracDelta_CCMEC`
but **never for `_Martini`/`_Empirical`** → `CalcWeightXSecShape_Martini` is called but always with
`fCCXSecShapeMartiniTwkDial=0` → early return → unity. Root cause is **upstream in nusystematics**:
`GENIEReWeightEngineConfig.cc::ConfigureMECWeightEngine` only wires a *subset* of the MEC dials to
the GReWeight engine — it lists `XSecShape_CCMEC` but **omits** `_Empirical`, `_Martini`,
`EnergyDependence_CCMEC`, `DecayAng2MEC` (all of which `GReWeightXSecMEC::IsHandled` accepts). The
omitted dials are declared (so they appear in the header with ids) but never propagated to
`SetSystematic`, so they stay at 1.0.

**Fix:** add the four missing `kXSecTwkDial_*` to the `AddIndependentParameters({...})` list in
`ConfigureMECWeightEngine` (patch: `nusystematics_MEC_dials.patch`). **VALIDATED** by rebuilding
nusystematics and rerunning: `XSecShape_CCMEC_Martini` now varies (11/200 events, weights
0.28–1.43) and `_Empirical` too (0.016–3.74). `EnergyDependence_CCMEC`/`DecayAng2MEC` stay unit —
those have a *further*, separate calc-level issue (and `FracDelta_CCMEC`, already wired, likewise:
`CalcWeightPNDelta` returns unity for it) — out of scope of this propagation fix.

→ **To make it stick in production:** the fix must land in the nusystematics fork that this project
FetchContent-pulls — `github.com/pgranger23/nusystematics`, branch `dune_atmospherics_fix_inf`
(file `src/nusystematics/systproviders/GENIEReWeightEngineConfig.cc`). Apply
`nusystematics_MEC_dials.patch` there and push; a clean `build_with_reweight.sh` will then pick it
up. (A local validated build with the fix is in `build/Linux/lib/libnusystematics_systproviders.so`.)
**PUSHED** to the fork's `dune_atmospherics_fix_inf` branch.

### 1d. The other inert MEC dials — investigated; each a distinct, deeper issue
Wiring the engine makes `XSecShape_CCMEC_{Martini,Empirical}` work (validated: weights ~0.28–1.43 /
0.015–3.74). The remaining inert MEC dials are NOT the same bug — instrumented (`std::cerr`) runs
pinned each:
- **`FracDelta_CCMEC`** — by design unavailable for SuSAv2. In `GReWeightXSecMEC::CalcWeightPNDelta`
  the delta-fraction term is **commented out** (~lines 860-861, "DeltaNotDelta only works for
  Valencia"); the SuSAv2 branch (~line 777) computes only `pn_frac`, never `delta_frac`. AR23_20i is
  SuSAv2 ⇒ no effect. Needs a SuSAv2 delta-fraction implementation (the file's own TODO).
- **`EnergyDependence_CCMEC`** — wiring is correct (its `SetSystematic` now fires), but
  `BuildEnergyDepRatioGraphs` yields a **trivial envelope**: the trace shows `r_upper = r_lower = 1`
  at every energy, so `weight = 1 + twk·(r−1) = 1`. The alt models used for the envelope don't differ
  from SuSAv2 in normalised energy shape — a model/config issue in the envelope construction, not the
  dial wiring. (Kept in the fix since the wiring is correct.)
- **`DecayAng2MEC`** — a **responseless** "frequency" param. In `CalcWeightAngularDist` the weight is
  `3·twk_dial·cos²(twk_dial2·θ) + (1−twk_dial)` with `twk_dial`=DecayAngMEC (amplitude),
  `twk_dial2`=DecayAng2MEC (frequency); when DecayAngMEC=0 (its nominal, as when DecayAng2MEC is
  varied alone) it collapses to 1 regardless of DecayAng2MEC. Its real effect is delivered jointly via
  the `DecayAngMECVariationResponse` response param (`AddResponseAndDependentDials`, not
  `AddIndependentParameters`). So it was **removed** from the independent-dial fix.

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

## Issue 4 — DecayAngMECVariationResponse stayed 1.0 (now fixed)
The MEC nucleon-cluster decay-angle systematic is delivered through the **response** parameter
`DecayAngMECVariationResponse`, whose 25 universes index the joint grid of its dependent
*responseless* dials `DecayAngMEC` (amplitude) and `DecayAng2MEC` (frequency). Two bugs kept it at
1.0:
1. **No wiring (nusystematics).** `ConfigureMECWeightEngine` only called `AddIndependentParameters`
   and never `AddResponseAndDependentDials` for any MEC response — every other channel (QE/RES/NCEL/
   COH/DIS/FSI) wires its `*VariationResponse`, but MEC did not. So the response was declared (header
   id) but never connected to a GReWeight engine. Also `DecayAngMEC` was (incorrectly) in the
   independent list, giving it only the standalone amplitude term.
   → **Fix:** add `AddResponseAndDependentDials(MECmd, "DecayAngMECVariationResponse",
   {kXSecTwkDial_DecayAngMEC, kXSecTwkDial_DecayAng2MEC}, "xsec_mec", ...)` and drop `DecayAngMEC`
   from the independent list (`nusystematics_DecayAng_response.patch`).
   Pushed to `pgranger23/nusystematics @ dune_atmospherics_fix_inf`.
2. **0/0 singularity (Reweight).** Wiring the response exposed a latent NaN in
   `GReWeightXSecMEC::CalcWeightAngularDist`: the normalization has a `(1 - 4*twk_dial2^2)`
   denominator that vanishes at `twk_dial2 = ±0.5` (a value in the frequency grid) → NaN. The limit
   of the singular term is 0.
   → **Fix:** guard the denominator and substitute the limit
   (`Reweight_DecayAng_nan_guard.patch`). Pushed to `pgranger23/Reweight @ dune_atmospherics_fix_inf`
   (forked from larsb-p/Reweight; the local Reweight origin was repointed to the fork, and the
   guarded `libGRwClc` installed into `local_install/lib`).

**Validated:** `DecayAngMECVariationResponse` now varies (no NaN/inf); `DecayAngMEC`/`DecayAng2MEC`
standalone branches correctly go to 1.0 (responseless). Note: at negative amplitude (`DecayAngMEC` =
−0.5/−1) the angular model dips negative (~−2) — an inherent property of the linear isotropic↔cos²
interpolation extrapolated to negative amplitude; these are clipped to 0 by the analysis tooling
(`systs.py` `clip(lower_bound=0)`).

## Issue 5 — XSecShape CCMEC morph weights: intermittent NaN/inf + extreme outliers (guarded)
After wiring the DecayAng response, a clean production build showed `XSecShape_CCMEC_Martini`/
`_Empirical` (and the generic `XSecShape_CCMEC`) emit **NaN/inf** at a couple of pathological MEC
events (53, 163), plus large finite outliers (generic up to ~162×). Investigation: the same events
on a *rebuild* came back finite — i.e. it's an **intermittent, build-sensitive numerical
instability** (the likelihood ratio `tweaked_prob_density / prob_density_def` degenerating when the
default-model density ~0), **not a deterministic regression** from the DecayAng work. NaN/inf is not
caught by the downstream `clip(lower_bound=0)` and would corrupt the histograms.
→ **Fix (`Reweight_XSecShape_guard_clamp.patch`, pushed to `pgranger23/Reweight`):** in all three
`CalcWeightXSecShape*` functions, substitute unity for non-finite weights and clamp to MEC-style
limits `[0.1, 1000]` (matching `MECq0q3InterpWeighting`'s `WeightLimits`). Validated: no NaN/inf;
DecayAng response unaffected. The guarded `libGRwClc` was installed into `local_install/lib` so a
clean `build_with_reweight.sh` ships it.

## Issue 6 — CCQETemplateReweight segfault on high-q0 events (grid job 69)
Grid job 69 (events 207000-210000) **segfaulted** in `CCQETemplateReweight::GetEventResponse`
(`CCQETemplateReweight_tool.cc:234`). Root cause: `GetQ0BinIndex` returns **`Nq0Bins`** for any
event with `q0 >= q0BinEdges[Nq0Bins]` (the top bin edge), but `ResponseParameterIndices` and the
`resp` vector are sized `Nq0Bins` (valid indices `0..Nq0Bins-1`). `GetEventResponse` then indexes
`ResponseParameterIndices[Nq0Bins]` / `resp[Nq0Bins]` **out of bounds → SIGSEGV**. Triggered by rare
high-q0 (high-energy atmospheric) events (~1 in 2e5), so only ~1 job in ~960 hits it per run.
**Severity is amplified by `job.sh` having no `set -e`:** the crash truncates the output yet the job
still `ifdh cp`s the partial file and exits 0 — so truncated outputs look "successful".
→ **Fix (`nusystematics_CCQETemplate_q0bin_oob.patch`, pushed to `pgranger23/nusystematics`,
`14c06d1`):** return the last valid bin (`Nq0Bins-1`) for the q0-overflow case, plus a defensive
clamp of `q0_bin_index` in `GetEventResponse`. **Validated:** re-running the job-69 range
(`-s 207000 -N 1000`, incl. entry 207567) completes with no segfault, all 1000 events written,
all CCQETemplate dials finite. The fixed `libnusystematics` is in the repackaged
`xsec_systs_calculation_wmec.tar.gz`. NOTE: any job whose output is short/undersized
(full ~4.8 MB / 3000 events) likely crashed and must be resubmitted.
