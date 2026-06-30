# Reproducing the atmospherics xsec-systematics variation production

This documents the repos/branches and the end-to-end steps to regenerate the
per-event systematic weight files (`merged_new_systs_wsbn[_wheader].root`).

Everything runs in the **SL7 container**:
`/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest`
(or `fnal-dev-sl7` for interactive builds).

## 1. Repositories — all on branch `dune_atmospherics_fix_inf`

| Repo | Personal fork | Role / what was changed |
| --- | --- | --- |
| **GENIE Generator** | `github.com/pgranger23/Generator` | Custom GENIE with the **Martini MEC model** + MEC nucleon-cluster decay-angle generation (= `larsb-p/Generator`). Provides the alternative MEC model the reweight dials morph to. *Unmodified by us — forked for provenance.* |
| **GENIE Reweight** | `github.com/pgranger23/Reweight` | `GReWeightXSecMEC`: **DecayAng 0/0 `norm` NaN guard** at `twk_dial2=±0.5`; **`XSecShape` non-finite guard + `[0.1,1000]` clamp**; idempotent `Makefile` install. (Forked from `larsb-p/Reweight`.) |
| **nusystematics** | `github.com/pgranger23/nusystematics` | `GENIEReWeightEngineConfig.cc`: wire the **MEC morph dials** (`XSecShape_CCMEC_{Martini,Empirical}`, `EnergyDependence_CCMEC`) + the **`DecayAngMECVariationResponse`** response; **CCQETemplate q0-bin out-of-bounds segfault fix**. *Auto-fetched by CMake (see below).* |
| **atmospherics-tools** | `github.com/pgranger23/atmospherics-tools` | The production wrapper (this repo): fcl with the 16 unsupported **SBN-FSI dials removed**, `job.sh` that ships the **Martini hadron tensors** and builds a merged `genie_config/data` tree, `tools/systs.py` large-file header writer, build/submit/resubmit scripts, `ROOT_CAUSE_ANALYSIS.md`. |

Current tips (branch `dune_atmospherics_fix_inf`): Generator `5bcbee40b`,
Reweight `328501d3`, nusystematics `14c06d1`, atmospherics-tools `3539693`.

> nusystematics does **not** need a manual checkout — `xsec_systs_calculation/CMakeLists.txt`
> pulls it via `FetchContent` from `pgranger23/nusystematics @ dune_atmospherics_fix_inf`.
> `systematicstools` and the other deps are unmodified (default upstreams).

## 2. Build (in the SL7 container)

1. **Generator** — clone `pgranger23/Generator @ dune_atmospherics_fix_inf`, configure/build
   and `make install` into a GENIE install prefix (`local_install/`).
2. **Reweight** — clone `pgranger23/Reweight @ dune_atmospherics_fix_inf`, build against that
   Generator and `make install` into the same `local_install/` (provides the guarded `libGRwClc`).
3. **atmospherics-tools** — clone `pgranger23/atmospherics-tools @ dune_atmospherics_fix_inf`.
   In `xsec_systs_calculation/build_with_reweight.sh` **edit the hardcoded paths**
   (`GENIE_REWEIGHT=…/Reweight` and the `…/local_install/lib` copy lines) to your workspace.
   Then:
   ```bash
   cd xsec_systs_calculation
   ./make_tarball.sh        # runs build_with_reweight.sh in the SL7 container, then packages
   ```
   `make_tarball.sh` → `build_with_reweight.sh` CMake-fetches nusystematics, links the local
   Reweight, copies the `local_install` GENIE+Reweight libs into `build/Linux/lib`, builds the
   `UpdateReweight` app, and writes `../xsec_systs_calculation_wmec.tar.gz`.

The Martini MEC hadron tensors (absent from base ups GENIE) are shipped in
`genie_config/hadron_tensors_custom/martini/` and merged into `genie_config/data/` at
runtime by `job.sh`.

## 3. Submit (grid)

Needs a valid grid proxy/token in `jobscript-proxy.pem` (regenerate before submitting).
```bash
bash jobsub.sh           # STEP_SIZE=3000, -N 1000  -> /pnfs/.../xsec_systs_outputs_wsbn/
```

## 4. Recover failed / truncated jobs

`job.sh` has no `set -e`, so a crash still exits 0 with a truncated output. After the run, find
bad jobs and resubmit only those:
```bash
# missing or undersized (full ~2.94 MB / 3000 events) -> failed_jobs.txt
# edit the FAILED_JOBS array in job_resubmit.sh, then:
bash jobsub_resubmit.sh           # -N <n>, file:// job_resubmit.sh
# job_resubmit2.sh / jobsub_resubmit2.sh: same but --expected-lifetime=8h for slow ranges
```

## 5. Merge + add the self-describing header

```bash
hadd merged_new_systs_wsbn.root /pnfs/.../xsec_systs_output_*.root
cp merged_new_systs_wsbn.root merged_new_systs_wsbn_wheader.root
# in the SL7 env (has ROOT + uproot + fhicl_parser):
python -c "import sys; sys.path.insert(0,'tools'); \
  from systs import SystConfig, add_headers_to_root_file; \
  add_headers_to_root_file(SystConfig('fcl/systs_atmospherics_v2.fcl'), 'merged_new_systs_wsbn_wheader.root')"
```
`add_headers_to_root_file` appends the `systsHeader` tree (one row per dial: name, cv,
variations, id, isCorrection). It auto-uses a PyROOT `TFile`-UPDATE path for files >2 GB
(uproot's writer overflows a 32-bit key offset on large files).
