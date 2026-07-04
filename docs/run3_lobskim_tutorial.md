# Run 3 LobSkim tutorial

## Scope and assumptions

This tutorial documents the current Run 3 `LobsterSkimming` parent repository at `/users/apiccine/work/LobSkim/LobsterSkimming`. In this workspace, the local working branch is `run3-mva-anpicci`; the upstream branch this work is intended to land in is `run3_mva`.

The source descriptions below are based on inspection of the current files. The validation-state section separately records bounded evidence from explicit-input smoke preparation and two one-sample worker-backed `TESTING` runs. That evidence does not establish full production readiness or campaign-scale behavior.

Run 3 and Run 2 live in separate parent repositories and can have different CMSSW releases, wrapper behavior, input-mode semantics, cfg naming, and campaign paths. Do not copy Run 2 values into this repository unless the current Run 3 source proves the same value.

## Repository layout

| path | ownership and purpose |
| --- | --- |
| `README.md` | Concise Run 3 landing page. |
| `docs/run3_lobskim_tutorial.md` | This branch-specific operational reference. |
| `setup.py` | No-argument setup entry point; dispatches sparse `topeft` setup and CMSSW setup when their target paths are absent. |
| `scripts/install_configs.sh` | Sparse-checkout helper for selected cfg and sample-JSON paths from `TopEFT/topeft`. |
| `scripts/install_cmssw.sh` | Builds the Run 3 CMSSW/payload area when intentionally executed by a user. |
| `skimmer/lobster_config.py` | User knobs, validation, derived labels/paths, storage, workflows, and Lobster `Config`. |
| `skimmer/skim_wrapper.py` | Worker payload for stage-in, NanoAOD postprocessing, and per-task output merge. |
| `skimmer/tools/utils.py` | Reads cfg entries and referenced sample JSON files. |
| `skimmer/factory_skim.json` | Source-controlled Run 3 factory resource profile; not a complete factory command. |
| `topeft/` | Independent nested Git repository supplying cfg and JSON metadata; protected from parent-repository cleanup. |
| `CMSSW_14_0_6/` | Expected Run 3 CMSSW sandbox path when present locally; directory presence alone does not prove build completeness. |

The workspace also contains separate repositories for Run 2, Lobster, WMCore, and nested Run 2 `topeft`. Confirm the repository path before editing, staging, or running anything.

## Branch and checkout policy

The Run 3 upstream branch intended to receive and host this work is `run3_mva`. The local workspace branch used for integration and local validation is `run3-mva-anpicci`; do not treat that local branch name as the long-term public branch.

Fresh setup is configured to obtain `https://github.com/TopEFT/topeft.git` at branch/tag `run3_test_mmerged`. That fresh-setup policy and an existing local checkout are different facts. The recorded workspace checkout is on `run3_test_mmerged_anpicci` at `06f84d838f3aed6cc18fdd1a301b1fb5fb651624` and has a local cfg modification. Do not switch it as part of normal tutorial use. Reconciling `run3_test_mmerged` and `run3_test_mmerged_anpicci` is a separate maintenance decision requiring explicit authorization and a clean plan for local changes.

Do not fetch, switch, merge, rebase, reset, clean, or stash merely to make local names match this tutorial. First identify the repository boundary and preserve local work.

## Protected nested topeft repository

`topeft/` is a complete, independent Git repository. The Run 3 parent may display it as untracked because it is not a submodule. That output does not make it disposable.

- Run Git inspection for `topeft` from `LobsterSkimming/topeft`, not from the parent.
- Do not clean or delete it from the parent.
- Do not stage it into the parent.
- Do not assume authorization to edit the parent also permits editing `topeft`, or vice versa.
- Preserve an existing checkout even when `setup.py` would choose a different branch for a fresh install.

`setup.py` only checks whether the `topeft` path exists. If it exists, setup skips installation and does not verify its branch or content.

## Environment policy

Keep setup/build and Lobster-control concerns explicit:

1. CMSSW release creation and compilation belong in a clean shell outside conda/mamba.
2. Lobster control belongs in the site-approved Lobster environment and must use the activation sequence approved for this workspace.
3. The worker sandbox is fixed by `sandbox_location` in `skimmer/lobster_config.py`; it is not inferred from whichever shell happens to be active.
4. Do not mutate a shared environment to repair a campaign. Record the missing dependency and handle environment maintenance separately.

The manager environment does not prove which interpreter the worker payload resolves inside its sandbox. Record both environments separately when validating a campaign.

## Setup overview

`setup.py` has a no-argument entry point and immediately calls `main()` when executed by a user. It finds the parent Git top level, then:

- skips `topeft` when that path already exists;
- otherwise calls `scripts/install_configs.sh` for the two sparse paths;
- skips the CMSSW installation when `CMSSW_14_0_6` already exists;
- otherwise calls `scripts/install_cmssw.sh` with the parent path, release, and SCRAM architecture.

Both shell installers use `set -euo pipefail` and positional-argument checks. This describes source-visible behavior only; no GitHub, CVMFS, SCRAM, or storage execution is recorded for the current setup.

## Conda/mamba environment

Do not run CMSSW creation or compilation from an active conda/mamba environment. Use a separate clean shell for setup. When controlling Lobster, activate only the established Lobster environment; do not install or upgrade packages ad hoc during a campaign.

Before production, record the exact activation sequence used at the site, the resulting Python executable, relevant environment variables, proxy/credential state, and the CMSSW shell state. This repository does not currently provide a canonical user-facing activation command, so this tutorial does not invent one.

## CMSSW setup

The Run 3 setup values are:

| item | current source value |
| --- | --- |
| release | `CMSSW_14_0_6` |
| SCRAM architecture | `el9_amd64_gcc12` |
| CVMFS setup root | `/cvmfs/cms.cern.ch` |
| CMGTools repository | `https://github.com/anpicci/topEFT_ttHMVA_Run3.git` |
| CMGTools branch | `nd_run3` |
| commented historical payload references | `sscruz/cmgtools-lite`, `jdelrieg/topEFT_ttHMVA_Run3`, and `cms-nanoAOD/nanoAOD-tools` appear only as comments in the current installer. |

The installer refuses to continue when `CMSSW_BASE` is already set, sources CMS defaults from CVMFS, exports `SCRAM_ARCH`, creates the CMSSW release, clones the Run 3 CMGTools payload under `src/CMGTools`, evaluates the release runtime, and runs `scram b`. These operations are network- and build-heavy and require separate approval/validation. Do not treat the presence of a directory as proof that the release or payload is complete.

`skimmer/lobster_config.py` independently fixes the worker sandbox to `<Run3-parent>/CMSSW_14_0_6`; if the release location changes, setup and config must be reviewed together.

## topeft sample cfg/json setup

Fresh setup uses:

| item | current value |
| --- | --- |
| repository | `https://github.com/TopEFT/topeft.git` |
| branch/tag | `run3_test_mmerged` |
| sparse path 1 | `input_samples/cfgs` |
| sparse path 2 | `input_samples/sample_jsons` |

These paths resolve under the independent checkout as `topeft/input_samples/cfgs` and `topeft/input_samples/sample_jsons`. A campaign cfg is resolved as `topeft/input_samples/cfgs/<CFG_NAME>`.

`scripts/install_configs.sh` initializes the target as a sparse Git checkout when needed, configures the requested sparse paths, verifies or adds `origin`, fetches the requested selector with `--depth=1`, and checks out either a local branch, remote branch, or tag named by the selector. It is not a general-purpose cleanup or reconciliation tool for an existing local checkout.

`read_cfg` ignores blank/commented cfg lines, treats `root:` lines as cfg metadata, applies every regex in `MATCH`, resolves JSON paths relative to the cfg directory, and loads each selected JSON. The sample key is the JSON basename without `.json`.

## Source files used by a skim

A campaign is assembled from:

1. `skimmer/lobster_config.py` for campaign knobs and workflow objects.
2. `topeft/input_samples/cfgs/<CFG_NAME>` for the active JSON list.
3. Each referenced JSON for `files`, `path`, `year`, `isData`, and analysis metadata.
4. `skimmer/skim_wrapper.py` as the worker command and extra input.
5. `CMSSW_14_0_6` as the Lobster sandbox.
6. Lobster framework classes imported from the local Lobster environment: `AdvancedOptions`, `Category`, `Config`, `Dataset`, `StorageConfiguration`, `Workflow`, and `cmssw.Dataset`.

## Quick start checklist

Before requesting a production launch:

1. Confirm the parent is the Run 3 clone and inspect its branch/status without altering it.
2. Inspect nested `topeft` separately; record its branch, commit, and dirty state.
3. Confirm the expected CMSSW release and payload paths exist; do not infer completeness from names alone.
4. Edit only the intended user knobs in `skimmer/lobster_config.py`.
5. Verify `TYPE`, `YEAR`, and the derived `CFG_NAME` identify an existing cfg.
6. Review every selected JSON for the metadata required by the chosen `INPUT_MODE`.
7. Calculate the derived label, workdir, plotdir, and output path before launch; ensure they do not collide with an existing campaign.
8. Obtain user-approved Work Queue and Lobster commands; the placeholders below are deliberately non-runnable.
9. Record credentials, storage reachability, quotas, factory resources, logs, and shutdown procedure.
10. Run production only in a separately authorized round and retain the exact command/output record.

## Editing skimmer/lobster_config.py

The `USER KNOBS` block is the intended campaign-editing surface. Avoid changing helpers, workflow construction, storage logic, or wrapper code just to launch a different sample set.

For common goals:

- To switch from background to signal samples, change `TYPE` to `"signal"` and confirm the derived cfg exists for the selected `YEAR`.
- To run data, change `TYPE` to `"data"`; the config derives and validates `<YEAR>_data.cfg`.
- To select a different Run 3 period, change `YEAR` and review the resulting `TAG`, `CFG_NAME`, and output paths.
- To narrow a retry or debugging subset, edit `MATCH` with regexes that match cfg lines.
- To choose DBS discovery instead of explicit file lists, change `INPUT_MODE` to `"dbs"` and verify every selected JSON has a valid `path`.
- To choose per-sample input discovery, change `INPUT_MODE` to `"auto"` and review the resolved DBS/files counts and mixed-storage profile.
- To run SR-like skim cuts, change `TARGET` to `"SR"`; leave `TARGET="CR"` for the looser control-region cut.

After edits, review the complete resolved configuration output produced when Lobster eventually evaluates the config in an authorized run; do not rely only on filenames.

## User knobs reference

| knob | default | allowed_or_expected_values | what_it_controls | when_to_edit | risk_or_caveat | source_file |
| --- | --- | --- | --- | --- | --- | --- |
| `TESTING` | `False` | Boolean | Switches labels and paths to timestamped `test/lobster_skimtest_*` locations. | For a deliberately bounded test campaign. | A test-labelled path does not make submission safe or lightweight. | `skimmer/lobster_config.py` |
| `INPUT_MODE` | `"files"` | `dbs`, `files`, or `auto` after lower-case normalization. | Per-sample dataset object type and campaign storage profile. | Use `files` for explicit lists, `dbs` for dataset discovery, or `auto` for metadata-based per-sample selection. | `auto` gives a valid DBS dataset path precedence over files and can produce a mixed campaign. | same |
| `PROTOCOL_LOCAL` | `"file://"` | Storage protocol understood by Lobster; current design expects `file://`. | Local input/output prefix. | Only for an intentional storage-profile change. | `file://` requires worker-visible paths and existing destination parents. | same |
| `PROTOCOL_REMOTE` | `"root://"` | Current design expects `root://`. | Remote XRootD prefix. | Only when changing remote transport. | Endpoint reachability was not validated. | same |
| `TARGET` | `"CR"` | `SR`, `CR` after upper-case normalization | Changes the skim cut and path/label segment. | Choose signal-region or control-region skim semantics. | `SR` adds the opposite-sign two-lepton veto; `CR` does not. | same |
| `YEAR` | `"2023BPix"` | Validated as `2022`, `2022EE`, `2023`, or `2023BPix`, with case-insensitive normalization to the canonical spelling. | `TAG`, `CFG_NAME`, paths, and module validation. | Choose the Run 3 run period. | Adding a year requires all canonical cfgs, updated source validation, and metadata review. | same |
| `STEP` | `"skimmed"` | Nonempty campaign label chosen by operator. | Label/path segment. | When naming a distinct processing step. | Changing it creates different paths; it does not change physics behavior. | same |
| `TYPE` | `"background"` | Validated as `background`, `signal`, or `data` after lower-case normalization. | `TAG` and `CFG_NAME`. | Match the campaign kind. | Unsupported values fail before paths and workflows are built. | same |
| `TAG` | `f"{TYPE}/NAOD_ULv9_lepMVA-run3/{YEAR}"` | Derived string unless deliberately edited. | Campaign/output namespace. | Rarely; usually edit `TYPE` or `YEAR`. | Manual edits can make output identity misleading. | same |
| `CFG_NAME` | `ND_2023BPix_background_samples.cfg` for defaults | Derived as `ND_<YEAR>_<TYPE>_samples.cfg` for background/signal and `<YEAR>_data.cfg` for data, then checked against the cfg directory. | Selects the sample list. | Not an independent user knob in the current source; change `TYPE` and `YEAR`. | A missing canonical cfg fails with the expected name and a bounded available-file preview. | same |
| `MATCH` | `[]` | List of regular-expression strings | Filters cfg JSON lines; empty selects all active lines. | For a narrow subset or retry. | Regex matches the cfg line; verify selected count and names. | same |
| `WRAPPER` | `"skim_wrapper.py"` | Worker wrapper present in the skimmer working context. | Payload command and `extra_inputs`. | Only when intentionally replacing the worker contract. | Wrapper changes require separate validation. | same |
| `SRC_REMOTE` | `"cms-xrd-global.cern.ch"` | XRootD host, without protocol. | Remote stage-in prefix and DBS `xrootd_servers`. | When site/storage policy changes. | Network access and fallback order must be validated. | same |
| `SRC_LOCAL` | `"/cms/cephfs/data"` | Worker-visible absolute path. | Local file-input prefix. | When local storage mount changes. | Must exist on workers; no check is performed statically. | same |
| `DST_REMOTE` | `"cmsxrootd.crc.nd.edu"` | XRootD destination host. | Remote stage-out prefix. | When output endpoint changes. | Credentials, quota, and permissions must be checked. | same |
| `DST_LOCAL` | `"/cms/cephfs/data"` | Worker-visible absolute path. | Local stage-out prefix. | When local destination mount changes. | Destination parents must already be reachable for `file://`. | same |
| `WORKDIR_BASE` | `"/tmpscratch/users/$USER"` | Writable manager-local base path. | Lobster state/database directory. | When project state belongs on another filesystem. | Never reuse a workdir for an unrelated campaign. Shell variable expansion behavior must be confirmed in the actual Lobster environment. | same |

## Derived labels and paths

`TSTAMP1` is `YYYYMMDD_HHMM`; `startingday` is `YYMMDD`; `ver` is `vYYMMDD`.

| label_or_path | source_expression_summary | default_or_pattern | meaning | when_it_changes | risk_or_caveat |
| --- | --- | --- | --- | --- | --- |
| `TAG` | `<TYPE>/NAOD_ULv9_lepMVA-run3/<YEAR>` | `background/NAOD_ULv9_lepMVA-run3/2023BPix` | Campaign/output namespace. | Validated `TYPE` and `YEAR`. | It is derived after type/year validation but before sample metadata is loaded. |
| `master_label` | `<STEP>_<TARGET>_<INPUT_MODE>_lobPY3_<timestamp>` | `skimmed_CR_files_lobPY3_YYYYMMDD_HHMM` | Lobster config label. | Step, target, mode, evaluation minute; test has `testlobPY3`. | Label says `lobPY3`; actual payload command starts with `python3`. |
| `sandbox_location` | `<git-top>/CMSSW_14_0_6` | Absolute Run 3 release path | Worker CMSSW sandbox. | Parent checkout location or source edit. | Must match a valid built release. |
| `cfg_fpath` | `<git-top>/topeft/input_samples/cfgs/<CFG_NAME>` | Default ends in `ND_2023BPix_background_samples.cfg` | Selected cfg. | `CFG_NAME`, `TYPE`, `YEAR`, or checkout path. | Nested checkout ownership applies. |
| `workdir_path` | `<WORKDIR_BASE>/<STEP>_<TARGET>/<TAG>/vYYMMDD` | `/tmpscratch/users/$USER/skimmed_CR/background/NAOD_ULv9_lepMVA-run3/2023BPix/vYYMMDD` | Manager state. | Base, step, target, type, year, day; test uses timestamped test path. | Same-day repeated non-test launches collide unless intentionally resumed. |
| `plotdir_path` | `~/www/lobster/<STEP>_<TARGET>/<TAG>/vYYMMDD` | Tilde-based web path | Monitoring plots. | Step, target, type, year, day; test uses timestamp. | Confirm tilde handling and directory permissions. |
| `output_path` | `/store/user/$USER/<STEP>_<TARGET>/<TAG>/vYYMMDD` | Default Run 3 stage-out namespace | Logical output path appended to destination prefixes. | Step, target, type, year, day; test uses timestamp. | Confirm `$USER`, quota, permissions, and non-collision. |
| local source prefix | `PROTOCOL_LOCAL + SRC_LOCAL + "//"` | `file:///cms/cephfs/data//` | Files-mode local input prefix. | Local protocol/root. | Workers must see the path. |
| remote source prefix | `PROTOCOL_REMOTE + SRC_REMOTE + "//"` | `root://cms-xrd-global.cern.ch//` | Files-mode remote input prefix and DBS AAA host. | Remote protocol/host. | Not network-validated here. |
| local output URL | `file://` + local root + `//` + output path | `file:///cms/cephfs/data///store/user/$USER/...` | First files-mode stage-out method. | Protocol/root/output. | Repeated slashes are source-visible; worker filesystem must be accessible. |
| remote output URL | `root://` + remote host + `//` + output path | `root://cmsxrootd.crc.nd.edu///store/user/$USER/...` | DBS output and files-mode fallback. | Protocol/host/output. | Not network-validated here. |

## TYPE, YEAR, and CFG_NAME behavior

`TYPE` and `YEAR` drive both output identity and cfg selection:

```text
TAG = f"{TYPE}/NAOD_ULv9_lepMVA-run3/{YEAR}"
CFG_NAME = f"{YEAR}_data.cfg" if TYPE == "data" else f"ND_{YEAR}_{TYPE}_samples.cfg"
```

With the current defaults, `TYPE="background"` and `YEAR="2023BPix"` derive:

```text
TAG = background/NAOD_ULv9_lepMVA-run3/2023BPix
CFG_NAME = ND_2023BPix_background_samples.cfg
```

The source normalizes and validates `TYPE` and `YEAR` before deriving `TAG`. It then lists the cfg directory, derives exactly one canonical name, and requires that name to exist. Background and signal use `ND_<YEAR>_<TYPE>_samples.cfg`; data uses `<YEAR>_data.cfg`. There is no broad filename search and no independent `CFG_NAME` override in the current source. A missing cfg error reports the requested type/year, expected name, and a bounded preview of available cfgs.

The resolver intentionally chooses only canonical cfg names. It does not automatically discover alternate `NDSkim_*`, CR/SR-specific, central-signal, signal-subtype, synchronization, or general cfgs. Users should not expect `TYPE` and `YEAR` to discover those alternate cfg families automatically. Supporting those files requires an explicit selection design rather than a broad filename search.

`select_module_name` accepts current Run 3 sample names containing `2022` or `2023` and returns the only factory exposed by the wrapper's import path: `lepMVA`. It no longer returns unsupported Run 2 fallback names. An unmappable sample fails before workflow construction with its sample name, selected year, cfg, and supported module list. This config does not filter samples by JSON `year`; the selected cfg and `MATCH` determine the active sample list.

## TARGET and skim-cut behavior

Only `SR` and `CR` are valid. Both require at least two selected electrons, muons, and taus in total. Electron, muon, and tau object requirements are encoded directly in `build_skim_cut`.

`SR` additionally vetoes events with exactly two selected leptons whose summed charge is zero. `CR` uses only the at-least-two requirement. `TYPE` does not affect this cut. The generated expression has whitespace removed and is checked for balanced parentheses and accidental literal outer quotes.

## INPUT_MODE and sample input metadata

| mode | required metadata | dataset object | storage behavior |
| --- | --- | --- | --- |
| `files` | Nonempty `files` list containing usable file strings. | Lobster `Dataset(files=..., files_per_task=1, patterns=["*.root"])` | `storage_files`, with local and remote input prefixes and local and remote output prefixes. |
| `dbs` | `path` string suitable for `cmssw.Dataset(dataset=...)`. | `cmssw.Dataset(dataset=jsn["path"], lumis_per_task=1, file_based=True)` | `storage_dbs`, remote output only, and `AdvancedOptions.xrootd_servers=[SRC_REMOTE]`. |
| `auto` | A DBS-form `path`, or a nonempty usable `files` list. | Chosen per sample using the same dataset constructors as explicit modes. | DBS path takes precedence; files are the fallback. Mixed campaigns use `storage_files` and enable `xrootd_servers`. |

The current source accepts all three modes. A DBS path must have `/PrimaryDataset/ProcessedDataset/DataTier` form with a recognized CMS data tier; `/store/...`, `root://...`, `file://...`, and other absolute storage paths are not DBS dataset names. Explicit `files` and `dbs` validate their required metadata. In `auto`, a valid DBS path wins when both inputs are usable, matching the Run 2 policy; otherwise the config uses files. Missing or unusable metadata fails with the sample, cfg path, requested mode, available keys, and a bounded metadata summary.

The mixed `auto` policy is source-defined and covered by lightweight helper/static validation. No real same-context cfg containing both DBS-backed and files-backed samples was found, so mixed `auto` has not been exercised in a Lobster runtime. Before using mixed files/DBS `auto` campaigns for production, run a bounded validation that checks Lobster, Work Queue, storage, XRootD, and DBS behavior.

Each selected JSON should also carry analysis metadata such as `year`, `isData`, `xsec`, event counts, and weight sums. The current Run 3 config consumes `files` and `path` for input construction; it does not enforce or filter on JSON `year`.

The resolved summary reports the implemented auto status, per-sample input-mode counts, selected Config-level storage profile/counts, cfg-name source and naming rule, type/year validation status, XRootD server behavior, and supported/selected/rejected lepMVA module names. Sample and command previews remain bounded.

## Storage and XRootD behavior

`storage_dbs` has only the remote output URL and leaves input streaming enabled. `storage_files` tries local then remote input prefixes and local then remote output prefixes, also with streaming enabled.

Selection is:

- `INPUT_MODE="files"` uses `storage_files`.
- `INPUT_MODE="dbs"` uses `storage_dbs`.
- `INPUT_MODE="auto"` uses `storage_dbs` when every selected workflow is DBS-based; if any workflow uses explicit files, including a mixed campaign, it uses `storage_files`.

`AdvancedOptions.xrootd_servers` is set to `["cms-xrd-global.cern.ch"]` whenever at least one selected workflow is DBS-based, including mixed auto campaigns. File-based workflows receive explicit input file names and may stage in missing local basenames with `xrdcp -f`.

The framework `StorageFileURLs.md` notes that `file://` output is a local filesystem stage-out method: the worker must see the destination parent directory, and Lobster does not auto-create missing parents during task stage-out. The remote `root://` path is a fallback/output target but was not tested here.

## Workflow and AdvancedOptions behavior

Each selected sample becomes one workflow with hyphens replaced by underscores in its label.

| field | current value or rule | operational meaning |
| --- | --- | --- |
| category | `processing` | Shared resource category. |
| category resources | 1 core, 1500 MB memory, 4500 MB disk | Initial per-task request. |
| sandbox | `CMSSW_14_0_6` | Worker software release. |
| task splitting | one lumi per task for DBS; one file per task for files | Dataset-dependent unit construction. |
| `extra_inputs` | `[WRAPPER]` | Ships `skim_wrapper.py` with tasks. |
| `outputs` | `output.root` | Required wrapper product. |
| command | `python3 <wrapper> --cut <cut> --module <module> --out-dir . @inputfiles` | Actual Lobster payload string. |
| `merge_command` | `haddnano.py @outputfiles @inputfiles` | Lobster-level output merge. |
| `merge_size` | `537M` | Merge target size. |
| `globaltag` | `False` | No CMSSW global-tag setup requested by workflow. |
| `cleanup_input` | `False` | Input cleanup disabled. |

`AdvancedOptions` currently sets `bad_exit_codes=[127, 160]`, `log_level=1`, `payload=10`, `osg_version="3.6"`, `threshold_for_failure=10`, and `threshold_for_skipping=10`, plus `xrootd_servers` only for DBS mode. Changing thresholds affects retry behavior and should be based on diagnosed failures, not used to hide a persistent root cause.

## Wrapper command construction

| argument_or_behavior | source_file | Run_3_behavior | user_should_edit | risk_or_caveat |
| --- | --- | --- | --- | --- |
| interpreter | `skimmer/lobster_config.py` | Payload starts with `python3`. | No, unless separately validating the sandbox/runtime contract. | The worker environment must resolve `python3` correctly inside the sandbox. |
| wrapper | `skimmer/lobster_config.py`, `skimmer/skim_wrapper.py` | Default `skim_wrapper.py`; shipped as `extra_inputs`. | Normally no. | Replacement changes the worker contract. |
| `--cut` | config helper and wrapper parser | Receives whitespace-free generated skim expression. | Change `TARGET`, not command text. | Shell quoting is intentionally not inserted into the actual Lobster command. |
| `--module` | config helper and wrapper parser | `lepMVA` for sample names containing `2022` or `2023`; other names fail before workflow construction. | Normally no; fix sample naming only through authorized metadata changes. | The wrapper import path exposes only the validated `lepMVA` factory. |
| `--out-dir` | config helper and wrapper parser | `.` | Normally no. | Wrapper and expected output assume task working directory. |
| `--nevents` | wrapper parser only | Optional wrapper flag exists, but current config does not pass it. | Not through current user knobs. | Adding it requires workflow-command review. |
| `@inputfiles` | Lobster substitution | Expanded to task inputs. | No. | Expansion semantics belong to Lobster. |
| `output.root` | wrapper/workflow | Final required task output. | No. | Missing output causes task failure/stage-out failure. |
| workflow merge | config | `haddnano.py @outputfiles @inputfiles`. | Only in a separately reviewed workflow change. | Distinct from wrapper's per-task merge. |

The config also creates a `shlex.join` display command for logs. That display form is not passed to Lobster because parser handling of shell quoting has not been established.

## lepMVA module wiring and output branches

### Where lepMVA modules are defined

The Run 3 installer places the CMGTools payload in the `CMSSW_14_0_6` sandbox by cloning the `nd_run3` branch of `anpicci/topEFT_ttHMVA_Run3` in `scripts/install_cmssw.sh`. In the current local payload, `CMSSW_14_0_6/src/CMGTools/NanoProc/python/tools/nanoAOD/lepMVA_run3.py` defines the `lepMVA_run3` class and the `lepMVA` factory. The wrapper imports that file through the Python path `CMGTools.NanoProc.tools.nanoAOD.lepMVA_run3` when `nano_postproc.py` starts.

The imported file defines `lepMVA`; it does not define `lepMVA_2016_preVFP`, `lepMVA_2016`, `lepMVA_2017`, or `lepMVA_2018`. The config now treats those names as unsupported fallbacks and fails instead of returning them through this import path.

### How lepMVA is selected and wired into nano_postproc.py

The wiring is per sample and follows this path:

1. `select_module_name` in `skimmer/lobster_config.py` returns `lepMVA` when the sample name contains `2022` or `2023`; any other sample name fails with the selected year, cfg, and supported module list.
2. `build_payload_command` adds `--module <module_name>` to the exact string assigned to `Workflow(command=command_string)`.
3. `skimmer/skim_wrapper.py` requires and parses `--module`, then constructs `nano_postproc.py -I CMGTools.NanoProc.tools.nanoAOD.lepMVA_run3 <module>`.
4. For current Run 3 sample names, NanoAODTools imports and runs the `lepMVA` factory from `lepMVA_run3.py`.

Selection is based on the sample-name string, not JSON `year` or `INPUT_MODE`. Current sample names containing `2022` or `2023` map to `lepMVA`; other names fail early because the wrapper import path exposes only `lepMVA`. Supporting other Run 3 naming conventions requires a source-backed module-selection update. The selected cfg and `MATCH` determine which sample names reach this logic. The resolved summary records supported, selected, and intentionally rejected fallback module names.

### How to remove lepMVA from the processing chain

There is no supported user knob that disables lepMVA. Removing it requires a reviewed source change to the `-I CMGTools.NanoProc.tools.nanoAOD.lepMVA_run3 <module>` arguments assembled in `skimmer/skim_wrapper.py`. If `--module` is no longer needed, update its parser requirement and the `build_payload_command` call in `skimmer/lobster_config.py` in the same change so the wrapper contract remains consistent.

Do not remove lepMVA by editing sample cfg or JSON files. Sample names choose the module variant, but the processing chain itself is wired through the config-to-wrapper command path. Do not delete or edit the CMGTools module merely to skip it.

### Where lepMVA output branch names are defined

`lepMVA_run3.py` owns the output schema. The `lepMVA_run3` constructor defaults `suffix="run3"` and builds `self.branches` as `Electron_mvaTTH<suffix>` and `Muon_mvaTTH<suffix>`. The current `lepMVA` factory does not override that suffix, so it declares `Electron_mvaTTHrun3` and `Muon_mvaTTHrun3`. `beginFile` passes `self.branches` to `declareOutput`, and `analyze` fills matching keys. Neither `lobster_config.py` nor `skim_wrapper.py` passes a branch name or suffix.

### When to change lepMVA output branch names

Change these names only for an intentional output-schema change. Prefer changing the `suffix` supplied by the `lepMVA` factory in `lepMVA_run3.py`; if the branch expressions themselves change, keep `self.branches`, the `analyze` result keys, and object attributes consistent. Do not try to rename them in the Lobster control config. A rename can affect downstream plotting and analysis code, friend trees, selections, and compatibility with existing skim outputs. Update all downstream consumers and validate a small local output through a separately approved workflow before launching a campaign.

### Source checklist before changing lepMVA behavior

- `skimmer/lobster_config.py`: `select_module_name`, `build_payload_command`, and `Workflow(command=...)`.
- `skimmer/skim_wrapper.py`: required `--module` parsing and the `nano_postproc.py -I` argument list.
- `CMSSW_14_0_6/src/CMGTools/NanoProc/python/tools/nanoAOD/lepMVA_run3.py`: class, factory, model files, branch declarations, and branch fills.
- `CMSSW_14_0_6/src/CMGTools/NanoProc/python/tools/nanoAOD/friendVariableProducerTools.py`: `declareOutput` and `writeOutput` helpers used by the module.
- Downstream consumers of `Electron_mvaTTHrun3` and `Muon_mvaTTHrun3` before any schema change.

## skim_wrapper.py behavior

The wrapper parses one or more positional input files plus required `--cut`, required `--module`, optional `--out-dir` (default `.`), and optional `--nevents`.

It then:

1. Lists the task working directory.
2. Builds local basenames for every input.
3. Uses an already-local basename when present, including a basename found after replacing `file:` with an empty string.
4. Otherwise invokes `xrdcp -f <input> <local_name>` as fallback stage-in.
5. Runs `nano_postproc.py -c <cut> -I CMGTools.NanoProc.tools.nanoAOD.lepMVA_run3 <module> <out_dir> <local_files>`.
6. Adds `-N <nevents>` only when the optional wrapper argument is provided.
7. Lists the task working directory again.
8. Maps each input basename from `.root` to `_Skim.root`.
9. Invokes `haddnano.py output.root <skimmed_files>`.

The wrapper does not run `haddnano.py` through `python`; it calls `haddnano.py` directly. It does not validate that the expected `_Skim.root` files exist before merge. Recorded worker execution is limited to the two bounded cases below.

## Work Queue and Lobster runbook

These site-specific commands are operational examples, not evidence that the corresponding services, credentials, or paths are currently ready. Confirm them before use.

Run the factory and Lobster commands from the Run 3 `skimmer/` directory, where `lobster_config.py`, `factory_skim.json`, and `with_oasis_certs` are expected to be visible as working files. The documented `lobster process lobster_config.py` command depends on that working directory; do not rewrite it as `lobster process skimmer/lobster_config.py` without a separate source-backed change.

Preferred shell sequence:

1. Start from a fresh terminal session. A fresh terminal is preferred but not mandatory if the existing shell is clean.
2. Run the `preset_lobconda` cleanup before activating any conda/Lobster environment:

   ```bash
   unset PYTHONPATH
   unset PERL5LIB
   ```

3. Activate the appropriate conda/Lobster environment using the site-approved command:

   ```text
   <ACTIVATE_LOBSTER_CONDA_ENVIRONMENT>
   ```

4. Move to the Run 3 skimmer directory:

   ```bash
   cd /users/apiccine/work/LobSkim/LobsterSkimming/skimmer
   ```

5. Confirm the required files exist before launching:

   ```bash
   test -f factory_skim.json
   test -f with_oasis_certs
   ```

6. Confirm the Run 3 factory `--runos` is `al9-wq-7.11.1`, not the Run 2 CC7 value.

Run 3 Work Queue factory command:

```bash
nohup work_queue_factory -T condor -M lobster_apiccine.* -dall -o /tmp/apiccine_factoryskim.debug -C factory_skim.json -S /tmp/wq-apiccine-factoryskim --runos al9-wq-7.11.1 --wrapper ./with_oasis_certs --wrapper-input with_oasis_certs > /tmp/apiccine_factoryskim.log &
```

Factory command details:

- The command must run from inside the activated conda/Lobster environment.
- The command is backgrounded with `nohup` and `&`.
- Factory stdout is redirected to `/tmp/apiccine_factoryskim.log`.
- Factory debug output is written to `/tmp/apiccine_factoryskim.debug`.
- The worker project/sandbox directory is `/tmp/wq-apiccine-factoryskim`.
- The Work Queue master regex is `lobster_apiccine.*`.
- The factory profile is `factory_skim.json`.
- The wrapper file is `./with_oasis_certs`.
- The wrapper is also passed as wrapper input with `--wrapper-input with_oasis_certs`.

The factory is an external, operator-managed service. The Lobster config does
not start or stop it. The operator is responsible for its credentials,
resource profile, project-pattern match, logs, lifetime, and shutdown. Record
factory state separately from Lobster manager state.

Run 3 Lobster submission command:

```bash
lobster process lobster_config.py
```

Run the Lobster command from the same `skimmer/` directory assumption unless a future source-backed runbook revision documents a different launch directory.

Non-goals and safety caveats:

- Do not run `work_queue_factory`, `lobster process`, Condor, XRootD, DBS, setup, installer, CMSSW, or worker commands unless a separate production authorization explicitly permits those capabilities.
- Do not infer production readiness from these documented commands. The bounded validation record below covers two one-sample worker-backed `TESTING` cases, not setup, credentials, production storage, throughput, or a full campaign.
- Do not edit `factory_skim.json`, `with_oasis_certs`, wrapper code, setup scripts, sample cfgs, or nested `topeft` merely to match this runbook.
- Run 3 uses `--runos al9-wq-7.11.1`; Run 2 uses `--runos cc7-wq-7.11.1`.

`skimmer/factory_skim.json` is tracked in this repository and currently contains:

| field | value |
| --- | --- |
| `max-workers` | `500` |
| `min-workers` | `10` |
| `cores` | `4` |
| `memory` | `8000` |
| `disk` | `20000` |

Treat it as a source-controlled resource profile, not as a complete or validated factory launch command.

## Monitoring

Use the resolved project workdir when monitoring:

```bash
lobster status <workdir>
```

The operating record should identify `process.log`, `process.err`, `configure.log`, the resolved workdir and plotdir, factory logs, task exit codes, selected/failed/skipped counts, storage failures, worker resource exhaustion, and how campaign completion will be established. Submission acceptance is not campaign completion.

## Retry, recovery, and cleanup guidance

First diagnose the source of failures. Preserve the workdir database and its config. Do not reset, delete, or clean project state as routine recovery.

Request graceful shutdown with the resolved workdir:

```bash
lobster terminate <workdir>
```

The request is checkpoint-based. The manager notices it in the control loop,
stops creating new work, and exits at a later lifecycle checkpoint; tasks that
are already running may finish first. Confirm closure from `process.log`,
`process.err`, and process state before changing or removing project state.

The framework recovery note explains that `threshold_for_failure` controls unit retries and `threshold_for_skipping` controls file-access retries. For an existing project, changes belong in the workdir copy of `config.py` and must be applied through the supported Lobster configuration workflow. Recovery must define whether the manager is running, the current counters, the root-cause fix, new thresholds, resume/finalize intent, and verification. Never infer that cleanup is safe merely because processing stopped.

## Validation checklist before running

### static_preflight_checks

- Confirm Run 3 parent path, local branch `run3-mva-anpicci`, upstream target branch `run3_mva`, intended commit, and tracked state.
- Inspect nested `topeft` status separately; preserve its branch and local changes.
- Confirm only intended documentation/config edits are present and review their diff.

### configuration_source_checks

- Record all user-knob values and calculate normalized `INPUT_MODE`, normalized `TARGET`, `TYPE`, `YEAR`, `TAG`, and `CFG_NAME`.
- Confirm `CFG_NAME` exists under `topeft/input_samples/cfgs`.
- Confirm `sandbox_location`, wrapper name, workflow merge command, and output names.
- Review generated cut semantics and module mapping from source.

### environment_checks

- Confirm setup/build shell separation from conda/mamba.
- Record production Lobster environment, CMSSW environment, Python resolution, credentials, and site runtime.
- Confirm `CMSSW_14_0_6` and required CMGTools payload contents are complete.

### sample_metadata_checks

- Resolve the exact cfg and active JSON list after comments and `MATCH`.
- For `files`, verify every file string; for `dbs`, verify every `path`.
- Confirm selected sample names imply the intended module mapping.
- Confirm JSON `year`, `isData`, event counts, and weight metadata are plausible for the analysis, even though the current config does not validate all of them.

### path/output_checks

- Record `TAG`, master label, workdir, plotdir, logical output, local URL, and remote URL.
- Check for collisions, existing project state, quota, permissions, worker visibility, and endpoint policy.
- Confirm factory logs and monitoring destinations.

### production_run_checks

- Obtain explicit approval for factory, Lobster, Condor, network, storage, and worker execution.
- Record exact commands, host, working directory, external services, process IDs, and termination state.
- Validate selected workflows and outputs; do not equate static inspection, submission, or partial output with completion.

## Troubleshooting

| symptom | likely source-visible cause | first check |
| --- | --- | --- |
| Derived cfg is missing | `TYPE`, `YEAR`, or data cfg naming does not match available files. | Inspect `TYPE`, `YEAR`, `CFG_NAME`, and `topeft/input_samples/cfgs`. |
| Wrong campaign namespace | `TYPE` or `YEAR` changed without checking `TAG`. | Recalculate `TAG`, workdir, plotdir, and output path before launch. |
| No samples selected | `MATCH` excluded all active cfg lines, or the cfg only has commented entries. | Review active cfg lines and regexes. |
| DBS setup fails | Selected JSON lacks a usable `path` for `INPUT_MODE="dbs"`. | Read the selected JSON metadata. |
| Files setup fails | Selected JSON lacks a usable `files` list for `INPUT_MODE="files"`. | Read the selected JSON metadata and storage prefixes. |
| Wrong module | Sample name does not contain the expected `2022` or `2023` marker. | Check `select_module_name` behavior and selected sample names. |
| Sandbox failure | Missing/incomplete CMSSW release or incompatible environment. | Verify exact sandbox path and payload contents; do not rerun setup blindly. |
| Stage-in failure | Input not local and fallback XRootD transfer fails. | Check resolved input URL, credentials, endpoint, and wrapper log. |
| Local stage-out failure | Worker cannot see local destination or its parent is missing. | Check `file://` path visibility and remote fallback. |
| Remote stage-out failure | Endpoint, credentials, permissions, or quota issue. | Inspect task transfer errors and destination policy. |
| Workdir collision | Same step/target/type/year/day reused. | Compare derived non-test workdir before launch. |
| Repeated failure threshold | Root cause persists or counters reached limits. | Use framework recovery guidance; do not delete workdir state. |
| Config works statically but campaign fails | Static inspection cannot prove services, environment, or payload runtime. | Perform a separately authorized bounded validation with exact records. |

## Framework references

Consult the local Lobster framework clone for framework-level details:

- `/users/apiccine/work/LobSkim/lobster/docs/config.rst`
- `/users/apiccine/work/LobSkim/lobster/docs/run.rst`
- `/users/apiccine/work/LobSkim/lobster/docs/monitor.rst`
- `/users/apiccine/work/LobSkim/lobster/docs/trouble.rst`
- `/users/apiccine/work/LobSkim/lobster/docs/FailureThresholdRecovery.md`
- `/users/apiccine/work/LobSkim/lobster/docs/operational_guidance.md`
- `/users/apiccine/work/LobSkim/lobster/docs/StorageFileURLs.md`
- `/users/apiccine/work/LobSkim/lobster/examples/factory.json`

Some framework examples are historical or site-specific. This tutorial is the source of truth for current Run 3 branch-specific values; neither source supersedes an explicit current production authorization.

## Bounded Run 3 validation state

A sequence of bounded checks established the following scoped evidence for the current Run 3 configuration. These results are validation records, not a production campaign certification.

| behavior | evidence and scope | result |
| --- | --- | --- |
| Canonical cfg, input-mode, and lepMVA resolution | Static helper checks used representative Run 3 metadata. | Source-defined behavior passed the bounded matrix checks. |
| Explicit files runtime preparation | One isolated `TESTING=True` signal config used real files metadata with no workers. | Lobster prepared the files-backed workflow successfully. |
| Explicit DBS runtime preparation | One isolated `TESTING=True` data config used a single-dataset DBS query with no workers. | Lobster prepared the DBS-backed workflow successfully. |
| Auto DBS worker execution | The check selected only `EGamma_0_V1_Run2023D-22Sep2023` with `TYPE="data"`, `YEAR="2023BPix"`, and `INPUT_MODE="auto"`. | Auto resolved to DBS; a worker task completed successfully and produced a ROOT output in the dedicated validation area. |
| Auto files worker execution | The check selected only `tHq_2022EE` with `TYPE="signal"`, `YEAR="2022EE"`, and `INPUT_MODE="auto"`. | Auto resolved to files; worker tasks completed successfully and produced ROOT outputs in the dedicated validation area. |
| Stable validation config paths | Each generated config was imported twice in isolated processes and its state-relevant values were compared. | Tag/label, workdir, plotdir, output, storage output, selected sample, and expected `config.pkl` location were stable across re-imports. |
| Terminate by config | `lobster terminate <generated_config>` was invoked after task-completion/output evidence. | The first invocation returned zero for both cases; both managers exited gracefully without a local OS signal. |

The configs were temporary validation copies, each selected exactly one sample, and all active workdir, plot, output, and `file://` stageout paths were under a dedicated validation directory. Workers came from an already-running, operator-managed `work_queue_factory`. The factory was started and managed separately; the validation only used its available workers. Work Queue manager catalog advertisement occurred because no source-supported disable option was found.

For generated `TESTING` configs that must later be passed to `lobster terminate <config>`, keep every state-lookup value deterministic across imports. In particular, freeze the tag or label and the derived workdir, plotdir, output, storage output, and `config.pkl` location. Verify those resolved values in isolated re-imports before processing. A config that recomputes timestamped paths on import can cause the terminate command to look in a different workdir.

Termination is checkpoint-based rather than instantaneous. Terminate was requested shortly after the first success evidence, but the managers exited at a later lifecycle checkpoint and already-running tasks were allowed to finish.

A representative manual inspection of auto-DBS output `output_10.root` found `Events`, `LuminosityBlocks`, and `Runs` trees; `Events` had 7,230 entries; and the `Electron_mvaTTHrun3` and `Muon_mvaTTHrun3` branches were present. This was not an automated content scan of every output or a physics-correctness validation.

An earlier attempt exposed timestamp-dependent generated-config lookup: re-importing a config could resolve a different workdir and prevent config-based termination. The final checks used fixed state paths in isolated configs, confirmed graceful manager closure, and demonstrated successful config-based termination. Repository source was not changed to implement this validation setup.

## Known caveats and intentionally unresolved items

- Fresh setup selects `run3_test_mmerged`, while the observed nested checkout is `run3_test_mmerged_anpicci`; no reconciliation is prescribed here.
- Setup/install bodies, production storage, credentials, and campaign-scale services remain unvalidated.
- Work Queue factory and Lobster submission commands are site-specific and require operator review before use.
- Canonical cfg derivation intentionally does not select alternate `NDSkim_*`, CR/SR-specific, central-signal, signal-subtype, synchronization, or general cfg names. Users should not expect `TYPE` and `YEAR` to discover those files automatically; support requires a separate explicit source decision and selection design rather than a broad filename search.
- Auto DBS and auto files each passed one-sample worker-backed `TESTING` validation. No full production or broad multi-sample throughput campaign was run. Mixed storage remains source-defined/static-validated only because no real same-context mixed candidate was found; mixed files/DBS `auto` requires a separately approved runtime validation.
- Run 3 lepMVA selection is sample-name based. Current names containing `2022` or `2023` map to `lepMVA`; other names fail early because the wrapper import path exposes only `lepMVA`. Supporting other Run 3 sample-name conventions requires a source-backed module-selection update.
- `skim_wrapper.py` requires runtime tools such as `xrdcp`, `nano_postproc.py`, `CMGTools.NanoProc.tools.nanoAOD.lepMVA_run3`, and `haddnano.py`. The worker payload was exercised only in the two bounded one-sample cases.
- Physics correctness and final analysis-level event-content validation remain out of scope. The representative ROOT inspection does not replace an automated all-output content audit.
- Graceful termination is checkpoint-based; tasks already running when termination is requested may finish before the manager exits.
- The setup script's existing-check logic does not validate checkout branch or installation completeness.
