# Run 2 LobSkim tutorial

## 1. Scope and assumptions

This tutorial documents the Run 2 `LobsterSkimmingRun2` repository, its independent nested `topeft` checkout, and its integration with the Lobster framework. The examples assume a workspace containing sibling `LobsterSkimmingRun2`, `LobsterSkimming`, and `lobster` checkouts. Set a reusable workspace variable before following path examples:

```bash
export LOBSKIM_WORKSPACE=<your-work-area>/LobSkim
```

The source descriptions below are based on static inspection. The setup installers, configuration module, wrapper, Lobster commands, Work Queue, CMSSW build, network access, and production campaign have no recorded Run 2 runtime validation. Treat production behavior as unvalidated until bounded evidence is recorded.

Run 2 and Run 3 use different parent repositories, CMSSW releases, payloads, and campaign assumptions. Do not copy Run 3 setup values into this repository merely because the directory layouts look similar.

## 2. Repository layout

| path | ownership and purpose |
| --- | --- |
| `README.md` | Concise Run 2 landing page. |
| `docs/run2_lobskim_tutorial.md` | This repository-specific operational reference. |
| `setup.py` | No-argument setup entry point; dispatches the two installers when their targets are absent. |
| `scripts/install_configs.sh` | Sparse-checkout helper for cfg and sample-JSON paths. |
| `scripts/install_cmssw.sh` | Builds the Run 2 CMSSW/payload area when intentionally executed. |
| `skimmer/lobster_config.py` | User knobs, validation, workflow construction, storage, paths, and Lobster `Config`. |
| `skimmer/skim_wrapper.py` | Python-2-style worker payload for stage-in, NanoAOD postprocessing, and output merge. |
| `skimmer/tools/utils.py` | Reads cfg entries and referenced sample JSON files. |
| `postmortem/` | Separate Run 2 post-mortem EFT reweighting material with incomplete provenance/validation. |
| `topeft/` | Independent nested Git repository supplying cfg and JSON metadata; protected from parent-repository cleanup. |

The workspace may also contain a separate Run 3 parent at `${LOBSKIM_WORKSPACE}/LobsterSkimming`. Confirm the repository path before editing or running anything.

## 3. Protected nested topeft repository

`topeft/` is a complete, independent Git repository. The Run 2 parent may display it as untracked because it is not a submodule. That output does not make it disposable.

- Run Git inspection for `topeft` from `LobsterSkimmingRun2/topeft`, not from the parent.
- Do not clean or delete it from the parent.
- Do not stage it into the parent.
- Do not assume authorization to edit the parent also permits editing `topeft`, or vice versa.
- Preserve an existing checkout even when `setup.py` would choose a different branch for a fresh install.

`setup.py` only checks whether the `topeft` path exists. If it exists, setup skips installation and does not verify its branch or content.

## 4. Environment policy

Keep setup/build and Lobster-control concerns explicit:

1. CMSSW release creation and compilation belong in a clean shell outside conda/mamba.
2. The Lobster control process needs the site-approved Lobster environment and a compatible CMSSW environment in the order required by the local installation.
3. The worker sandbox is fixed by `sandbox_location` in `skimmer/lobster_config.py`; it is not inferred from whichever shell happens to be active.
4. Do not mutate a shared environment to repair a campaign. Record the missing dependency and handle environment maintenance separately.

The manager environment does not prove which interpreter the worker payload resolves inside its sandbox. Record both environments separately when validating a campaign.

## 5. Setup overview

`setup.py` has a no-argument entry point and immediately calls `main()` when executed by a user. It finds the parent Git top level, then:

- skips `topeft` when that path already exists;
- otherwise calls `scripts/install_configs.sh` for the two sparse paths;
- skips the CMSSW installation when `CMSSW_10_6_19_patch2` already exists;
- otherwise calls `scripts/install_cmssw.sh` with the parent path, release, and SCRAM architecture.

Both installers use `set -euo pipefail` and positional-argument checks. This describes source-visible behavior only; no GitHub, CVMFS, SCRAM, or storage execution is recorded for the current setup.

## 6. Conda/mamba environment

Do not run CMSSW creation or compilation from an active conda/mamba environment. Use a separate clean shell for setup. When controlling Lobster, activate only the established Lobster environment; do not install or upgrade packages ad hoc during a campaign.

Before production, record the exact activation sequence used at the site, the resulting Python executable, and relevant environment variables. This repository does not prescribe an environment name; examples use `conda activate <your-lobster-env>` as a placeholder.

## 7. CMSSW setup

The Run 2 setup values are:

| item | current source value |
| --- | --- |
| release | `CMSSW_10_6_19_patch2` |
| SCRAM architecture | `slc7_amd64_gcc700` |
| CVMFS setup root | `/cvmfs/cms.cern.ch` |
| CMGTools repository | `git@github.com:sscruz/cmgtools-lite.git` |
| CMGTools branch | `104X_dev_nano_lepMVA` |
| NanoAODTools repository | `https://github.com/cms-nanoAOD/nanoAOD-tools.git` |

The installer refuses to continue when `CMSSW_BASE` is already set, sources the CMS defaults, exports the Run 2 architecture, creates the release, clones the two payload repositories under `src`, evaluates the release runtime, and builds it. These operations are network- and build-heavy and require separate approval/validation. Do not treat the presence of a directory as proof that the release or payload is complete.

`skimmer/lobster_config.py` independently fixes the worker sandbox to `<Run2-parent>/CMSSW_10_6_19_patch2`; if the release location changes, setup and config must be reviewed together.

## 8. topeft sample cfg/json setup

Fresh setup uses:

| item | current value |
| --- | --- |
| repository | `https://github.com/TopEFT/topeft.git` |
| branch/tag | `run3_test_mmerged` |
| sparse path 1 | `input_samples/cfgs` |
| sparse path 2 | `input_samples/sample_jsons` |

These paths resolve under the independent checkout as `topeft/input_samples/cfgs` and `topeft/input_samples/sample_jsons`. A campaign cfg is resolved as `topeft/input_samples/cfgs/<CFG_NAME>`.

`read_cfg` ignores blank/commented cfg lines, treats `root:` lines as cfg metadata, applies every regex in `MATCH`, resolves JSON paths relative to the cfg directory, and loads each selected JSON. The sample key is the JSON basename without `.json`.

## 9. Source files used by a skim

A campaign is assembled from:

1. `skimmer/lobster_config.py` for campaign knobs and workflow objects.
2. `topeft/input_samples/cfgs/<CFG_NAME>` for the active JSON list.
3. Each referenced JSON for `year`, `path` and/or `files`, plus analysis metadata.
4. `skimmer/skim_wrapper.py` as the worker command.
5. `CMSSW_10_6_19_patch2` as the Lobster sandbox.
6. `CMSSW_10_6_19_patch2/src/PhysicsTools/NanoAODTools/scripts/haddnano.py` as an extra input and merge helper.

## 10. Quick start checklist

Before requesting a production launch:

1. Confirm the parent is the Run 2 clone and inspect its branch/status without altering it.
2. Inspect nested `topeft` separately; record its branch, commit, and dirty state.
3. Confirm the expected CMSSW release and payload paths exist; do not infer completeness from names alone.
4. Edit only the intended user knobs in `skimmer/lobster_config.py`.
5. Verify `TYPE` and `CFG_NAME` identify the same campaign kind.
6. Review every selected JSON for a supported `year` and the metadata required by `INPUT_MODE`.
7. Calculate the derived label, workdir, plotdir, and output path before launch; ensure they do not collide with an existing campaign.
8. Obtain user-approved Work Queue and Lobster commands; the placeholders below are deliberately non-runnable.
9. Record credentials, storage reachability, quotas, factory resources, logs, and shutdown procedure.
10. Run production only after explicit authorization and retain the exact command/output record.

## 11. Editing skimmer/lobster_config.py

The `USER KNOBS` block is the intended campaign-editing surface. Avoid changing helpers, workflow construction, storage logic, or wrapper code just to launch a different sample set. After edits, review the complete resolved configuration output produced when Lobster eventually evaluates the config; do not rely only on filenames.

For common goals:

- To switch campaign kind, change `TYPE` and set `CFG_NAME` to the matching data, background, or signal cfg.
- To limit a cfg to one Run 2 period, set `YEAR` to `2016APV`, `2016`, `2017`, or `2018`; leave it empty to retain all supported years.
- To narrow a retry or bounded check, edit `MATCH` with regexes that match cfg lines.
- To require explicit file lists, set `INPUT_MODE="files"` and verify every selected JSON has a usable `files` list.
- To require DBS discovery, set `INPUT_MODE="dbs"` and verify every selected JSON has a valid CMS dataset `path`.
- To choose per-sample discovery, set `INPUT_MODE="auto"` and review the resolved DBS/files counts and storage profile.
- To run SR-like skim cuts, change `TARGET` to `"SR"`; leave `TARGET="CR"` for the looser control-region cut.

## 12. User knobs reference

| knob | default | allowed_or_expected_values | what_it_controls | when_to_edit | risk_or_caveat | source_file |
| --- | --- | --- | --- | --- | --- | --- |
| `TESTING` | `False` | Boolean | Switches labels and paths to timestamped `test/lobster_skimtest_*` locations. | For a deliberately bounded test campaign. | A test-labelled path does not make submission safe or lightweight. Timestamp-derived paths can change on a later config import. | `skimmer/lobster_config.py` |
| `INPUT_MODE` | `"auto"` | `dbs`, `files`, `auto` (case-normalized) | Per-sample dataset type, storage profile, and `xrootd_servers`. | When all samples must use one metadata contract or automatic selection is desired. | Mixed auto campaigns select `storage_files` if any sample uses files; mixed behavior lacks production validation here. | same |
| `PROTOCOL_LOCAL` | `"file://"` | Storage protocol understood by Lobster; current design expects `file://`. | Local input/output prefix. | Only for an intentional storage-profile change. | `file://` requires worker-visible paths and existing destination parents. | same |
| `PROTOCOL_REMOTE` | `"root://"` | Current design expects `root://`. | Remote XRootD prefix. | Only when changing remote transport. | Endpoint reachability was not validated. | same |
| `TARGET` | `"CR"` | `SR`, `CR` (case-normalized) | Changes the skim cut and path/label segment. | Choose signal-region or control-region skim semantics. | `SR` adds the opposite-sign two-lepton veto; `CR` does not. | same |
| `YEAR` | `""` | empty, `2016APV`, `2016`, `2017`, `2018` | Optional filter after cfg/JSON loading. | Limit a cfg to one Run 2 year. | `multi-year` is rejected; every selected JSON still needs valid `year`. | same |
| `STEP` | `"skimmed"` | Nonempty campaign label chosen by operator. | Label/path segment. | When naming a distinct processing step. | Changing it creates different paths; it does not change physics behavior. | same |
| `TYPE` | `"background"` | `data`, `background`, `signal` (case-normalized) | `TAG` and output identity. | Match the selected cfg kind. | It does not change the skim cut. A mismatch is rejected by default. | same |
| `CFG_NAME` | `"mc_background_samples.cfg"` | Existing cfg basename under `topeft/input_samples/cfgs`; name must identify data/background/signal. | Selects the sample list. | For every campaign sample-set change. | Unknown or ambiguous kind tokens are always rejected. | same |
| `ALLOW_TYPE_CFG_MISMATCH` | `False` | Boolean | Allows a known, unambiguous cfg kind to differ from `TYPE`. | Only for a reviewed exceptional labeling choice. | Does not allow unknown/ambiguous cfg names and can make output identity misleading. | same |
| `MATCH` | `[]` | List of regular-expression strings | Filters cfg JSON lines; empty selects all active lines. | For a narrow subset or retry. | Regex matches the cfg line; verify selected count and names. | same |
| `WRAPPER` | `"skim_wrapper.py"` | Worker wrapper present in the skimmer working context. | Payload command and `extra_inputs`. | Only when intentionally replacing the worker contract. | Wrapper changes require separate validation. | same |
| `SRC_REMOTE` | `"cms-xrd-global.cern.ch"` | XRootD host, without protocol. | Remote stage-in prefix and DBS `xrootd_servers`. | When site/storage policy changes. | Network access and fallback order must be validated. | same |
| `SRC_LOCAL` | `"/cms/cephfs/data"` | Worker-visible absolute path. | Local file-input prefix. | When local storage mount changes. | Must exist on workers; no check is performed statically. | same |
| `DST_REMOTE` | `"cmsxrootd.crc.nd.edu"` | XRootD destination host. | Remote stage-out prefix. | When output endpoint changes. | Credentials, quota, and permissions must be checked. | same |
| `DST_LOCAL` | `"/cms/cephfs/data"` | Worker-visible absolute path. | Local stage-out prefix. | When local destination mount changes. | Destination parents must already be reachable for `file://`. | same |
| `WORKDIR_BASE` | `"/tmpscratch/users/$USER"` | Writable manager-local base path. | Lobster state/database directory. | When project state belongs on another filesystem. | Never reuse a workdir for an unrelated campaign. Shell variable expansion behavior must be confirmed in the actual Lobster environment. | same |

## 13. Derived labels and paths

`CAMPAIGN_TYPE` is normalized from `TYPE`. `TSTAMP1` is `YYYYMMDD_HHMM`; `ver` is `vYYMMDD`.

| label_or_path | source_expression_summary | default_or_pattern | meaning | when_it_changes | risk_or_caveat |
| --- | --- | --- | --- | --- | --- |
| `TAG` | `<CAMPAIGN_TYPE>/NAOD_ULv9_lepMVA-run2` | `background/NAOD_ULv9_lepMVA-run2` | Campaign/output namespace. | `TYPE`. | Does not include `YEAR`, `TARGET`, or input mode. |
| `master_label` | `<STEP>_<TARGET>_<INPUT_MODE>_lobPY3_<timestamp>` | `skimmed_CR_auto_lobPY3_YYYYMMDD_HHMM` | Lobster config label. | Step, target, mode, evaluation minute; test has `testlobPY3`. | The `lobPY3` label does not prove the Python version used by the wrapper. |
| `sandbox_location` | `<git-top>/CMSSW_10_6_19_patch2` | Absolute Run 2 release path | Worker CMSSW sandbox. | Parent checkout location or source edit. | Must match a valid built release. |
| `cfg_fpath` | `<git-top>/topeft/input_samples/cfgs/<CFG_NAME>` | Default ends in `mc_background_samples.cfg` | Selected cfg. | `CFG_NAME` or checkout path. | Nested checkout ownership applies. |
| `workdir_path` | `<WORKDIR_BASE>/<STEP>_<TARGET>/<TAG>/vYYMMDD` | `/tmpscratch/users/$USER/skimmed_CR/background/NAOD_ULv9_lepMVA-run2/vYYMMDD` | Manager state. | Base, step, target, type, day; test uses timestamped test path. | Same-day repeated non-test launches collide unless intentionally resumed. |
| `plotdir_path` | `~/www/lobster/<STEP>_<TARGET>/<TAG>/vYYMMDD` | Tilde-based web path | Monitoring plots. | Step, target, type, day; test uses timestamp. | Confirm tilde handling and directory permissions. |
| `output_path` | `/store/user/$USER/<STEP>_<TARGET>/<TAG>/vYYMMDD` | Default Run 2 stage-out namespace | Logical output path appended to destination prefixes. | Step, target, type, day; test uses timestamp. | Confirm `$USER`, quota, permissions, and non-collision. |
| local output URL | `file://` + local root + `//` + output path | `file:///cms/cephfs/data///store/user/$USER/...` | First files-mode stage-out method. | Protocol/root/output. | Repeated slashes are source-visible; worker filesystem must be accessible. |
| remote output URL | `root://` + remote host + `//` + output path | `root://cmsxrootd.crc.nd.edu///store/user/$USER/...` | DBS output and files-mode fallback. | Protocol/host/output. | Not network-validated here. |

## 14. TESTING paths and config-based control

`TESTING=True` changes the label, workdir, plotdir, and output path to
minute-stamped test locations. It does not disable submission, workers,
network access, DBS access, or storage writes.

Lobster commands can identify a project by its workdir or by importing the
original config and reading `Config.workdir`. For a generated config that will
later be passed to `lobster terminate <config.py>`, freeze every state lookup
value across imports, including the label, workdir, plotdir, output path,
storage outputs, and expected `config.pkl` location. Verify those values in a
separate import before launch. If the config can recompute a different
timestamp, use the known workdir for `lobster status`, `lobster configure`, and
`lobster terminate` instead.

## 15. TYPE and CFG_NAME consistency

`TYPE` accepts only `data`, `background`, or `signal`. The cfg basename is tokenized on non-alphanumeric characters:

- token `data` classifies data;
- `background` or `bkg` classifies background;
- `signal` or `sig` classifies signal.

No marker is `unknown`; multiple kinds are `ambiguous`. Both are rejected even when `ALLOW_TYPE_CFG_MISMATCH=True`. A known mismatch is rejected unless that override is explicitly enabled. Validation occurs before `TAG`, workdir, plotdir, and output paths are derived.

For common changes: use `TYPE="signal"` with `mc_signal_samples.cfg`, `TYPE="data"` with `data_samples.cfg`, and retain `TYPE="background"` with `mc_background_samples.cfg` unless a reviewed cfg says otherwise.

## 16. TARGET and skim-cut behavior

Only `SR` and `CR` are valid. Both require at least two selected electrons, muons, and taus in total. Electron, muon, and tau object requirements are encoded directly in `build_skim_cut`.

`SR` additionally vetoes events with exactly two selected leptons whose summed charge is zero. `CR` uses only the at-least-two requirement. `TYPE` does not affect this cut. The generated expression has whitespace removed and is checked for balanced parentheses and accidental literal outer quotes.

## 17. YEAR filtering and sample-year metadata

Allowed normalized sample years are `2016APV`, `2016`, `2017`, and `2018`. Every active selected JSON must contain a nonempty string `year` with one of those values, even when `YEAR=""`.

- Empty `YEAR` means no year filter, not an unspecified module.
- A concrete `YEAR` filters after all active cfg JSONs are loaded and validated.
- A filter that removes every sample raises an error.
- `YEAR="multi-year"` is explicitly rejected; use empty `YEAR` for all supported years.

The wrapper module is selected per sample: `lepMVA_2016APV`, `lepMVA_2016`, `lepMVA_2017`, or `lepMVA_2018`. Review the printed before/after year counts when configuration is eventually evaluated in an authorized run.

## 18. INPUT_MODE and sample input metadata

| mode | required metadata | dataset object | behavior |
| --- | --- | --- | --- |
| `files` | Nonempty `files` list containing nonempty strings. | Lobster `Dataset(files=..., files_per_task=1, patterns=["*.root"])` | Every sample uses explicit files. |
| `dbs` | `path` in `/PrimaryDataset/ProcessedDataset/DataTier` form; tier must be a recognized CMS tier. | `cmssw.Dataset(dataset=..., lumis_per_task=1, file_based=True)` | Every sample uses DBS-style discovery and AAA server configuration. |
| `auto` | Prefer a valid DBS `path`; otherwise require usable `files`. | Per-sample choice | Supports source-visible mixed decisions. |

`/store/...`, `root://...`, `file://...`, and arbitrary absolute paths are not DBS dataset names. Valid tier endings include `NANOAOD`, `NANOAODSIM`, `MINIAOD`, `MINIAODSIM`, `AOD`, `AODSIM`, and `USER`.

In addition to input metadata, each selected JSON requires `year`. Other fields such as cross section, data flags, event counts, and tree name are owned by the analysis metadata but are not consumed by this config's input-mode decision.

## 19. Storage and XRootD behavior

`storage_dbs` has only the remote output and leaves input streaming enabled. `storage_files` tries local then remote input prefixes and local then remote output prefixes, also with streaming enabled.

Selection is:

- global `files`: `storage_files`;
- global `dbs`: `storage_dbs`;
- `auto`: `storage_files` if at least one selected workflow uses explicit files, otherwise `storage_dbs`.

If at least one workflow uses DBS, `AdvancedOptions.xrootd_servers` becomes `["cms-xrd-global.cern.ch"]`. The historical Lobster troubleshooting docs caution about mixed DBS/local datasets, while this current config contains explicit per-sample auto logic. Treat mixed-mode production as a current source-visible intent that still needs campaign validation.

At worker level, `skim_wrapper.py` normalizes `file:` strings to `file://`, checks whether each basename already exists locally, and otherwise invokes `xrdcp -f` as a last-resort stage-in. No real Run 2 XRootD validation is recorded.

## 20. Workflow and AdvancedOptions behavior

Each selected sample becomes one workflow with hyphens replaced by underscores in its label.

| field | current value or rule | operational meaning |
| --- | --- | --- |
| category | `processing` | Shared resource category. |
| category resources | 1 core, 1500 MB memory, 4500 MB disk | Initial per-task request. |
| sandbox | `CMSSW_10_6_19_patch2` | Worker software release. |
| task splitting | one lumi per task for DBS; one file per task for files | Dataset-dependent unit construction. |
| `extra_inputs` | wrapper plus NanoAODTools `haddnano.py` | Files shipped with tasks. |
| `outputs` | `output.root` | Required wrapper product. |
| command | `python <wrapper> --cut <cut> --module <module> --out-dir . @inputfiles` | Actual Lobster payload string. |
| `merge_command` | `haddnano.py @outputfiles @inputfiles` | Lobster-level output merge. |
| `merge_size` | `537M` | Merge target size. |
| `globaltag` | `False` | No CMSSW global-tag setup requested by workflow. |
| `cleanup_input` | `False` | Input cleanup disabled. |

`AdvancedOptions` currently sets `bad_exit_codes=[127, 160]`, `log_level=1`, `payload=10`, `osg_version="3.6"`, `threshold_for_failure=10`, and `threshold_for_skipping=10`, plus `xrootd_servers` when DBS inputs exist. Changing thresholds affects retry behavior and should be based on diagnosed failures, not used to hide a persistent root cause.

## 21. Wrapper command construction

| argument_or_behavior | source_file | Run_2_behavior | user_should_edit | risk_or_caveat |
| --- | --- | --- | --- | --- |
| interpreter | `lobster_config.py` | Payload starts with `python`. | No, unless separately validating the sandbox/runtime contract. | `skim_wrapper.py` uses Python-2 print syntax. |
| wrapper | both | Default `skim_wrapper.py`; shipped as `extra_inputs`. | Normally no. | Replacement changes the worker contract. |
| `--cut` | config helper | Receives whitespace-free generated skim expression. | Change `TARGET`, not the constructed command. | Shell quoting is intentionally not inserted into the actual Lobster command. |
| `--module` | config helper | Per-sample `lepMVA_<year>`. | Change/fix JSON `year`, not command text. | Wrong year selects wrong module. |
| `--out-dir` | config helper | `.` | Normally no. | Wrapper and expected output assume task working directory. |
| `@inputfiles` | Lobster substitution | Expanded to task inputs. | No. | Expansion semantics belong to Lobster. |
| `output.root` | wrapper/workflow | Final required task output. | No. | Missing output causes task failure/stage-out failure. |
| workflow merge | config | `haddnano.py @outputfiles @inputfiles`. | Only in a separately reviewed workflow change. | Distinct from wrapper's per-task merge. |

The config also creates a `shlex.join` display command for logs. That display form is not passed to Lobster because parser handling of shell quoting has not been established.

## 22. lepMVA module wiring and output branches

### Where lepMVA modules are defined

The Run 2 installer places the CMGTools payload in the `CMSSW_10_6_19_patch2` sandbox by cloning the `104X_dev_nano_lepMVA` branch of `sscruz/cmgtools-lite` in `scripts/install_cmssw.sh`. In the current local payload:

- `CMSSW_10_6_19_patch2/src/CMGTools/TTHAnalysis/python/tools/nanoAOD/LepMVAULFriend.py` defines `LepMVAFriend` and the factories `lepMVA_2016`, `lepMVA_2016APV`, `lepMVA_2017`, and `lepMVA_2018`.
- `CMSSW_10_6_19_patch2/src/CMGTools/TTHAnalysis/python/tools/nanoAOD/ttH_modules.py` imports those four factories, making them available through `CMGTools.TTHAnalysis.tools.nanoAOD.ttH_modules` when `nano_postproc.py` processes `-I`.

The factories are selected by name when `nano_postproc.py` starts the postprocessing job; they are not defined by sample cfg or JSON metadata.

### How lepMVA is selected and wired into nano_postproc.py

The wiring is per sample and follows this path:

1. `select_module_name` in `skimmer/lobster_config.py` maps the normalized JSON `year` to `lepMVA_2016APV`, `lepMVA_2016`, `lepMVA_2017`, or `lepMVA_2018`.
2. `build_payload_command` adds `--module <module_name>` to the exact string assigned to `Workflow(command=command_string)`.
3. `skimmer/skim_wrapper.py` parses `--module` and constructs `nano_postproc.py -I CMGTools.TTHAnalysis.tools.nanoAOD.ttH_modules lepJetBTagDeepFlav,<module>`.
4. NanoAODTools imports the named factories from `ttH_modules.py` and runs `lepJetBTagDeepFlav` followed by the selected lepMVA factory.

Module selection is year-based, independent of `INPUT_MODE`, and performed separately for every selected sample. The `YEAR` campaign filter can reduce the selected samples, but each remaining sample's JSON `year` still determines its module.

### How to remove lepMVA from the processing chain

There is no supported user knob that disables lepMVA. Removing it requires a reviewed source change to the module list assembled in `skimmer/skim_wrapper.py`: remove the selected `<module>` from `lepJetBTagDeepFlav,<module>` while deciding explicitly whether `lepJetBTagDeepFlav` remains. If `--module` is no longer needed, update its parser requirement and the `build_payload_command` call in `skimmer/lobster_config.py` in the same change so the wrapper contract remains consistent.

Do not remove lepMVA by editing sample cfg or JSON files. Those files select samples and provide the year used to choose a variant; the processing chain itself is wired through the config-to-wrapper command path. Do not delete or edit CMGTools files merely to skip the module.

### Where lepMVA output branch names are defined

`LepMVAULFriend.py` owns the output schema. `LepMVAFriend.beginFile` declares `Electron_mvaTTHUL` and `Muon_mvaTTHUL` from the format string `<collection>_mvaTTHUL`; `analyze` fills those same names. Neither `lobster_config.py` nor `skim_wrapper.py` passes an output branch name.

### When to change lepMVA output branch names

Change these names only for an intentional output-schema change. Make the declaration and fill keys agree in `LepMVAULFriend.py`; do not try to rename them in the Lobster control config. A rename can affect downstream plotting and analysis code, friend trees, selections, and compatibility with existing skim outputs. Update all downstream consumers and validate a small local output through a separately approved workflow before launching a campaign.

### Source checklist before changing lepMVA behavior

- `skimmer/lobster_config.py`: `select_module_name`, `build_payload_command`, and `Workflow(command=...)`.
- `skimmer/skim_wrapper.py`: `--module` parsing and the `nano_postproc.py -I` argument list.
- `CMSSW_10_6_19_patch2/src/CMGTools/TTHAnalysis/python/tools/nanoAOD/ttH_modules.py`: imported module factories and the additional `lepJetBTagDeepFlav` module.
- `CMSSW_10_6_19_patch2/src/CMGTools/TTHAnalysis/python/tools/nanoAOD/LepMVAULFriend.py`: factory variants, branch declarations, and branch fills.
- Downstream consumers of `Electron_mvaTTHUL` and `Muon_mvaTTHUL` before any schema change.

## 23. skim_wrapper.py behavior

The wrapper parses one or more positional input files plus `--cut`, `--module`, and optional `--out-dir` (default `.`). It:

1. Lists the task working directory.
2. Converts `file:` occurrences to `file://`.
3. Uses an already-local basename when present; otherwise performs fallback `xrdcp -f` stage-in.
4. Sleeps for 10 seconds.
5. Runs `nano_postproc.py` with the cut and imports `CMGTools.TTHAnalysis.tools.nanoAOD.ttH_modules`, selecting `lepJetBTagDeepFlav,<module>`.
6. Sleeps again and lists the working directory.
7. Maps each input basename from `.root` to `_Skim.root` and invokes `python haddnano.py output.root ...`.

The script is Python-2-style source. It does not validate that `--cut` or `--module` was supplied before constructing the command, and output filename derivation assumes `.root` inputs. No recorded Run 2 worker execution establishes this behavior end to end.

## 24. Work Queue and Lobster runbook

These site-specific commands are operational examples, not evidence that the corresponding services, credentials, or paths are currently ready. Confirm them before use.

Run the factory and Lobster commands from the Run 2 `skimmer/` directory, where `lobster_config.py`, `factory_skim.json`, and `with_oasis_certs` are expected to be visible as working files. The documented `lobster process lobster_config.py` command depends on that working directory; do not rewrite it as `lobster process skimmer/lobster_config.py` without a separate source-backed change.

Preferred shell sequence:

1. Start from a fresh terminal session. A fresh terminal is preferred but not mandatory if the existing shell is clean.
2. Run the `preset_lobconda` cleanup before activating any conda/Lobster environment:

   ```bash
   unset PYTHONPATH
   unset PERL5LIB
   ```

3. Activate the appropriate conda/Lobster environment using the site-approved command:

   ```bash
   conda activate <your-lobster-env>
   ```

4. Move to the Run 2 skimmer directory:

   ```bash
   cd ${LOBSKIM_WORKSPACE}/LobsterSkimmingRun2/skimmer
   ```

5. Confirm the required files exist before launching:

   ```bash
   test -f factory_skim.json
   test -f with_oasis_certs
   ```

6. Confirm the Run 2 factory `--runos` is `cc7-wq-7.11.1`, not the Run 3 AL9 value.

Run 2 Work Queue factory command:

```bash
nohup work_queue_factory -T condor -M "lobster_${USER}.*" -dall -o /tmp/${USER}_factoryskim.debug -C factory_skim.json -S /tmp/wq-${USER}-factoryskim --runos cc7-wq-7.11.1 --wrapper ./with_oasis_certs --wrapper-input with_oasis_certs > /tmp/${USER}_factoryskim.log &
```

Factory command details:

- The command must run from inside the activated conda/Lobster environment.
- The command is backgrounded with `nohup` and `&`.
- Factory stdout is redirected to `/tmp/${USER}_factoryskim.log`.
- Factory debug output is written to `/tmp/${USER}_factoryskim.debug`.
- The worker project/sandbox directory is `/tmp/wq-${USER}-factoryskim`.
- The Work Queue manager pattern is `lobster_${USER}.*`, giving each user a distinct project namespace.
- The factory profile is `factory_skim.json`.
- The wrapper file is `./with_oasis_certs`.
- The wrapper is also passed as wrapper input with `--wrapper-input with_oasis_certs`.

The factory is an external, operator-managed service. The Lobster config does
not start or stop it. The operator is responsible for its credentials,
resource profile, project-pattern match, logs, lifetime, and shutdown. Record
factory state separately from Lobster manager state.

Run 2 Lobster submission command:

```bash
lobster process lobster_config.py
```

Run the Lobster command from the same `skimmer/` directory assumption unless a future source-backed runbook revision documents a different launch directory.

Non-goals and safety caveats:

- Do not run `work_queue_factory`, `lobster process`, Condor, XRootD, DBS, setup, installer, CMSSW, or worker commands unless a separate production procedure explicitly authorizes those capabilities.
- Do not infer production readiness from these documented commands; setup, services, credentials, storage, and worker behavior remain unvalidated here.
- Do not edit `factory_skim.json`, `with_oasis_certs`, wrapper code, setup scripts, sample cfgs, or nested `topeft` merely to match this runbook.
- Run 2 uses `--runos cc7-wq-7.11.1`; Run 3 uses `--runos al9-wq-7.11.1`.

## 25. Monitoring

Use the resolved project workdir when monitoring:

```bash
lobster status <workdir>
```

The operating record should identify `process.log`, `process.err`, `configure.log`, the resolved workdir and plotdir, factory logs, task exit codes, selected/failed/skipped counts, storage failures, worker resource exhaustion, and how campaign completion will be established. Submission acceptance is not campaign completion.

## 26. Retry, recovery, and cleanup guidance

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

## 27. Validation checklist before running

### static_preflight_checks

- Confirm the Run 2 repository path, intended commit, and tracked state.
- Inspect nested `topeft` status separately; preserve its branch and local changes.
- Confirm only intended configuration changes are present and review their diff.

### configuration_source_checks

- Record all user-knob values and calculate normalized `TYPE`, cfg kind, `TARGET`, `YEAR`, and `INPUT_MODE`.
- Confirm cfg classification is known/unambiguous and matches `TYPE` unless a reviewed override is intended.
- Confirm sandbox, wrapper, and `haddnano.py` paths.
- Review generated cut semantics and per-year module mapping from source.

### environment_checks

- Confirm setup/build shell separation from conda/mamba.
- Record production Lobster environment, CMSSW environment, Python resolution, credentials, and site runtime.
- Confirm `CMSSW_10_6_19_patch2` and required CMGTools/NanoAODTools payloads are complete.

### sample_metadata_checks

- Resolve the exact cfg and active JSON list after comments and `MATCH`.
- Require supported `year` for every selected JSON.
- For `files`, verify every file string; for `dbs`, verify dataset syntax; for `auto`, record every per-sample decision.
- Confirm the selected year filter does not remove all samples.

### path/output_checks

- Record `TAG`, master label, workdir, plotdir, logical output, local URL, and remote URL.
- Check for collisions, existing project state, quota, permissions, worker visibility, and endpoint policy.
- Confirm factory logs and monitoring destinations.

### production_run_checks

- Obtain explicit approval for factory, Lobster, Condor, network, storage, and worker execution.
- Record exact commands, host, working directory, external services, process IDs, and termination state.
- Validate selected workflows and outputs; do not equate static inspection, submission, or partial output with completion.

## 28. Troubleshooting

| symptom | likely source-visible cause | first check |
| --- | --- | --- |
| TYPE/cfg consistency error | Cfg basename tokens do not match `TYPE`, or are unknown/ambiguous. | Inspect `TYPE`, `CFG_NAME`, and basename tokens before paths are derived. |
| No samples remain | `MATCH` or `YEAR` excluded everything. | Review active cfg lines, regexes, JSON years, and printed counts. |
| Missing/invalid year | JSON lacks a supported string `year`. | Fix metadata in the independently owned source through an authorized change. |
| Input contract error | DBS path is malformed or `files` is empty/invalid. | Read the detailed metadata summary and selected `INPUT_MODE`. |
| Wrong module | JSON `year` is wrong. | Check per-sample year decision and `lepMVA_<year>` mapping. |
| Sandbox failure | Missing/incomplete CMSSW release or incompatible environment. | Verify exact sandbox path and payload contents; do not rerun setup blindly. |
| Stage-in failure | Input not local and fallback XRootD transfer fails. | Check resolved input URL, credentials, endpoint, and wrapper log. |
| Local stage-out failure | Worker cannot see local destination or its parent is missing. | Check `file://` path visibility and remote fallback. |
| Remote stage-out failure | Endpoint, credentials, permissions, or quota issue. | Inspect task transfer errors and destination policy. |
| Workdir collision | Same step/target/type/day reused. | Compare derived non-test workdir before launch. |
| Repeated failure threshold | Root cause persists or counters reached limits. | Use framework recovery guidance; do not delete workdir state. |
| Config works statically but campaign fails | Static inspection cannot prove services, environment, or payload runtime. | Perform a separately authorized bounded validation with exact records. |

## 29. Post-mortem EFT reweighting

`postmortem/globalEFTreWeighting.py` is separate Run 2-only material. The current notes describe a historical `CMSSW_10_6_26`/CMGTools 80X setup, which differs from the main skim sandbox. Exact upstream provenance, local delta, maintained example config, and production validation are incomplete. See `postmortem/README.md`; do not apply the procedure to Run 3 or infer readiness without a dedicated review.

## 30. Framework references

Consult the Lobster framework documentation for framework-level details:

- [Configuration](https://github.com/NDCMS/lobster/blob/lobster-python3/docs/config.rst)
- [Running jobs](https://github.com/NDCMS/lobster/blob/lobster-python3/docs/run.rst)
- [Monitoring](https://github.com/NDCMS/lobster/blob/lobster-python3/docs/monitor.rst)
- [Troubleshooting](https://github.com/NDCMS/lobster/blob/lobster-python3/docs/trouble.rst)
- [Failure threshold recovery](https://github.com/NDCMS/lobster/blob/lobster-python3/docs/FailureThresholdRecovery.md)
- [Operational guidance](https://github.com/NDCMS/lobster/blob/lobster-python3/docs/operational_guidance.md)
- [Storage file URLs](https://github.com/NDCMS/lobster/blob/lobster-python3/docs/StorageFileURLs.md)
- [Example factory configuration](https://github.com/NDCMS/lobster/blob/lobster-python3/examples/factory.json)

Some framework examples are historical or site-specific. This tutorial is the source of truth for current Run 2 repository-specific values; neither source supersedes an explicit current production authorization.

## 31. Run 2 validation scope

The current Run 2 documentation records static source behavior only. There is
no Run 2 counterpart to the one-sample worker-backed Run 3 evidence, no
recorded Run 2 config-based termination check, and no recorded Run 2 ROOT
content inspection. Use the stable-path guidance above when preparing a
bounded `TESTING` config, but do not treat the Run 3 results as Run 2 evidence.

## 32. Known caveats and intentionally unresolved items

- Fresh setup selects `run3_test_mmerged`; an existing nested checkout may differ, and setup does not reconcile it automatically.
- Setup/install bodies and production services remain unvalidated.
- Work Queue factory and Lobster submission commands are site-specific and require operator review before use.
- Run 2 does not include a source-controlled canonical factory JSON; verify the profile supplied at launch.
- Mixed DBS/files behavior is implemented but has no recorded Run 2 runtime validation.
- `skim_wrapper.py` is Python-2-style and relies on sandbox command resolution.
- The setup script's existing-check logic does not validate checkout branch or installation completeness.
- Post-mortem payload provenance and validation remain incomplete.
