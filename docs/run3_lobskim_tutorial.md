# Run 3 LobSkim tutorial

## Scope and assumptions

This tutorial documents the current Run 3 `LobsterSkimming` parent repository at `/users/apiccine/work/LobSkim/LobsterSkimming`. In this workspace, the local working branch is `run3-mva-anpicci`; the upstream branch this work is intended to land in is `run3_mva`.

The statements below are based on static inspection of the current source files. The setup entry point, installer scripts, `skimmer/lobster_config.py`, `skimmer/skim_wrapper.py`, Lobster commands, Work Queue factory, CMSSW setup, network access, storage access, and production campaign were not executed during this documentation update. Treat production behavior as unvalidated until an explicitly authorized validation or production round records exact commands and results.

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

Fresh setup is configured to obtain `https://github.com/TopEFT/topeft.git` at branch/tag `run3_test_mmerged`. That fresh-setup policy and an existing local checkout are different facts. At the time of this documentation update, the nested Run 3 checkout was on `run3_test_mmerged_anpicci` at `06f84d838f3aed6cc18fdd1a301b1fb5fb651624` and had a pre-existing cfg modification. Do not switch it as part of normal tutorial use. Reconciling `run3_test_mmerged` and `run3_test_mmerged_anpicci` is a separate maintenance decision requiring explicit authorization and a clean plan for local changes.

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

The canonical interpreter used for workspace diagnostics is `/users/apiccine/work/miniconda3/envs/LobSkim/bin/python`, through `/users/apiccine/work/LobSkim/codex-run.sh`. That diagnostic policy does not prove the worker payload interpreter or production control-shell environment.

## Setup overview

`setup.py` has a no-argument entry point and immediately calls `main()` when executed by a user. It finds the parent Git top level, then:

- skips `topeft` when that path already exists;
- otherwise calls `scripts/install_configs.sh` for the two sparse paths;
- skips the CMSSW installation when `CMSSW_14_0_6` already exists;
- otherwise calls `scripts/install_cmssw.sh` with the parent path, release, and SCRAM architecture.

Both shell installers use `set -euo pipefail` and positional-argument checks. This describes source-visible behavior only. The setup and installers were read statically; they were not run against GitHub, CVMFS, SCRAM, or storage in this documentation round.

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
- To run data, change `TYPE` to `"data"` and confirm the derived `ND_<YEAR>_data.cfg` exists.
- To select a different Run 3 period, change `YEAR` and review the resulting `TAG`, `CFG_NAME`, and output paths.
- To narrow a retry or debugging subset, edit `MATCH` with regexes that match cfg lines.
- To choose DBS discovery instead of explicit file lists, change `INPUT_MODE` to `"dbs"` and verify every selected JSON has a valid `path`.
- To run SR-like skim cuts, change `TARGET` to `"SR"`; leave `TARGET="CR"` for the looser control-region cut.

After edits, review the complete resolved configuration output produced when Lobster eventually evaluates the config in an authorized run; do not rely only on filenames.

## User knobs reference

| knob | default | allowed_or_expected_values | what_it_controls | when_to_edit | risk_or_caveat | source_file |
| --- | --- | --- | --- | --- | --- | --- |
| `TESTING` | `False` | Boolean | Switches labels and paths to timestamped `test/lobster_skimtest_*` locations. | For a deliberately bounded test campaign. | A test-labelled path does not make submission safe or lightweight. | `skimmer/lobster_config.py` |
| `INPUT_MODE` | `"files"` | Current source accepts `dbs` or `files` after lower-case normalization; intended design also includes `auto` after a source update. | Dataset object type and storage profile. | Use `files` for explicit JSON file lists; use `dbs` for dataset discovery from `path`; use `auto` only after it is implemented and validated. | `auto` input mode is not yet implemented in the current Run 3 source. Support is intended and should be added in a follow-up source update before relying on automatic per-sample input-mode selection. | same |
| `PROTOCOL_LOCAL` | `"file://"` | Storage protocol understood by Lobster; current design expects `file://`. | Local input/output prefix. | Only for an intentional storage-profile change. | `file://` requires worker-visible paths and existing destination parents. | same |
| `PROTOCOL_REMOTE` | `"root://"` | Current design expects `root://`. | Remote XRootD prefix. | Only when changing remote transport. | Endpoint reachability was not validated. | same |
| `TARGET` | `"CR"` | `SR`, `CR` after upper-case normalization | Changes the skim cut and path/label segment. | Choose signal-region or control-region skim semantics. | `SR` adds the opposite-sign two-lepton veto; `CR` does not. | same |
| `YEAR` | `"2023BPix"` | Source does not validate explicitly; practical current cfg periods are `2022`, `2022EE`, `2023`, `2023BPix`. | `TAG`, `CFG_NAME`, paths, and module selection through sample names. | Choose the Run 3 run period. | Invalid values fail later when the derived cfg is absent or samples do not match expected names. | same |
| `STEP` | `"skimmed"` | Nonempty campaign label chosen by operator. | Label/path segment. | When naming a distinct processing step. | Changing it creates different paths; it does not change physics behavior. | same |
| `TYPE` | `"background"` | Expected `background`, `signal`, or `data` because of cfg-name derivation. | `TAG` and `CFG_NAME`. | Match the campaign kind. | Source does not validate `TYPE`; unexpected values derive cfg names that may not exist. | same |
| `TAG` | `f"{TYPE}/NAOD_ULv9_lepMVA-run3/{YEAR}"` | Derived string unless deliberately edited. | Campaign/output namespace. | Rarely; usually edit `TYPE` or `YEAR`. | Manual edits can make output identity misleading. | same |
| `CFG_NAME` | `ND_2023BPix_background_samples.cfg` for defaults | `ND_<YEAR>_<TYPE>_samples.cfg` except data uses `ND_<YEAR>_data.cfg`; must exist under `topeft/input_samples/cfgs`. | Selects the sample list. | For custom cfg selection or when source derivation is insufficient. | Source does not separately check consistency between `TYPE` and a manual cfg override. | same |
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
| `TAG` | `<TYPE>/NAOD_ULv9_lepMVA-run3/<YEAR>` | `background/NAOD_ULv9_lepMVA-run3/2023BPix` | Campaign/output namespace. | `TYPE`, `YEAR`, or manual edit. | It is derived before sample validation; misleading values can route outputs incorrectly. |
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
CFG_NAME = f"ND_{YEAR}_{TYPE}_samples.cfg" if TYPE != "data" else f"ND_{YEAR}_data.cfg"
```

With the current defaults, `TYPE="background"` and `YEAR="2023BPix"` derive:

```text
TAG = background/NAOD_ULv9_lepMVA-run3/2023BPix
CFG_NAME = ND_2023BPix_background_samples.cfg
```

The source validates `INPUT_MODE` and `TARGET`, but it does not explicitly validate `TYPE`, `YEAR`, or `CFG_NAME`. Current available `ND_<year>_*` cfgs show practical Run 3 periods `2022`, `2022EE`, `2023`, and `2023BPix`, with background and signal cfgs following the derived pattern and data using `<year>_data.cfg` names in the observed checkout rather than `ND_<year>_data.cfg` for some periods. If `TYPE="data"` is intended, confirm the derived cfg actually exists before launching; if not, a separate source/config maintenance decision is needed.

`select_module_name` uses sample-name string matching. Samples containing `2022` or `2023` use module `lepMVA`; older UL names map to Run 2-style modules. This Run 3 config does not filter samples by JSON `year`; the selected cfg and `MATCH` determine the active sample list.

## TARGET and skim-cut behavior

Only `SR` and `CR` are valid. Both require at least two selected electrons, muons, and taus in total. Electron, muon, and tau object requirements are encoded directly in `build_skim_cut`.

`SR` additionally vetoes events with exactly two selected leptons whose summed charge is zero. `CR` uses only the at-least-two requirement. `TYPE` does not affect this cut. The generated expression has whitespace removed and is checked for balanced parentheses and accidental literal outer quotes.

## INPUT_MODE and sample input metadata

| mode | required metadata | dataset object | storage behavior |
| --- | --- | --- | --- |
| `files` | Nonempty `files` list containing usable file strings. | Lobster `Dataset(files=..., files_per_task=1, patterns=["*.root"])` | `storage_files`, with local and remote input prefixes and local and remote output prefixes. |
| `dbs` | `path` string suitable for `cmssw.Dataset(dataset=...)`. | `cmssw.Dataset(dataset=jsn["path"], lumis_per_task=1, file_based=True)` | `storage_dbs`, remote output only, and `AdvancedOptions.xrootd_servers=[SRC_REMOTE]`. |

The current Run 3 source accepts only `INPUT_MODE="files"` and `INPUT_MODE="dbs"`. `auto` input mode is not yet implemented in the current Run 3 source. Support is intended and should be added in a follow-up source update before relying on automatic per-sample input-mode selection. Until then, `/store/...`, `root://...`, and `file://...` file entries are consumed only by `files` mode. DBS mode uses the JSON `path` field and relies on Lobster/CMSSW dataset handling.

Each selected JSON should also carry analysis metadata such as `year`, `isData`, `xsec`, event counts, and weight sums. The current Run 3 config consumes `files` and `path` for input construction; it does not enforce or filter on JSON `year`.

## Storage and XRootD behavior

`storage_dbs` has only the remote output URL and leaves input streaming enabled. `storage_files` tries local then remote input prefixes and local then remote output prefixes, also with streaming enabled.

Selection is global:

- `INPUT_MODE="files"` uses `storage_files`.
- `INPUT_MODE="dbs"` uses `storage_dbs`.
- Intended `INPUT_MODE="auto"` support is a follow-up source update; do not rely on mixed per-sample selection until it is implemented and validated.

For DBS mode, `AdvancedOptions.xrootd_servers` is set to `["cms-xrd-global.cern.ch"]`. For files mode, the wrapper receives explicit input file names and may stage in missing local basenames with `xrdcp -f`.

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
| `--module` | config helper and wrapper parser | `lepMVA` for sample names containing `2022` or `2023`; older UL names map to year-specific modules. | Normally no; fix sample naming only through authorized metadata changes. | Wrong sample naming can select the wrong module. |
| `--out-dir` | config helper and wrapper parser | `.` | Normally no. | Wrapper and expected output assume task working directory. |
| `--nevents` | wrapper parser only | Optional wrapper flag exists, but current config does not pass it. | Not through current user knobs. | Adding it requires workflow-command review. |
| `@inputfiles` | Lobster substitution | Expanded to task inputs. | No. | Expansion semantics belong to Lobster. |
| `output.root` | wrapper/workflow | Final required task output. | No. | Missing output causes task failure/stage-out failure. |
| workflow merge | config | `haddnano.py @outputfiles @inputfiles`. | Only in a separately reviewed workflow change. | Distinct from wrapper's per-task merge. |

The config also creates a `shlex.join` display command for logs. That display form is not passed to Lobster because parser handling of shell quoting has not been established.

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

The wrapper does not run `haddnano.py` through `python`; it calls `haddnano.py` directly. It does not validate that the expected `_Skim.root` files exist before merge. It was not executed in this documentation round.

## Work Queue and Lobster runbook

These commands are user-supplied and source-visible documentation, not commands executed during this documentation update. They still require an explicitly authorized production round before use.

Run the factory and Lobster commands from the Run 3 `skimmer/` directory, where `lobster_config.py`, `factory_skim.json`, and `with_oasis_certs` are expected to be visible as working files. This records the working-directory assumption for the user-supplied `lobster process lobster_config.py` command; do not rewrite it as `lobster process skimmer/lobster_config.py` without a separate source-backed change.

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

Run 3 Lobster submission command:

```bash
lobster process lobster_config.py
```

Run the Lobster command from the same `skimmer/` directory assumption unless a future source-backed runbook revision documents a different launch directory.

Non-goals and safety caveats:

- Do not run `work_queue_factory`, `lobster process`, Condor, XRootD, DBS, setup, installer, CMSSW, or worker commands unless a separate production authorization explicitly permits those capabilities.
- Do not infer production readiness from these documented commands; setup, services, credentials, storage, and worker behavior remain unvalidated here.
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

## Monitoring placeholder

The exact status command remains unresolved. Use the user-approved status command only after it is supplied:

```text
<USER_SUPPLIED_RUN3_LOBSTER_STATUS_COMMAND>
```

The operating record should identify `process.log`, `process.err`, `configure.log`, the resolved workdir and plotdir, factory logs, task exit codes, selected/failed/skipped counts, storage failures, worker resource exhaustion, and how campaign completion will be established. Submission acceptance is not campaign completion.

## Retry, recovery, and cleanup guidance

First diagnose the source of failures. Preserve the workdir database and its config. Do not reset, delete, or clean project state as routine recovery.

Use only a reviewed command in place of:

```text
<USER_SUPPLIED_RUN3_LOBSTER_RECOVERY_OR_CLEANUP_COMMAND>
```

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

### production_execution_checks

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
- `/users/apiccine/work/LobSkim/lobster/docs/StorageFileURLs.md`
- `/users/apiccine/work/LobSkim/lobster/examples/factory.json`

Some framework examples are historical or site-specific. This tutorial is the source of truth for current Run 3 branch-specific values; neither source supersedes an explicit current production authorization.

## Known caveats and intentionally unresolved items

- Fresh setup selects `run3_test_mmerged`, while the observed nested checkout is `run3_test_mmerged_anpicci`; no reconciliation is prescribed here.
- Setup/install bodies and production services remain unvalidated.
- Exact Work Queue factory and Lobster submission commands are now documented from user-supplied text; status, recovery, cleanup, and shutdown commands remain user-supplied placeholders.
- Data cfg naming may need source/config review: current `CFG_NAME` derivation expects `ND_<YEAR>_data.cfg`, while observed Run 3 cfgs include names such as `2023BPix_data.cfg`.
- `TYPE` and `YEAR` are not explicitly validated by current Run 3 source.
- `auto` input mode is not yet implemented in the current Run 3 source. Support is intended and should be added in a follow-up source update before relying on automatic per-sample input-mode selection.
- Mixed DBS/files behavior is not implemented in current Run 3 source.
- `skim_wrapper.py` requires runtime tools such as `xrdcp`, `nano_postproc.py`, `CMGTools.NanoProc.tools.nanoAOD.lepMVA_run3`, and `haddnano.py`; none were executed here.
- The setup script's existing-check logic does not validate checkout branch or installation completeness.
