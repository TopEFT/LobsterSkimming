# Run 3 LobSkim / LobsterSkimming

This repository defines Run 3 NanoAOD skimming workflows for the Lobster framework.

Read [the Run 3 LobSkim tutorial](docs/run3_lobskim_tutorial.md) before changing or launching a campaign. It documents the current source-visible setup, configuration fields, sample metadata contracts, wrapper behavior, storage paths, and validation checklist.

The tutorial's [lepMVA module wiring and output branches](docs/run3_lobskim_tutorial.md#lepmva-module-wiring-and-output-branches) section traces module selection into `nano_postproc.py`, explains how removal from the processing chain must be handled, and identifies the source of the output branch names.

Important boundaries:

- `topeft/` is an independent nested Git repository, not a disposable untracked directory of this parent repository. Do not clean, reset, delete, or absorb it from the parent.
- Fresh setup selects the `TopEFT/topeft` branch/tag `run3_test_mmerged`. An existing checkout may differ; inspect and preserve it rather than replacing it automatically.
- Exact Work Queue factory and Lobster production commands are documented in [docs/run3_lobskim_tutorial.md](docs/run3_lobskim_tutorial.md). Do not infer alternate commands from examples or old README text.
- Setup and installer behavior described in the tutorial is based on static source inspection. Setup/install execution and production validation remain unrecorded.
- Older references to `lobster_config_T3.py` are stale. The current Run 3 configuration file is `skimmer/lobster_config.py`.
