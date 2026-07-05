# Run 2 LobSkim / LobsterSkimmingRun2

This repository defines Run 2 NanoAOD skimming workflows for the Lobster framework.

Read [the Run 2 LobSkim tutorial](docs/run2_lobskim_tutorial.md) before changing or launching a campaign. It documents the current source-visible setup, configuration fields, sample metadata contracts, wrapper behavior, storage paths, and validation checklist.

The tutorial's [lepMVA module wiring and output branches](docs/run2_lobskim_tutorial.md#lepmva-module-wiring-and-output-branches) section traces module selection into `nano_postproc.py`, explains how removal from the processing chain must be handled, and identifies the source of the output branch names.

Important boundaries:

- `topeft/` is an independent nested Git repository, not a disposable untracked directory of this parent repository. Do not clean, reset, delete, or absorb it from the parent.
- Fresh setup selects the `TopEFT/topeft` branch/tag `run3_test_mmerged`. An existing checkout may differ; inspect and preserve it rather than replacing it automatically.
- Exact Work Queue factory and Lobster production commands are documented in [docs/run2_lobskim_tutorial.md](docs/run2_lobskim_tutorial.md).
- Setup and installer behavior described in the tutorial is based on static source inspection. Production and network behavior remain unvalidated.

Run 2-only post-mortem EFT reweighting notes are in [postmortem/README.md](postmortem/README.md).
