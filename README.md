# Run 2 LobSkim / LobsterSkimmingRun2

This repository defines the Run 2 NanoAOD skimming workflows used with the local Lobster framework. The current integration branch is `run2_el9_skims`.

Read [the Run 2 LobSkim tutorial](docs/run2_lobskim_tutorial.md) before changing or launching a campaign. It documents the current source-visible setup, configuration fields, sample metadata contracts, wrapper behavior, storage paths, and validation checklist.

Important boundaries:

- `topeft/` is an independent nested Git repository, not a disposable untracked directory of this parent repository. Do not clean, reset, delete, or absorb it from the parent.
- Fresh setup selects `TopEFT/topeft` branch/tag `run3_test_mmerged`. Existing local checkouts may be on `run3_test_mmerged_anpicci`; reconciling them is a separate authorized maintenance decision.
- Exact Work Queue factory and Lobster production commands are documented in [docs/run2_lobskim_tutorial.md](docs/run2_lobskim_tutorial.md). Local untracked factory JSON files are not canonical documentation.
- Setup and installer behavior described in the tutorial was established by static source inspection; production/network execution was not validated by the documentation round.

Run 2-only post-mortem EFT reweighting notes are in [postmortem/README.md](postmortem/README.md).
