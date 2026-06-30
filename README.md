# Run 3 LobSkim / LobsterSkimming

This repository defines the Run 3 NanoAOD skimming workflows used with the local Lobster framework. The current integration branch for this workspace is `run3-mva-anpicci`.

Read [the Run 3 LobSkim tutorial](docs/run3_lobskim_tutorial.md) before changing or launching a campaign. It documents the current source-visible setup, configuration fields, sample metadata contracts, wrapper behavior, storage paths, and validation checklist.

Important boundaries:

- `topeft/` is an independent nested Git repository, not a disposable untracked directory of this parent repository. Do not clean, reset, delete, or absorb it from the parent.
- Fresh setup selects `TopEFT/topeft` branch/tag `run3_test_mmerged`. Existing local checkouts may be on `run3_test_mmerged_anpicci`; reconciling them is a separate authorized maintenance decision.
- Exact Work Queue factory and Lobster production commands are documented in [docs/run3_lobskim_tutorial.md](docs/run3_lobskim_tutorial.md). Do not infer alternate commands from examples or old README text.
- Setup and installer behavior described in the tutorial was established by static source inspection; setup/install execution and production validation were not run by this documentation round.
- Older README references such as `lobster_config_T3.py` are stale for this branch. The current Run 3 configuration file is `skimmer/lobster_config.py`.
