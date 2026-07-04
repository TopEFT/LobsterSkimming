# Post-mortem EFT reweighting

## Scope and status

This material is Run 2-specific. It has not been validated for Run 3. Start with the [main Run 2 LobSkim tutorial](../docs/run2_lobskim_tutorial.md) for repository, environment, campaign, and nested-`topeft` policy.

`globalEFTreWeighting.py` is a local post-mortem payload derived from CMGTools-era functionality. Its exact upstream provenance and local delta are not yet fully documented, and it has no recorded runtime validation. Review the script and establish provenance before production use.

## Historical setup outline

The existing notes use a separate `CMSSW_10_6_26` area and CMGTools 80X content. Perform CMSSW setup outside conda/mamba, then activate the intended runtime environment only after the CMSSW area has been prepared. The historical sequence involves `cmsrel`, `cmsenv`, `git cms-init`, the `heppy_80X` CMSSW branch, a CMGTools 80X checkout, copying `globalEFTreWeighting.py` into CMGTools, and a SCRAM build.

These are historical requirements, not commands validated by this repository documentation. They also differ from the main skimming sandbox (`CMSSW_10_6_19_patch2`). Confirm release compatibility, upstream commits, local changes, and output validation before running.

## Running

The previously referenced `skimmer/lobster_config_p3_post-mortem.py` example is not present in the current repository and therefore is not maintained guidance. Before this procedure is treated as supported, provide a tracked, reviewed example configuration and record:

- the exact source commit and local delta for `globalEFTreWeighting.py`;
- the supported CMSSW and CMGTools commits;
- the input and output data contracts;
- a bounded validation result;
- user-approved Work Queue and Lobster commands.

Until that work is complete, do not infer production readiness from this README.
