# LobsterSkimming
CMS skimming code that's been lobsterized for use on ND's computing resources


### Setup dependencies
Run the python setup script
```
python setup.py
```
This script will set up a basic CMSSW release and check out the necessary [NanoAOD skimming](https://github.com/cms-nanoAOD/nanoAOD-tools) package. The script also automatically gets the cfg and json directories from the [TopEFT/topeft Run3 branch](https://github.com/TopEFT/topeft/tree/run3_test_mmerged) repo.

The repo is currently setup to automatically generate and use a `CMSSW_14_0_6` release, but can be changed to use a different release. The simplest method would be to modify the [setup.py](https://github.com/Andrew42/LobsterSkimming/blob/master/setup.py) script, re-run it, and then also update the [lobster_config_T3.py](https://github.com/anpicci/LobsterSkimming/blob/run3-mva/skimmer/skim_wrapper_T3.py) config by changing the `sandbox_location` variable.

The CMSSW release that is used to construct the environment that actually runs the skimming job does not necessarily have to match the CMSSW release that is used to provide the dependencies for the lobster program itself. So changing what CMSSW release is used for doing the actual skimming should not have to result in a change to the CMSSW release used for your lobster environment.

### Running the code
Assuming that you installed lobster in a dedicated conda/mamba environment (see [this lobster project branch](https://github.com/anpicci/lobster/tree/lobster-python3-run3))
- Make sure to have a working CMSSW release in your local `LobsterSkimming` area:
```
cd CMSSW_X_Y_Z/src
cmsenv
```
- In a new clean shell, activate your lobster conda/mamba environment
```
conda activate lobster_env
## mamba activate lobster_env
```
- Run the code
```
lobster process lobster_config_T3.py
```

### Configuring the lobster job
The main configuration option in `lobster_config_T3.py` is the choice of the `cfg_name` variable. This should correspond to one of the cfg names from the topeft cfg files and will determine which dataset samples should be run over. The `step`, `tag`, and `ver` variables are used to define the directory names of the lobster job. If `testing=True`, then only the `step` variable will be used and the output directories will instead be changed to something of the form `$USER/{step}/test/lobster_test_{tstamp}`, where `{tstamp}` will correspond to a datetime timestamp of when the lobster job was started.
