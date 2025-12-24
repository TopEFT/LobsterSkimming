import datetime
import os
import sys
import subprocess

from lobster import cmssw
from lobster.core import AdvancedOptions, Category, Config, Dataset, StorageConfiguration, Workflow

sys.path.append(os.path.split(__file__)[0])
from tools.utils import read_cfg


# =============================================================================
# USER KNOBS (edit here)
# =============================================================================
TESTING = True

# Choose one:
#   "dbs"   -> cmssw.Dataset(dataset=..., file_based=True) + advanced.xrootd_servers
#   "files" -> Dataset(files=...) + StorageConfiguration(input=[root://...//])
INPUT_MODE = "dbs"

STEP = "skims"
TAG = "data/NAOD_ULv12_lepMVA-run3/2022/Full"  # not used in TESTING mode
CFG_NAME = "ND_2022_background_samples.cfg"

# Only process json files that match these regexs (empty list matches everything)
MATCH = [r".*ZG_MLL.*600.*\.json"]

SKIM_CUT = (
    " nMuon+nElectron+nTau >=2 && "
    "Sum$( Muon_looseId && Muon_miniPFRelIso_all<0.4 && Muon_sip3d<8 ) + "
    "Sum$( Electron_miniPFRelIso_all<0.4 && Electron_sip3d<8 ) + "
    "Sum$( Tau_idDeepTau2018v2p5VSe>1 && Tau_idDeepTau2018v2p5VSmu>0 && Tau_idDeepTau2018v2p5VSjet>1 ) >=2 "
)

WRAPPER = "skim_wrapper.py"  # keep T3 wrapper behavior, single file

# XRootD endpoints (T3 defaults)
XRD_SRC = "cms-xrd-global.cern.ch"          # reads (DBS mode redirector, or files mode prefix)
XRD_DST = "skynet013.crc.nd.edu:1094"       # writes

# Workdir base (your choice)
WORKDIR_BASE = "/tmpscratch/users/$USER"


# =============================================================================
# VALIDATION
# =============================================================================
INPUT_MODE = INPUT_MODE.strip().lower()
if INPUT_MODE not in ("dbs", "files"):
    raise ValueError(f"INPUT_MODE must be 'dbs' or 'files', got: {INPUT_MODE!r}")


# =============================================================================
# DERIVED PATHS / LABELS
# =============================================================================
TSTAMP1 = datetime.datetime.now().strftime("%Y%m%d_%H%M")
startingday = datetime.datetime.now().strftime("%y%m%d")
ver = f"v{startingday}"

top_dir = subprocess.check_output(["git", "rev-parse", "--show-toplevel"]).decode("utf-8").strip()

# Portable sandbox
sandbox_location = os.path.join(top_dir, "CMSSW_14_0_6")
print(f"Sandbox location: {sandbox_location}")

cfg_fpath = os.path.join(top_dir, "topeft", "input_samples", "cfgs", CFG_NAME)
print(f"Where is your cfg?\t {cfg_fpath}")

mode_stamp = INPUT_MODE
master_label = f"{STEP}_{mode_stamp}_lobPY3_{TSTAMP1}"
workdir_path = f"{WORKDIR_BASE}/{STEP}/{TAG}/{ver}"
plotdir_path = f"~/www/lobster/{STEP}/{TAG}/{ver}"
output_path  = f"/store/user/$USER/{STEP}/{TAG}/{ver}"

if TESTING:
    master_label = f"{STEP}_{mode_stamp}_testlobPY3_{TSTAMP1}"
    workdir_path = f"{WORKDIR_BASE}/{STEP}/test/lobster_skimtest_{TSTAMP1}"
    plotdir_path = f"~/www/lobster/{STEP}/test/lobster_skimtest_{TSTAMP1}"
    output_path  = f"/store/user/$USER/{STEP}/test/lobster_skimtest_{TSTAMP1}"

print(f"INPUT_MODE = {INPUT_MODE}")
print(f"XRD_SRC    = {XRD_SRC}")
print(f"XRD_DST    = {XRD_DST}")


# =============================================================================
# STORAGE (mode-dependent)
# =============================================================================
storage_dbs = StorageConfiguration(
    output=[
        f"file:///cms/cephfs/data/{output_path}",
        f"root://{XRD_DST}/{output_path}",
    ],
    disable_input_streaming=False,
)

storage_files = StorageConfiguration(
    input=[
        f"root://{XRD_SRC}//",  # note the double slash after host
    ],
    output=[
        f"file:///cms/cephfs/data/{output_path}",
        f"root://{XRD_DST}/{output_path}",
    ],
    disable_input_streaming=False,
)

storage = storage_dbs if INPUT_MODE == "dbs" else storage_files


# =============================================================================
# WORKFLOWS
# =============================================================================
cfg = read_cfg(cfg_fpath, match=MATCH)
print("cfg jsons:", list(cfg["jsons"].keys()))

cat = Category(name="processing", cores=1, memory=1500, disk=4500)

skim_cut = SKIM_CUT.replace(" ", "")

wf = []
for sample in sorted(cfg["jsons"]):
    jsn = cfg["jsons"][sample]
    print(f"Processing sample: {sample}")
    for fn in jsn["files"]:
        print(f"  {fn}")

    # keep only this (your decision h)
    files = [x for x in jsn["files"]]

    # module selection (bug-fixed)
    if "HIPM_UL2016" in sample:
        module_name = "lepMVA_2016_preVFP"
    elif "UL2017" in sample:
        module_name = "lepMVA_2017"
    elif "UL2018" in sample:
        module_name = "lepMVA_2018"
    elif ("2022" in sample) or ("2023" in sample):
        module_name = "lepMVA"
    else:
        module_name = "lepMVA_2016"

    # dataset selection (mode-dependent)
    if INPUT_MODE == "dbs":
        ds = cmssw.Dataset(
            dataset=jsn["path"],
            lumis_per_task=1,   # since file_based=True, this is effectively files_per_task
            file_based=True,
        )
    else:
        ds = Dataset(
            files=files,
            files_per_task=1,
            patterns=["*.root"],
        )

    cmd = ["python3", WRAPPER]
    cmd += ["--cut", skim_cut]
    cmd += ["--module", module_name]
    cmd += ["--out-dir", "."]
    cmd += ["@inputfiles"]

    print("\nRemote command:\n", " ".join(cmd), "\n")

    skim_wf = Workflow(
        label=sample.replace("-", "_"),
        sandbox=cmssw.Sandbox(release=sandbox_location),
        dataset=ds,
        category=cat,
        extra_inputs=[WRAPPER],
        outputs=["output.root"],
        command=" ".join(cmd),
        merge_command="haddnano.py @outputfiles @inputfiles",
        merge_size="537M",
        globaltag=False,
        cleanup_input=False,
    )
    wf.append(skim_wf)


# =============================================================================
# ADVANCED OPTIONS (mode-dependent)
# =============================================================================
adv_kwargs = dict(
    bad_exit_codes=[127, 160],
    log_level=1,
    payload=10,
    osg_version="3.6",
    threshold_for_failure=1,
    threshold_for_skipping=1,
)

if INPUT_MODE == "dbs":
    adv_kwargs["xrootd_servers"] = [XRD_SRC]

config = Config(
    label=master_label,
    workdir=workdir_path,
    plotdir=plotdir_path,
    storage=storage,
    workflows=wf,
    advanced=AdvancedOptions(**adv_kwargs),
)