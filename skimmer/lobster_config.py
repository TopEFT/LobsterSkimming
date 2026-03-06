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
TESTING = False

# Choose one:
#   "dbs"   -> cmssw.Dataset(dataset=..., file_based=True) + advanced.xrootd_servers
#   "files" -> Dataset(files=...) + StorageConfiguration(input=[root://...//])
INPUT_MODE = "files"

# Choose one protocol:
#   "root://"  # used for remote XRootD paths
#   "file://"  # used for local CephFS paths
PROTOCOL_LOCAL  = "file://"
PROTOCOL_REMOTE = "root://"

TARGET = "SR" # impact the inclusion of the two-lepton veto in the skim cut
YEAR = "2023BPix"
STEP = "skimmed"  # for labeling only, no functional impact
TYPE = "background"  # "data" "signal" for labeling only, no functional impact 
TAG = f"{TYPE}/NAOD_ULv12_lepMVA-run3/{YEAR}"  # not used in TESTING mode
CFG_NAME = f"ND_{YEAR}_{TYPE}_samples.cfg" if TYPE != "data" else f"{YEAR}_data.cfg"  # config file listing input samples (relative to top-level dir)

# Only process json files that match these regexs (empty list matches everything)
# MATCH = [r".*TTLL\_MLL-50.*\.json"]
MATCH = []

el_req = "Electron_pt>10 && Electron_miniPFRelIso_all<0.4 && Electron_sip3d<8 && (abs(Electron_eta + Electron_deltaEtaSC) < 1.4442 || abs(Electron_eta + Electron_deltaEtaSC) > 1.566))"
mu_req = "Muon_pt>10 && Muon_looseId && Muon_miniPFRelIso_all<0.4 && Muon_sip3d<8"
ta_req = "Tau_pt>20 && Tau_dxy<0.05 && Tau_dz<1 && Tau_idDeepTau2018v2p5VSe>1 && Tau_idDeepTau2018v2p5VSmu>3 && Tau_idDeepTau2018v2p5VSjet>2"

num_leps = f"(Sum$({el_req}) + Sum$({mu_req}) + Sum$({ta_req}))"
net_charge = f"(Sum$(Electron_charge*({el_req})) + Sum$(Muon_charge*({mu_req})) + Sum$(Tau_charge*({ta_req})))"

os2l_veto = f"(({num_leps} != 2) || ({net_charge} != 0))"
nlep_req = f"({num_leps} >= 2)"

SKIM_CUT = f"{nlep_req} && {os2l_veto}" if TARGET == "SR" else f"{nlep_req}"

WRAPPER = "skim_wrapper.py"  # keep T3 wrapper behavior, single file

# XRootD endpoints (T3 defaults)
## Different xrd src redirectors depending on where the inputs are stored
# SRC = "cmsxcache.crc.nd.edu"          # To read ND T3 files from outside of ND T3 (like the opportunistic resources)
# SRC = "cms-xrd-global.cern.ch"          # To read remote files (global redirector)
SRC_REMOTE = "cmsxrootd.crc.nd.edu"            # To read ND T3 files from ND T3 via XRootD
SRC_LOCAL = "/cms/cephfs/data" # if PROTOCOL == "file://" else SRC  # Local CephFS access
DST_REMOTE = "cmsxrootd.crc.nd.edu"            # To write to ND T3
DST_LOCAL = "/cms/cephfs/data" # if PROTOCOL == "file://" else DST  # Local CephFS access
# DST_LOCAL = "/project01/ndcms" # VAST local storage access (make sure to set PROTOCOL = "file://")

# if SRC.startswith("/cms/cephfs/data") and PROTOCOL != "file://":
#     raise ValueError("When using CephFS local path for SRC, PROTOCOL must be 'file://'")
# elif not SRC.startswith("/cms/cephfs/data") and PROTOCOL != "root://":
#     raise ValueError("When using remote XRootD path for SRC, PROTOCOL must be 'root://'")

# if DST.startswith("/cms/cephfs/data") and PROTOCOL != "file://":
#     raise ValueError("When using CephFS local path for DST, PROTOCOL must be 'file://'")
# elif not DST.startswith("/cms/cephfs/data") and PROTOCOL != "root://":
#     raise ValueError("When using remote XRootD path for DST, PROTOCOL must be 'root://'")

SRC_PREFIX_LOCAL = PROTOCOL_LOCAL + SRC_LOCAL + "//"
DST_PREFIX_LOCAL = PROTOCOL_LOCAL + DST_LOCAL + "//"

SRC_PREFIX_REMOTE = PROTOCOL_REMOTE + SRC_REMOTE + "//"
DST_PREFIX_REMOTE = PROTOCOL_REMOTE + DST_REMOTE + "//"

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

try:
    top_dir = subprocess.check_output(["git", "rev-parse", "--show-toplevel"]).decode("utf-8").strip()

    # Portable sandbox
    sandbox_location = os.path.join(top_dir, "CMSSW_14_0_6")
    print(f"Sandbox location: {sandbox_location}")

    cfg_fpath = os.path.join(top_dir, "topeft", "input_samples", "cfgs", CFG_NAME)
    print(f"Where is your cfg?\t {cfg_fpath}")
except Exception as e:
    print("Error while determining top-level directory or cfg file path. Please check your git repository and configuration. Error details:")
    print(e)
    raise Exception("Failed to set up configuration due to the above error.") # to comment out if you are modifying from the lobster workdir

mode_stamp = INPUT_MODE
master_label = f"{STEP}_{TARGET}_{mode_stamp}_lobPY3_{TSTAMP1}"
workdir_path = f"{WORKDIR_BASE}/{STEP}_{TARGET}/{TAG}/{ver}"
plotdir_path = f"~/www/lobster/{STEP}_{TARGET}/{TAG}/{ver}"
output_path  = f"/store/user/$USER/{STEP}_{TARGET}/{TAG}/{ver}"

if TESTING:
    master_label = f"{STEP}_{TARGET}_{mode_stamp}_testlobPY3_{TSTAMP1}"
    workdir_path = f"{WORKDIR_BASE}/{STEP}_{TARGET}/test/lobster_skimtest_{TSTAMP1}"
    plotdir_path = f"~/www/lobster/{STEP}_{TARGET}/test/lobster_skimtest_{TSTAMP1}"
    output_path  = f"/store/user/$USER/{STEP}_{TARGET}/test/lobster_skimtest_{TSTAMP1}"

print(f"INPUT_MODE = {INPUT_MODE}")
print(f"SRC_PREFIX_LOCAL = {SRC_PREFIX_LOCAL}")
print(f"SRC_PREFIX_REMOTE = {SRC_PREFIX_REMOTE}")
print(f"DST_PREFIX_LOCAL = {DST_PREFIX_LOCAL}")
print(f"DST_PREFIX_REMOTE = {DST_PREFIX_REMOTE}")
print(f"{DST_PREFIX_LOCAL}{output_path}", "\n")
print(f"{DST_PREFIX_REMOTE}{output_path}", "\n")


# =============================================================================
# STORAGE (mode-dependent)
# =============================================================================
storage_dbs = StorageConfiguration(
    output=[
        f"{DST_PREFIX_REMOTE}{output_path}",
    ],
    disable_input_streaming=False,
)

storage_files = StorageConfiguration(
    input=[
        f"{SRC_PREFIX_LOCAL}",
        f"{SRC_PREFIX_REMOTE}"
    ],
    output=[
        f"{DST_PREFIX_LOCAL}{output_path}",
        f"{DST_PREFIX_REMOTE}{output_path}",
    ],
    disable_input_streaming=False,
)

storage = storage_dbs if INPUT_MODE == "dbs" else storage_files


# =============================================================================
# WORKFLOWS
# =============================================================================
try:
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
except Exception as e:
    print("Error while setting up workflows. Please check your configuration and try again. Error details:")
    print(e)
    wf = []  # set to empty list to avoid further errors in lobster execution
    raise Exception("Failed to set up workflows due to the above error.") # to comment out if you are modifying from the lobster workdir

# =============================================================================
# ADVANCED OPTIONS (mode-dependent)
# =============================================================================
adv_kwargs = dict(
    bad_exit_codes=[127, 160],
    log_level=1,
    payload=10,
    osg_version="3.6",
    threshold_for_failure=10,
    threshold_for_skipping=10,
)

if INPUT_MODE == "dbs":
    adv_kwargs["xrootd_servers"] = [SRC_REMOTE]

config = Config(
    label=master_label,
    workdir=workdir_path,
    plotdir=plotdir_path,
    storage=storage,
    workflows=wf,
    advanced=AdvancedOptions(**adv_kwargs),
)