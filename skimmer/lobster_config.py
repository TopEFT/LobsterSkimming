import datetime
import os
import re
import shlex
import subprocess
import sys

from lobster import cmssw
from lobster.core import AdvancedOptions, Category, Config, Dataset, StorageConfiguration, Workflow

sys.path.append(os.path.split(__file__)[0])
from tools.utils import read_cfg


# =============================================================================
# USER KNOBS
# =============================================================================

TESTING = False

# Choose one:
#   "dbs"   -> cmssw.Dataset(dataset=..., file_based=True) + advanced.xrootd_servers
#   "files" -> Dataset(files=...) + StorageConfiguration(input=[root://...//])
INPUT_MODE = "files"

# Choose one protocol:
#   "root://"  # used for remote XRootD paths
#   "file://"  # used for local CephFS paths
PROTOCOL_LOCAL = "file://"
PROTOCOL_REMOTE = "root://"

TARGET = "CR"  # impacts the inclusion of the two-lepton veto in the skim cut
YEAR = "2023BPix"
STEP = "skimmed"
TYPE = "background"
TAG = f"{TYPE}/NAOD_ULv9_lepMVA-run3/{YEAR}"
CFG_NAME = f"ND_{YEAR}_{TYPE}_samples.cfg" if TYPE != "data" else f"ND_{YEAR}_data.cfg"

# Empty list matches everything.
# For first retry, strongly consider something like:
# MATCH = [r".*ZG.*PTG_200to400.*\.json"]
MATCH = []

WRAPPER = "skim_wrapper.py"

# XRootD endpoints
SRC_REMOTE = "cms-xrd-global.cern.ch"
SRC_LOCAL = "/cms/cephfs/data"

DST_REMOTE = "cmsxrootd.crc.nd.edu"
DST_LOCAL = "/cms/cephfs/data"

WORKDIR_BASE = "/tmpscratch/users/$USER"


# =============================================================================
# HELPERS
# =============================================================================

def remove_whitespace(expr):
    return re.sub(r"\s+", "", expr)


def assert_balanced_parentheses(expr, label):
    depth = 0

    for idx, char in enumerate(expr):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1

        if depth < 0:
            raise ValueError(
                f"{label} has an unmatched closing parenthesis at position {idx}: {expr}"
            )

    if depth != 0:
        raise ValueError(
            f"{label} has unbalanced parentheses with final depth {depth}: {expr}"
        )


def assert_no_literal_outer_quotes(expr, label):
    if len(expr) >= 2 and expr[0] == expr[-1] and expr[0] in ("'", '"'):
        raise ValueError(
            f"{label} starts and ends with literal quote characters. "
            f"Do not include shell quotes in the Lobster argument value: {expr!r}"
        )


def build_skim_cut(target):
    # Keep each per-object requirement wrapped in one explicit outer pair of parentheses.
    # This makes Sum$({req}) and charge-weighted Sum$ expressions unambiguous.
    el_req = (
        "(Electron_pt>10"
        " && Electron_miniPFRelIso_all<0.4"
        " && Electron_sip3d<8"
        " && (abs(Electron_eta + Electron_deltaEtaSC)<1.4442"
        " || abs(Electron_eta + Electron_deltaEtaSC)>1.566))"
    )

    mu_req = (
        "(Muon_pt>10"
        " && Muon_looseId"
        " && Muon_miniPFRelIso_all<0.4"
        " && Muon_sip3d<8)"
    )

    ta_req = (
        "(Tau_pt>20"
        " && abs(Tau_dxy)<0.05"
        " && abs(Tau_dz)<1"
        " && Tau_idDeepTau2018v2p5VSe>1"
        " && Tau_idDeepTau2018v2p5VSmu>3"
        " && Tau_idDeepTau2018v2p5VSjet>2)"
    )

    for label, expr in [
        ("el_req", el_req),
        ("mu_req", mu_req),
        ("ta_req", ta_req),
    ]:
        assert_balanced_parentheses(expr, label)

    num_leps = f"(Sum$({el_req}) + Sum$({mu_req}) + Sum$({ta_req}))"

    net_charge = (
        f"(Sum$(Electron_charge*{el_req})"
        f" + Sum$(Muon_charge*{mu_req})"
        f" + Sum$(Tau_charge*{ta_req}))"
    )

    os2l_veto = f"(({num_leps} != 2) || ({net_charge} != 0))"
    nlep_req = f"({num_leps} >= 2)"

    if target == "SR":
        skim_cut = f"{nlep_req} && {os2l_veto}"
    else:
        skim_cut = nlep_req

    skim_cut = remove_whitespace(skim_cut)

    for label, expr in [
        ("num_leps", remove_whitespace(num_leps)),
        ("net_charge", remove_whitespace(net_charge)),
        ("os2l_veto", remove_whitespace(os2l_veto)),
        ("nlep_req", remove_whitespace(nlep_req)),
        ("skim_cut", skim_cut),
    ]:
        assert_balanced_parentheses(expr, label)
        assert_no_literal_outer_quotes(expr, label)

    return skim_cut


def select_module_name(sample):
    if "HIPM_UL2016" in sample:
        return "lepMVA_2016_preVFP"
    if "UL2017" in sample:
        return "lepMVA_2017"
    if "UL2018" in sample:
        return "lepMVA_2018"
    if "2022" in sample or "2023" in sample:
        return "lepMVA"

    return "lepMVA_2016"


def build_payload_command(wrapper, skim_cut, module_name, out_dir):
    cmd = ["python3", wrapper]
    cmd += ["--cut", skim_cut]
    cmd += ["--module", module_name]
    cmd += ["--out-dir", out_dir]
    cmd += ["@inputfiles"]

    # Important:
    # - command_string is the actual string passed to Lobster.
    # - display_command is only for human logs/copy-paste.
    # Do not use shlex.join as the actual Lobster command unless Lobster's parser is
    # confirmed to handle shell quoting exactly as expected.
    command_string = " ".join(cmd)
    display_command = shlex.join(cmd)

    return command_string, display_command


def sample_names_preview(sample_names, limit=30):
    names = [str(name) for name in sample_names]
    preview = names[:limit]
    truncated = max(len(names) - len(preview), 0)
    return preview, truncated


def count_values(values):
    counts = {}
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    return dict(sorted(counts.items()))


def format_preview_value(values, limit=30):
    if values is None:
        items = []
    elif isinstance(values, (list, tuple)):
        items = list(values)
    else:
        items = [values]

    preview, truncated = sample_names_preview(items, limit=limit)
    preview_text = ", ".join(preview) if preview else "none"
    if truncated:
        preview_text = f"{preview_text} (+{truncated} more)"
    return preview_text


def format_resolved_config_lines(
    header,
    fields,
    sample_names,
    payload_command_actual,
    payload_command_display,
    preview_limit=30,
):
    sample_preview, sample_truncated = sample_names_preview(
        sample_names,
        limit=preview_limit,
    )
    lines = [header]
    for key, value in fields:
        lines.append(f"  {key}: {value}")
    lines.append(f"  selected_sample_count: {len(sample_names)}")
    lines.append(
        "  selected_sample_names_preview: "
        f"{', '.join(sample_preview) if sample_preview else 'none'}"
    )
    lines.append(f"  selected_sample_names_truncated: {sample_truncated}")
    lines.append(
        "  payload_command_actual: "
        f"{format_preview_value(payload_command_actual, limit=5)}"
    )
    lines.append(
        "  payload_command_display: "
        f"{format_preview_value(payload_command_display, limit=5)}"
    )
    return lines


def print_resolved_config(
    header,
    fields,
    sample_names,
    payload_command_actual,
    payload_command_display,
    preview_limit=30,
):
    for line in format_resolved_config_lines(
        header=header,
        fields=fields,
        sample_names=sample_names,
        payload_command_actual=payload_command_actual,
        payload_command_display=payload_command_display,
        preview_limit=preview_limit,
    ):
        print(line)


# =============================================================================
# VALIDATION
# =============================================================================

INPUT_MODE = INPUT_MODE.strip().lower()
if INPUT_MODE not in ("dbs", "files"):
    raise ValueError(f"INPUT_MODE must be 'dbs' or 'files', got: {INPUT_MODE!r}")

TARGET = TARGET.strip().upper()
if TARGET not in ("SR", "CR"):
    raise ValueError(f"TARGET must be 'SR' or 'CR', got: {TARGET!r}")


# =============================================================================
# DERIVED PATHS / LABELS
# =============================================================================

TSTAMP1 = datetime.datetime.now().strftime("%Y%m%d_%H%M")
startingday = datetime.datetime.now().strftime("%y%m%d")
ver = f"v{startingday}"

try:
    top_dir = subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"]
    ).decode("utf-8").strip()

    sandbox_location = os.path.join(top_dir, "CMSSW_14_0_6")
    print(f"Sandbox location: {sandbox_location}")

    cfg_fpath = os.path.join(top_dir, "topeft", "input_samples", "cfgs", CFG_NAME)
    print(f"Where is your cfg?\t {cfg_fpath}")

except Exception as err:
    print(
        "Error while determining top-level directory or cfg file path. "
        "Please check your git repository and configuration. Error details:"
    )
    print(err)
    raise RuntimeError("Failed to set up configuration due to the above error.") from err


mode_stamp = INPUT_MODE
master_label = f"{STEP}_{TARGET}_{mode_stamp}_lobPY3_{TSTAMP1}"
workdir_path = f"{WORKDIR_BASE}/{STEP}_{TARGET}/{TAG}/{ver}"
plotdir_path = f"~/www/lobster/{STEP}_{TARGET}/{TAG}/{ver}"
output_path = f"/store/user/$USER/{STEP}_{TARGET}/{TAG}/{ver}"

if TESTING:
    master_label = f"{STEP}_{TARGET}_{mode_stamp}_testlobPY3_{TSTAMP1}"
    workdir_path = f"{WORKDIR_BASE}/{STEP}_{TARGET}/test/lobster_skimtest_{TSTAMP1}"
    plotdir_path = f"~/www/lobster/{STEP}_{TARGET}/test/lobster_skimtest_{TSTAMP1}"
    output_path = f"/store/user/$USER/{STEP}_{TARGET}/test/lobster_skimtest_{TSTAMP1}"


# =============================================================================
# STORAGE
# =============================================================================

SRC_PREFIX_LOCAL = PROTOCOL_LOCAL + SRC_LOCAL + "//"
DST_PREFIX_LOCAL = PROTOCOL_LOCAL + DST_LOCAL + "//"

SRC_PREFIX_REMOTE = PROTOCOL_REMOTE + SRC_REMOTE + "//"
DST_PREFIX_REMOTE = PROTOCOL_REMOTE + DST_REMOTE + "//"

storage_dbs = StorageConfiguration(
    output=[
        f"{DST_PREFIX_REMOTE}{output_path}",
    ],
    disable_input_streaming=False,
)

storage_files = StorageConfiguration(
    input=[
        SRC_PREFIX_LOCAL,
        SRC_PREFIX_REMOTE,
    ],
    output=[
        f"{DST_PREFIX_LOCAL}{output_path}",
        f"{DST_PREFIX_REMOTE}{output_path}",
    ],
    disable_input_streaming=False,
)

storage = storage_dbs if INPUT_MODE == "dbs" else storage_files
MERGE_COMMAND = "haddnano.py @outputfiles @inputfiles"


# =============================================================================
# SKIM CUT
# =============================================================================

skim_cut = build_skim_cut(TARGET)

print(f"INPUT_MODE = {INPUT_MODE}")
print(f"TARGET = {TARGET}")
print(f"SRC_PREFIX_LOCAL = {SRC_PREFIX_LOCAL}")
print(f"SRC_PREFIX_REMOTE = {SRC_PREFIX_REMOTE}")
print(f"DST_PREFIX_LOCAL = {DST_PREFIX_LOCAL}")
print(f"DST_PREFIX_REMOTE = {DST_PREFIX_REMOTE}")
print(f"Local output path:  {DST_PREFIX_LOCAL}{output_path}")
print(f"Remote output path: {DST_PREFIX_REMOTE}{output_path}")
print(f"skim_cut repr: {skim_cut!r}")
print(f"skim_cut: {skim_cut}")

assert_no_literal_outer_quotes(skim_cut, "skim_cut")
assert_balanced_parentheses(skim_cut, "skim_cut")


# =============================================================================
# WORKFLOWS
# =============================================================================

try:
    cfg = read_cfg(cfg_fpath, match=MATCH)
    print(f"cfg_json_count: {len(cfg['jsons'])}")
    print(
        "cfg_json_names_preview: "
        f"{format_preview_value(list(cfg['jsons'].keys()), limit=30)}"
    )

    category_name = "processing"
    category_cores = 1
    category_memory = 1500
    category_disk = 4500
    cat = Category(
        name=category_name,
        cores=category_cores,
        memory=category_memory,
        disk=category_disk,
    )

    workflows = []
    payload_command_actual_values = []
    payload_command_display_values = []

    for sample in sorted(cfg["jsons"]):
        jsn = cfg["jsons"][sample]

        print(f"Processing sample: {sample}")
        print(f"  file_count: {len(jsn['files'])}")

        files = list(jsn["files"])
        module_name = select_module_name(sample)

        if INPUT_MODE == "dbs":
            ds = cmssw.Dataset(
                dataset=jsn["path"],
                lumis_per_task=1,
                file_based=True,
            )
        else:
            ds = Dataset(
                files=files,
                files_per_task=1,
                patterns=["*.root"],
            )

        command_string, display_command = build_payload_command(
            wrapper=WRAPPER,
            skim_cut=skim_cut,
            module_name=module_name,
            out_dir=".",
        )
        if command_string not in payload_command_actual_values:
            payload_command_actual_values.append(command_string)
        if display_command not in payload_command_display_values:
            payload_command_display_values.append(display_command)

        print("\nActual Lobster command string:")
        print(command_string)

        print("\nCopy-paste-safe display command:")
        print(display_command)

        skim_wf = Workflow(
            label=sample.replace("-", "_"),
            sandbox=cmssw.Sandbox(release=sandbox_location),
            dataset=ds,
            category=cat,
            extra_inputs=[WRAPPER],
            outputs=["output.root"],
            command=command_string,
            merge_command=MERGE_COMMAND,
            merge_size="537M",
            globaltag=False,
            cleanup_input=False,
        )

        workflows.append(skim_wf)

except Exception as err:
    print(
        "Error while setting up workflows. "
        "Please check your configuration and try again. Error details:"
    )
    print(err)
    raise RuntimeError("Failed to set up workflows due to the above error.") from err


# =============================================================================
# ADVANCED OPTIONS
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

selected_input_mode_counts = {INPUT_MODE: len(workflows)}
selected_storage_profile = "dbs" if storage is storage_dbs else "files"
selected_storage_mode_counts = {
    "files": 1 if storage is storage_files else 0,
    "dbs": 1 if storage is storage_dbs else 0,
}
selected_year_counts = {YEAR: len(workflows)}

print_resolved_config(
    header="resolved_config_header: Run 3 Lobster resolved configuration",
    fields=[
        ("repo_or_campaign_label", "Run 3 Lobster skimming"),
        ("testing", TESTING),
        ("input_mode", INPUT_MODE),
        ("input_mode_allowed_values_current_source", "dbs, files"),
        ("auto_mode_status", "intended_not_implemented"),
        ("auto_mode_status_not_implemented", "intended_not_implemented"),
        ("target", TARGET),
        ("type_or_campaign_type", TYPE),
        ("type", TYPE),
        ("year", YEAR),
        ("run_period", YEAR),
        ("cfg_name", CFG_NAME),
        ("derived_cfg_name_from_type_year", CFG_NAME),
        ("cfg_path", cfg_fpath),
        ("tag", TAG),
        ("workdir_path", workdir_path),
        ("plotdir_path", plotdir_path),
        ("output_path", output_path),
        ("testing_mode_path_behavior", "timestamped test paths" if TESTING else "date-stamped campaign paths"),
        ("sandbox_location", sandbox_location),
        ("wrapper", WRAPPER),
        ("master_label", master_label),
        ("workflow_count", len(workflows)),
        ("storage_profile_or_mode", selected_storage_profile),
        ("selected_year_counts", selected_year_counts),
        ("selected_input_mode_counts", selected_input_mode_counts),
        ("selected_storage_mode_counts", selected_storage_mode_counts),
        ("xrootd_servers", adv_kwargs.get("xrootd_servers", "disabled")),
        ("category_name", category_name),
        ("category_cores", category_cores),
        ("category_memory", category_memory),
        ("category_disk", category_disk),
        ("advanced_options_bad_exit_codes", adv_kwargs["bad_exit_codes"]),
        ("advanced_options_payload", adv_kwargs["payload"]),
        ("advanced_options_osg_version", adv_kwargs["osg_version"]),
        ("advanced_options_log_level", adv_kwargs["log_level"]),
        ("advanced_options_threshold_for_failure", adv_kwargs["threshold_for_failure"]),
        ("advanced_options_threshold_for_skipping", adv_kwargs["threshold_for_skipping"]),
        ("payload_command_display_or_none", format_preview_value(payload_command_display_values, limit=5)),
        ("merge_command", MERGE_COMMAND),
    ],
    sample_names=list(sorted(cfg["jsons"].keys())),
    payload_command_actual=payload_command_actual_values,
    payload_command_display=payload_command_display_values,
)


config = Config(
    label=master_label,
    workdir=workdir_path,
    plotdir=plotdir_path,
    storage=storage,
    workflows=workflows,
    advanced=AdvancedOptions(**adv_kwargs),
)
