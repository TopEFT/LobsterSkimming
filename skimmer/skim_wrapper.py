import argparse
import os
import shlex
import subprocess


parser = argparse.ArgumentParser(description="")
parser.add_argument("infiles", nargs="+", help="")
parser.add_argument("--cut", "-c", type=str, required=True, help="")
parser.add_argument("--module", "-m", type=str, required=True, help="")
parser.add_argument("--out-dir", "-o", type=str, default=".", help="")
parser.add_argument("--nevents", "-n", type=str, default=None, help="")

args = parser.parse_args()

skim_cut = args.cut

if (
    skim_cut
    and len(skim_cut) >= 2
    and skim_cut[0] == skim_cut[-1]
    and skim_cut[0] in ("'", '"')
):
    print("WARNING: stripping one accidental layer of literal outer quotes from --cut")
    skim_cut = skim_cut[1:-1]

module = args.module
out_dir = args.out_dir
infiles = args.infiles
nevents = args.nevents

print("Full list of input arguments:", args)
print("Cut argument:", args.cut)
print("Cut argument repr:", repr(args.cut))
print("Effective cut repr:", repr(skim_cut))


indent = " " * 4


def list_cwd():
    lines = ["Current working directory:"]
    for f_name in os.listdir("."):
        lines.append(indent + f_name)
    print("\n".join(lines))


list_cwd()

local_files = []
cwd_files = set(os.listdir("."))

for inf in infiles:
    local_name = inf.rsplit("/")[-1]
    local_files.append(local_name)

    if local_name.replace("file:", "") in cwd_files or local_name in cwd_files:
        continue

    print("inf, local_name:", inf, local_name)

    cmd_args = ["xrdcp", "-f", inf, local_name]
    print("Copy command:", shlex.join(cmd_args))
    subprocess.check_call(cmd_args)

    cwd_files.add(local_name)

to_skim = local_files

cmd_args = ["nano_postproc.py"]
cmd_args.extend(["-c", skim_cut.replace(" ", "")])
cmd_args.extend(["-I", "CMGTools.NanoProc.tools.nanoAOD.lepMVA_run3", module])
cmd_args.extend([out_dir])
cmd_args.extend(to_skim)

if nevents:
    cmd_args.extend(["-N", nevents])

print("Skim command:", shlex.join(cmd_args))
subprocess.check_call(cmd_args)

list_cwd()

to_merge = [x.rsplit("/")[-1].replace(".root", "_Skim.root") for x in to_skim]

cmd_args = ["haddnano.py", "output.root"]
cmd_args.extend(to_merge)

print("Merge command:", shlex.join(cmd_args))
subprocess.check_call(cmd_args)