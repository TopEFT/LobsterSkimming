import os
import subprocess



def main():
    top_dir,fname = os.path.split(__file__)
    if not top_dir:
        top_dir = "."

    os.chdir(top_dir)
    abs_path = subprocess.check_output(["git","rev-parse","--show-toplevel"])
    abs_path = abs_path.strip()

    cmssw_release = "CMSSW_10_6_19_patch2"
    scram_arch = "slc7_amd64_gcc700"


    if os.path.exists("topeft"):
        print("topeft directory already installed, skipping this part\n")
    else:
        print("Installing topeft input sample cfg and JSON directories")
        topeft_url = "https://github.com/TopEFT/topeft.git"
        topeft_branch = "run3_test_mmerged"
        prj_head = "{}/topeft".format(abs_path)
        # Sparse paths are relative to prj_head. These install under
        # topeft/input_samples/cfgs and topeft/input_samples/sample_jsons.
        cfg_dir = "input_samples/cfgs"
        sample_jsons_dir = "input_samples/sample_jsons"
        subprocess.check_call(["./scripts/install_configs.sh",topeft_url,prj_head,topeft_branch,cfg_dir,sample_jsons_dir])
        print("")

    if os.path.exists(cmssw_release):
        print("CMSSW release {} detected, skipping this part".format(cmssw_release))
    else:
        print("Setting up CMSSW release and getting NanoAODTools")
        subprocess.check_call(["./scripts/install_cmssw.sh",abs_path,cmssw_release,scram_arch])

    print("\nDone!\nMake sure to do a cmsenv before activating and/or using lobster!")
main()
