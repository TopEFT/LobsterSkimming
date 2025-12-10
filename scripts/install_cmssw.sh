#!/usr/bin/env bash
# Usage: ./install_cmssw.sh <install-dir> <cmssw-version> <scram-arch>

set -euo pipefail

setup_cmssw() {
    local dir="$1"
    local cmssw_ver="$2"
    local scram_arch="$3"

    cd "${dir}" || { echo "ERROR: directory '${dir}' does not exist"; return 1; }

    if [[ -n "${CMSSW_BASE:-}" ]]; then
        echo "Already in a CMSSW release (${CMSSW_BASE})!"
        echo "Must start from a fresh environment, exiting now."
        return 1
    fi

    local cvmfs_dir=/cvmfs/cms.cern.ch

    echo "dir: ${dir}"
    echo "cmssw_release: ${cmssw_ver}"
    echo "SCRAM_ARCH: ${scram_arch}"

    if [[ -d "${cvmfs_dir}" ]]; then
        echo "Found CVMFS!"
        # Makes cmsrel available to the environment
        # shellcheck source=/dev/null
        source "${cvmfs_dir}/cmsset_default.sh"
    else
        echo "Couldn't find CVMFS directory (${cvmfs_dir}), exiting now"
        return 1
    fi

    export SCRAM_ARCH="${scram_arch}"

    echo "Setting up ${cmssw_ver}"
    scram p CMSSW "${cmssw_ver}"

    cd "${cmssw_ver}/src"
    #git clone git@github.com:sscruz/cmgtools-lite.git -b 104X_dev_nano_lepMVA CMGTools
    #git clone https://github.com/jdelrieg/topEFT_ttHMVA_Run3.git -b newcmgtools_python3 CMGTools
    git clone https://github.com/anpicci/topEFT_ttHMVA_Run3.git -b nd_run3 CMGTools
    #git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools

    echo "Getting CMS ENV from ${PWD}"
    eval "$(scramv1 runtime -sh)"
    scram b
}

setup_cmssw "$@"