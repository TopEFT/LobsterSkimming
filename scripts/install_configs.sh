#!/usr/bin/env bash

set -euo pipefail

# Usage:
#   install_configs.sh <repo-url> <target-dir> <branch-or-tag> <path1> <path2> [more paths...]

# See: https://stackoverflow.com/questions/60190759/how-do-i-clone-fetch-or-sparse-checkout-a-single-directory-or-a-list-of-directo/60190760#60190760
# See Also: https://github.com/frgomes/bash-scripts/blob/master/bin/git_sparse_checkout
function git_sparse_checkout {
    if [[ $# -lt 4 ]]; then
        echo "Usage: git_sparse_checkout <url> <target-dir> <branch/tag> <path1> [path2 ...]" >&2
        return 1
    fi

    local url="$1"
    local dir="$2"
    local tag="$3"
    shift 3
    local paths=("$@")

    echo "url: ${url}"
    echo "dir: ${dir}"
    echo "tag: ${tag}"
    mkdir -p "${dir}"
    cd "${dir}" || { echo "Failed to cd to ${dir}" >&2; return 1; }
    git init
    git config core.sparseCheckout true
    mkdir -p .git/info
    : > .git/info/sparse-checkout
    for path in "${paths[@]}"; do
        echo "Getting ${path}"
        echo "${path}" >> .git/info/sparse-checkout
    done
    echo "Adding remote url"
    git remote add origin "${url}"
    echo "Fetching"
    git fetch origin
    echo "Checking out tag"
    git checkout "${tag}"
}

# Installs the topeft input sample cfg and JSON directories
function install_topeft_configs {
    local url=https://github.com/TopEFT/topeft.git
    local topeft_branch=run3_test_mmerged
    local prj_head
    prj_head="$(git rev-parse --show-toplevel)/topeft"
    # Sparse paths are relative to prj_head.
    local cfg_dir=input_samples/cfgs
    local sample_jsons_dir=input_samples/sample_jsons

    git_sparse_checkout "${url}" "${prj_head}" "${topeft_branch}" "${cfg_dir}" "${sample_jsons_dir}"
}

if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <repo-url> <target-dir> <branch-or-tag> <path1> <path2> [more paths...]" >&2
    exit 1
fi

url="$1"
dir="$2"
tag="$3"
shift 3

git_sparse_checkout "${url}" "${dir}" "${tag}" "$@"
