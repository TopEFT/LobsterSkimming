#!/usr/bin/env bash
set -euo pipefail

########################################
# Sparse checkout helper
# Usage:
#   git_sparse_checkout <url> <target-dir> <branch-or-tag> <path1> [path2 ...]
########################################
git_sparse_checkout() {
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

    # Initialize repo if needed
    if [[ ! -d .git ]]; then
        git init
    fi

    # Enable sparse checkout
    git config core.sparseCheckout true
    mkdir -p .git/info

    # Reset sparse-checkout file and add desired paths
    : > .git/info/sparse-checkout
    for path in "${paths[@]}"; do
        echo "Getting ${path}"
        echo "${path}" >> .git/info/sparse-checkout
    done

    # Handle remote origin
    if git remote get-url origin &>/dev/null; then
        local existing_url
        existing_url=$(git remote get-url origin)
        if [[ "${existing_url}" != "${url}" ]]; then
            echo "Remote 'origin' already exists with different URL:" >&2
            echo "  existing: ${existing_url}" >&2
            echo "  requested: ${url}" >&2
            return 1
        fi
    else
        echo "Adding remote url"
        git remote add origin "${url}"
    fi

    echo "Fetching from origin"
    git fetch --depth=1 origin "${tag}"

    echo "Checking out '${tag}'"
    # Try to create/update a local branch tracking the remote branch
    if git rev-parse --verify "${tag}" &>/dev/null; then
        git checkout "${tag}"
    else
        git checkout -b "${tag}" "origin/${tag}" 2>/dev/null || git checkout "${tag}"
    fi
}

########################################
# Install topeft cfg and json directories
########################################
install_topeft_configs() {
    local url="https://github.com/TopEFT/topeft.git"
    local tag="run3_test_mmerged"
    local topdir
    topdir=$(git rev-parse --show-toplevel)
    local dir="${topdir}/topeft"

    local cfg_dir="input_samples/cfgs"
    local json_dir="input_samples/sample_jsons"

    echo "Installing topeft cfg and json directories"
    git_sparse_checkout "${url}" "${dir}" "${tag}" "${cfg_dir}" "${json_dir}"
}

########################################
# Main
# If args are given, behave as a generic helper:
#   ./script.sh <url> <dir> <branch> <path1> [path2 ...]
# If no args, install topeft configs into current repo.
########################################
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    if [[ $# -gt 0 ]]; then
        git_sparse_checkout "$@"
    else
        install_topeft_configs
    fi
fi