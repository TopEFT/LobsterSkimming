#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   install_configs.sh <repo-url> <target-dir> <branch-or-tag> <path1> <path2> [more paths...]

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
    # Try local branch, then remote branch, then tag
    if git show-ref --verify --quiet "refs/heads/${tag}"; then
        git checkout "${tag}"
    elif git show-ref --verify --quiet "refs/remotes/origin/${tag}"; then
        git checkout -b "${tag}" "origin/${tag}"
    else
        git checkout "${tag}"
    fi
}

if [[ "$#" -lt 5 ]]; then
    echo "Usage: $0 <repo-url> <target-dir> <branch-or-tag> <path1> <path2> [more paths...]" >&2
    exit 1
fi

url="$1"
dir="$2"
tag="$3"
shift 3

git_sparse_checkout "${url}" "${dir}" "${tag}" "$@"
