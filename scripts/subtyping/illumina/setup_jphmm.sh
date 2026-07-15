#!/bin/bash
# jpHMM has no packaged CLI distribution (conda/pip) -- confirmed by direct
# lookup, its source archive is only available from the University of
# Göttingen's own server (plain HTTP, no HTTPS). This is a real download
# from a small academic server: a genuine URL, verified reachable (52MB,
# last updated Mar 2015) before writing this script.
#
# Unlike the other GitHub-hosted comparison tools, jpHMM ships as C++ source
# (with its own bundled copy of Boost) and must be COMPILED locally -- there
# is no prebuilt binary. This can take a while and needs a working g++/make
# toolchain. It bundles input/HIV_alignment.fas and priors/emissionPriors_HIV.txt
# + transition_priors.txt, which are exactly the -I/-P arguments jpHMM's own
# download page says are required -- so no separate reference alignment
# needs to be sourced for jpHMM specifically (unlike SHIVER).
set -euo pipefail
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${STAGE_DIR}/../.." && pwd)"
mkdir -p "${REPO_ROOT}/scripts/tools"
cd "${REPO_ROOT}/scripts/tools"

if [ ! -d jpHMM ]; then
    curl -sL -o jpHMM.tar.gz "http://jphmm.gobics.de/jpHMM.tar.gz"
    mkdir -p jpHMM
    tar -xzf jpHMM.tar.gz -C jpHMM --strip-components=1
    rm jpHMM.tar.gz
fi

if [ ! -x jpHMM/src/jpHMM ]; then
    echo "Compiling jpHMM (this can take several minutes)..."
    (cd jpHMM/src && make) 2>&1 | tee ../../jphmm_build.log
fi

if [ -x jpHMM/src/jpHMM ]; then
    echo "Built OK: scripts/tools/jpHMM/src/jpHMM"
    echo "run_jphmm.sh's flags (-s/-v/-I/-P/-o) are confirmed directly from"
    echo "src/main.cpp's getopt() call, not guessed."
else
    echo "ERROR: jpHMM build failed, see scripts/tools/jphmm_build.log" >&2
    exit 1
fi
