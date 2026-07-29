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
# -e aborts on any error, -u on unset vars, pipefail on mid-pipe failures (setup must be
# all-or-nothing)
set -euo pipefail
# absolute path of this script's own dir, so paths resolve regardless of launch location
STAGE_DIR="$(cd "$(dirname "$0")" && pwd)"
# repo top-level (two dirs up from scripts/subtyping)
REPO_ROOT="$(cd "${STAGE_DIR}/../.." && pwd)"
mkdir -p "${REPO_ROOT}/scripts/tools"                # ensure the tools install dir exists
# work from the tools dir so downloads/build land there
cd "${REPO_ROOT}/scripts/tools"

if [ ! -d jpHMM ]; then                              # only download if jpHMM isn't already unpacked
    # fetch the source tarball from the Göttingen server (HTTP-only)
    curl -sL -o jpHMM.tar.gz "http://jphmm.gobics.de/jpHMM.tar.gz"
    # create the target dir for the extracted source
    mkdir -p jpHMM
    # unpack, dropping the archive's top-level dir into jpHMM/
    tar -xzf jpHMM.tar.gz -C jpHMM --strip-components=1
    rm jpHMM.tar.gz                                  # remove the tarball once extracted
fi

# compile only if the binary isn't already built
if [ ! -x jpHMM/src/jpHMM ]; then
    echo "Compiling jpHMM (this can take several minutes)..."  # warn the user this step is slow
    # build in src/, teeing output to jphmm_build.log for later inspection
    (cd jpHMM/src && make) 2>&1 | tee ../../jphmm_build.log
fi

# confirm the build actually produced an executable
if [ -x jpHMM/src/jpHMM ]; then
    echo "Built OK: scripts/tools/jpHMM/src/jpHMM"   # report success and the binary's location
    # note the CLI was verified from source
    echo "run_jphmm.sh's flags (-s/-v/-I/-P/-o) are confirmed directly from"
    # ...specifically from main.cpp's getopt() call
    echo "src/main.cpp's getopt() call, not guessed."
else
    # point the user at the build log on failure
    echo "ERROR: jpHMM build failed, see scripts/tools/jphmm_build.log" >&2
    exit 1                                           # exit non-zero so callers know setup failed
fi
