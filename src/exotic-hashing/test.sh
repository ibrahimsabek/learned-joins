#!/bin/bash

# setup script
set -e
cd "$(dirname "$0")"

# build and run tests
./build.sh eh_tests RELEASE
cmake-build-release/src/eh_tests $@
