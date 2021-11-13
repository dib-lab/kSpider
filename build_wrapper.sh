#!/usr/bin/env bash

# exit when any command fails
set -o nounset
set -o errexit
set -x
set -euox pipefail
set -e

# force clean it up
function cleanup() {
    echo "REMOVING OLD FILES IF EXISTS";
    rm -rf build/temp*
    rm -rf build/lib.linux*
    rm -rf dist/*
    rm -rf __pycache__/
    rm -rf *cxx
    rm -rf *pyc
    rm -rf src/swig_interfaces/kSpider_internal_wrap.cpp
    rm -rf pykSpider/internal/kSpider_internal.py
    rm -rf *so
    rm -rf pykSpider/kSpider.egg-info/
    rm -rf build/bdist.linux-x86_64
}

trap cleanup EXIT
cleanup

# Build the project if not already built
BUILD_DIR="build"
[[ -d ${BUILD_DIR} ]] || cmake -Bbuild && cmake --build build -j4


echo "BDIST WHEEL"
$(which python) setup.py bdist_wheel

cd dist/

$(which python) -m pip uninstall kSpider -y

$(which python) -m pip install kSpider*cp*.whl

rm -rf build/temp
rm -rf build/lib.linux*
rm -rf __pycache__/
rm -rf *cxx
rm -rf *pyc
rm -rf src/swig_interfaces/kSpider_internal_wrap.cpp
rm -rf kSpider_internal.py
rm -rf *so
rm -rf pykSpider/kSpider.egg-info/
rm -rf build/bdist.linux-x86_64