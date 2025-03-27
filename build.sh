#! /bin/bash
set -e

# check dependencies
if ! [ -x "$(command -v uv)" ]; then
    echo 'Error: uv is not installed.' >&2
    exit 1
fi

pushd $(dirname $0) > /dev/null
    proot=$PWD
    uv sync
    source .venv/bin/activate
    pushd src/homological_threading/fortran > /dev/null
        make clean
        make
    popd
    export PYTHONPATH=$PYTHONPATH:$proot/src
popd

echo "Build complete."
