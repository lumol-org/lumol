#!/bin/bash -ex
# Setup environement variables for Travis CI build

export PATH=~/.local/bin:$PATH:~/local/bin:~/.cargo/bin:~/Library/Python/2.7/bin
export RUSTFLAGS="-C link-dead-code"

if [[ "$TRAVIS_RUST_VERSION" == "stable" && "$TRAVIS_OS_NAME" == "linux" ]]; then
    export DO_COVERAGE=true
else
    export DO_COVERAGE=false
fi

if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    export CC=gcc-4.9
    export CXX=g++-4.9
fi
