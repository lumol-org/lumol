#!/bin/bash -ex
# Install tools for Travis CI build

# Build and cache kcov
if [ -f $HOME/local/bin/kcov ]; then
    echo "Using cached kcov from ~/local/bin/kcov"
else
    if [[ "${TRAVIS_OS_NAME}" == "osx" ]]; then
        export OPENSSL_ROOT_DIR=$(brew --prefix openssl)
    fi

    mkdir -p $HOME/local
    cd $HOME/local
    wget https://github.com/SimonKagstrom/kcov/archive/master.tar.gz
    tar xzf master.tar.gz && mkdir kcov-master/build && cd kcov-master/build
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/local ..
    make install
    cd $TRAVIS_BUILD_DIR
fi
