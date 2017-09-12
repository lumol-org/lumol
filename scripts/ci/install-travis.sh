#!/bin/bash -ex
# Install tools for Travis CI build

# Build and cache kcov
if [ -f $HOME/local/bin/kcov ]; then
    echo "Using cached kcov from ~/local/bin/kcov"
else
    mkdir -p $HOME/local
    cd $HOME/local
    wget https://github.com/SimonKagstrom/kcov/archive/master.tar.gz
    tar xzf master.tar.gz && mkdir kcov-master/build && cd kcov-master/build
    cmake -DCMAKE_INSTALL_PREFIX=$HOME/local ..
    make install
    cd $TRAVIS_BUILD_DIR
fi

# Build and cache mdbook
if [ -f $HOME/.cargo/bin/mdbook ]; then
    echo "Using cached mdbook from ~/.cargo/bin/mdbook"
else
    cargo install mdbook --vers "0.0.25"
fi

cd $TRAVIS_BUILD_DIR

pip install --user -r scripts/ci/requirements.txt
