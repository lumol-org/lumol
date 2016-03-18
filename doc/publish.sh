#!/bin/bash

# Build the docs and the user manual, and push them to github pages.

set -e
if test "${TRAVIS_RUST_VERSION}" != "stable" &&
        "${TRAVIS_PULL_REQUEST}" != "false"  &&
        "${TRAVIS_BRANCH}" != "master"
then
    exit
fi

cargo doc
cargo install mdbook
mdbook build doc
cp -r doc/book target/doc
cp doc/index.html target/doc

cd target
mkdir website
mv doc website/latest

cd website
# Redirect to latest docs
cat <<EOF > index.html
<meta http-equiv=refresh content=0;url=latest/index.html>
EOF

git init
git config user.name "Travis CI"
git config user.email "luthaf@luthaf.fr"
git add .
git commit -m "Update documentation"

git push --force --quiet "https://${GH_TOKEN}@github.com/Luthaf/cymbalum.git" master:gh-pages > /dev/null 2>&1
