#!/bin/bash -ex
# Build the docs and the user manual

# Build the doc
cargo doc --all --no-deps
mdbook build doc

# Move it to the right place
rm -rf target/gh-pages/latest
mkdir -p target/gh-pages/latest
cd target/gh-pages/
cp -r ../../doc/book latest/book
cp ../../doc/index.html latest/index.html
mv ../doc/* latest/

cat <<EOF > index.html
<meta http-equiv=refresh content=0;url=latest/>
EOF
